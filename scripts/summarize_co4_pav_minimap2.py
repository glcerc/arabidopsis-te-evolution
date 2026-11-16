#!/usr/bin/env python3
"""Summarize Co-4 candidate-gene minimap2 PAF alignments."""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
BASE = ROOT / "results/comparative_genomics/GENESPACE/run_5_with_TAIR10"
OUT = BASE / "process_pangenome/Co4_specific_functional_analysis"
PAV = OUT / "PAV_validation"
ACCESSIONS = ["Col-0", "Altai-5", "Are-6", "Etna-2"]

candidates = list(csv.DictReader(
    (BASE / "process_pangenome/tables/Co4_only_from_Co_4_reference.tsv").open(),
    delimiter="\t"
))
annotation = {
    r["gene_id"]: r for r in csv.DictReader(
        (OUT / "Co4_specific_functional_annotations.tsv").open(),
        delimiter="\t"
    )
}

def union_length(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    start, end = intervals[0]
    total = 0
    for a, b in intervals[1:]:
        if a <= end:
            end = max(end, b)
        else:
            total += end - start
            start, end = a, b
    return total + end - start

def classify(coverage, identity, found):
    if not found:
        return "no_detected_alignment"
    if coverage >= 90 and identity >= 90:
        return "strong_genomic_match"
    return "partial_or_divergent_match"

all_hits = {}
per_reference = []
summary = {}

for accession in ACCESSIONS:
    path = PAV / "minimap2_results" / f"Co4_vs_{accession}.paf"
    if not path.exists():
        raise FileNotFoundError(path)

    hits = defaultdict(list)

    with path.open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                qname, qlen = p[0], int(p[1])
                qs, qe = int(p[2]), int(p[3])
                target = p[5]
                ts, te = int(p[7]), int(p[8])
                matches, block = int(p[9]), int(p[10])
                mapq = int(p[11])
            except ValueError:
                continue

            if block <= 0 or qlen <= 0:
                continue

            hits[qname].append({
                "qlen": qlen,
                "qs": qs,
                "qe": qe,
                "target": target,
                "ts": ts,
                "te": te,
                "identity": 100 * matches / block,
                "mapq": mapq,
                "strand": p[4],
            })

    counts = Counter()

    for candidate in candidates:
        gene = candidate["gene_id"]
        alignments = hits.get(gene, [])

        if alignments:
            # Primary screening: rank individual alignments by
            # covered query length and estimated identity.
            best = max(
                alignments,
                key=lambda h: (
                    (h["qe"] - h["qs"]) * h["identity"],
                    h["mapq"]
                )
            )

            qlen = best["qlen"]
            coverage = 100 * (best["qe"] - best["qs"]) / qlen
            identity = best["identity"]

            # Diagnostic only: same chromosome can contain multiple loci.
            same_target = [
                h for h in alignments
                if h["target"] == best["target"]
            ]
            target_union = 100 * union_length([
                (h["qs"], h["qe"]) for h in same_target
            ]) / qlen

            status = classify(coverage, identity, True)
            target = best["target"]
            mapq = best["mapq"]
        else:
            coverage = 0
            identity = 0
            target_union = 0
            status = classify(0, 0, False)
            target = ""
            mapq = ""

        row = {
            "pgID": candidate["pgID"],
            "gene_id": gene,
            "accession": accession,
            "alignment_count": len(alignments),
            "best_target": target,
            "best_query_coverage_pct": round(coverage, 2),
            "best_identity_pct": round(identity, 2),
            "same_target_union_coverage_pct": round(
                target_union, 2
            ),
            "best_mapq": mapq,
            "status": status
        }

        per_reference.append(row)
        all_hits[(gene, accession)] = row
        counts[status] += 1

    summary[accession] = dict(counts)

columns = list(per_reference[0])
with (PAV / "Co4_PAV_per_accession.tsv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(per_reference)

combined = []
overall = Counter()

for candidate in candidates:
    gene = candidate["gene_id"]
    statuses = [
        all_hits[(gene, accession)]["status"]
        for accession in ACCESSIONS
    ]

    if "strong_genomic_match" in statuses:
        group = "strong_match_in_at_least_one_reference"
    elif "partial_or_divergent_match" in statuses:
        group = "ambiguous_genomic_presence"
    else:
        group = "no_alignment_in_four_references"

    a = annotation.get(gene, {})

    row = {
        "pgID": candidate["pgID"],
        "gene_id": gene,
        "PAV_screen_class": group,
        "best_AED": a.get("best_AED", ""),
        "keyword_categories": a.get("keyword_categories", ""),
        "UniProt": a.get("UniProt", ""),
        "TAIR10": a.get("TAIR10", ""),
        "Pfam": a.get("Pfam", ""),
    }

    for accession in ACCESSIONS:
        hit = all_hits[(gene, accession)]
        row[f"{accession}_status"] = hit["status"]
        row[f"{accession}_coverage"] = hit[
            "best_query_coverage_pct"
        ]
        row[f"{accession}_identity"] = hit[
            "best_identity_pct"
        ]

    combined.append(row)
    overall[group] += 1

with (PAV / "Co4_PAV_combined.tsv").open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(combined[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(combined)

report = {
    "total_candidate_genes": len(candidates),
    "per_accession": summary,
    "combined": dict(overall),
    "screening_thresholds": {
        "strong_query_coverage_pct": 90,
        "strong_identity_pct": 90
    },
    "limitations": [
        "A strong genomic DNA match does not establish a functional gene.",
        "No minimap2 alignment does not establish biological absence.",
        "The best single alignment is used for screening classification.",
        "The same-target union is diagnostic, not locus-resolved coverage.",
        "Further sensitive search, synteny and assembly validation are needed."
    ]
}

(PAV / "Co4_PAV_summary.json").write_text(
    json.dumps(report, indent=2)
)

print("=== Per-accession results ===")
for name, counts in summary.items():
    print(name, counts)

print("\n=== Combined screening ===")
for name, count in overall.items():
    print(name, count)

print("\nSaved results to:", PAV)
