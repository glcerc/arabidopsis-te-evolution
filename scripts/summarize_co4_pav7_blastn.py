#!/usr/bin/env python3
"""Summarize sensitive BLASTN screening for seven Co-4 candidates."""

import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
PAV = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/"
    "Co4_specific_functional_analysis/PAV_validation"
)
ACCS = ["Col-0", "Altai-5", "Are-6", "Etna-2"]

def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

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

candidates = read_tsv(
    PAV / "Co4_higher_priority_PAV_candidates.tsv"
)

records = []
counts = {}

for acc in ACCS:
    hits = defaultdict(list)
    path = PAV / "validation_7" / f"{acc}_blastn.tsv"

    with path.open() as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            hits[p[0]].append({
                "qlen": int(p[1]),
                "target": p[2],
                "identity": float(p[4]),
                "length": int(p[5]),
                "qstart": int(p[6]),
                "qend": int(p[7]),
                "sstart": int(p[8]),
                "send": int(p[9]),
                "evalue": float(p[10]),
                "bitscore": float(p[11]),
            })

    counts[acc] = {}

    for candidate in candidates:
        gene = candidate["gene_id"]
        gene_hits = hits.get(gene, [])
        if gene_hits:
            best = max(gene_hits, key=lambda h: h["bitscore"])
            qlen = best["qlen"]
            qcov = 100 * (
                abs(best["qend"] - best["qstart"]) + 1
            ) / qlen

            # Only combine HSPs on the same target sequence.
            # This is exploratory, not locus-resolved coverage.
            same_target = [
                h for h in gene_hits
                if h["target"] == best["target"]
            ]
            union_cov = 100 * union_length([
                (min(h["qstart"], h["qend"]) - 1,
                 max(h["qstart"], h["qend"]))
                for h in same_target
            ]) / qlen

            identity = best["identity"]
            target = best["target"]
            bitscore = best["bitscore"]
        else:
            qcov = 0
            union_cov = 0
            identity = 0
            target = ""
            bitscore = 0

        if not gene_hits:
            status = "no_blastn_hit"
        elif qcov >= 90 and identity >= 90:
            status = "strong_single_alignment"
        elif qcov >= 50 and identity >= 80:
            status = "moderate_single_alignment"
        else:
            status = "weak_or_partial_alignment"

        counts[acc][status] = counts[acc].get(status, 0) + 1

        records.append({
            "gene_id": gene,
            "pgID": candidate["pgID"],
            "accession": acc,
            "hit_count": len(gene_hits),
            "best_target": target,
            "best_identity_pct": round(identity, 2),
            "best_query_coverage_pct": round(qcov, 2),
            "same_target_union_coverage_pct": round(
                union_cov, 2
            ),
            "best_bitscore": bitscore,
            "status": status,
        })

target = PAV / "Co4_7_Blastn_validation_summary.tsv"

with target.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=list(records[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(records)

report = {
    "candidates": len(candidates),
    "reference_genomes": ACCS,
    "per_accession_classes": counts,
    "limitations": [
        "Best HSP coverage is not cumulative locus coverage.",
        "Same-target union may combine distinct genomic loci.",
        "No BLASTN hit does not prove biological absence.",
        "Synteny and assembly quality must be checked."
    ]
}

(PAV / "Co4_7_Blastn_validation_summary.json").write_text(
    json.dumps(report, indent=2)
)

print("=== BLASTN classifications ===")
for acc, values in counts.items():
    print(acc, values)

print("\n=== Per-gene comparisons ===")
for candidate in candidates:
    gene = candidate["gene_id"]
    print("\n", gene)
    for row in records:
        if row["gene_id"] == gene:
            print(
                row["accession"],
                "identity:", row["best_identity_pct"],
                "coverage:", row["best_query_coverage_pct"],
                "status:", row["status"]
            )

print("\nSaved:", target)
