#!/usr/bin/env python3
"""Screen Co-4-only pan-gene candidates using existing TAIR10 BLASTP hits."""

import csv
import json
from pathlib import Path
from collections import Counter

ROOT = Path("/data/users/gercan/Gene_annotation")
GS = ROOT / "results/comparative_genomics/GENESPACE/run_5_with_TAIR10"
OUT = GS / "process_pangenome/Co4_specific_functional_analysis"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES = GS / "process_pangenome/tables/Co4_only_from_Co_4_reference.tsv"
ANNOTATIONS = OUT / "Co4_specific_functional_annotations.tsv"
PROTEINS = ROOT / "results/merged/Co4_maker.all.maker.proteins.fasta"
BLAST = ROOT / "results/functional_annotation/BLASTP/blastp_TAIR10.besthits.tsv"

def load_fasta_lengths(path):
    lengths = {}
    current = None

    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())

    return lengths

with CANDIDATES.open() as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

with ANNOTATIONS.open() as f:
    annotations = {
        r["gene_id"]: r
        for r in csv.DictReader(f, delimiter="\t")
    }

lengths = load_fasta_lengths(PROTEINS)
hits = {}

with BLAST.open() as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        try:
            record = {
                "query": p[0],
                "subject": p[1],
                "identity": float(p[2]),
                "alignment_length": int(p[3]),
                "qstart": int(p[6]),
                "qend": int(p[7]),
                "evalue": float(p[10]),
                "bitscore": float(p[11]),
            }
        except ValueError:
            continue

        if (
            p[0] not in hits
            or record["bitscore"] > hits[p[0]]["bitscore"]
        ):
            hits[p[0]] = record

rows = []

for candidate in candidates:
    gene = candidate["gene_id"]
    annotation = annotations.get(gene, {})
    tids = [
        t for t in annotation.get("transcript_ids", "").split(";")
        if t
    ]

    best = None

    for tid in tids:
        protein_ids = [tid, tid + "-p1", tid + "-protein-1"]
        protein_length = next(
            (lengths[p] for p in protein_ids if p in lengths),
            None
        )

        hit = hits.get(tid)

        if not hit:
            continue

        query_span = abs(hit["qend"] - hit["qstart"]) + 1
        qcov = (
            100.0 * query_span / protein_length
            if protein_length
            else None
        )

        entry = {
            **hit,
            "transcript_id": tid,
            "protein_length": protein_length,
            "query_coverage_pct": qcov,
        }

        if best is None or entry["bitscore"] > best["bitscore"]:
            best = entry

    if best is None:
        classification = "no_TAIR10_besthit"
    elif best["query_coverage_pct"] is None:
        classification = "protein_length_unresolved"
    elif (
        best["identity"] >= 95
        and best["query_coverage_pct"] >= 80
    ):
        classification = "strong_reference_similarity"
    elif (
        best["identity"] >= 70
        and best["query_coverage_pct"] >= 70
    ):
        classification = "moderate_reference_similarity"
    else:
        classification = "weak_or_partial_reference_similarity"

    row = {
        "pgID": candidate["pgID"],
        "gene_id": gene,
        "transcript_id": best["transcript_id"] if best else "",
        "protein_length": best["protein_length"] if best else "",
        "TAIR10_hit": best["subject"] if best else "",
        "identity_pct": best["identity"] if best else "",
        "query_coverage_pct": (
            round(best["query_coverage_pct"], 2)
            if best and best["query_coverage_pct"] is not None
            else ""
        ),
        "evalue": best["evalue"] if best else "",
        "bitscore": best["bitscore"] if best else "",
        "best_AED": annotation.get("best_AED", ""),
        "keyword_categories": annotation.get("keyword_categories", ""),
        "screening_class": classification,
    }
    rows.append(row)

columns = list(rows[0].keys())
target = OUT / "Co4_candidate_TAIR10_protein_screen.tsv"

with target.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

counts = Counter(r["screening_class"] for r in rows)

summary = {
    "total_candidate_genes": len(rows),
    "protein_fasta_records": len(lengths),
    "screening_classes": dict(counts),
    "note": (
        "Query coverage is estimated from the best BLASTP alignment "
        "query coordinates. This is not reference coverage. "
        "No hit in besthits.tsv is not evidence of true gene absence."
    )
}

(OUT / "TAIR10_protein_screen_summary.json").write_text(
    json.dumps(summary, indent=2)
)

print(json.dumps(summary, indent=2))
print("Saved:", target)

print("\n=== High-identity examples ===")
example_ids = [
    "augustus_masked-contig_132-processed-gene-3.30",
    "augustus_masked-contig_133-processed-gene-36.25",
    "maker-contig_148-augustus-gene-14.29",
]

for row in rows:
    if row["gene_id"] in example_ids:
        print(row)
