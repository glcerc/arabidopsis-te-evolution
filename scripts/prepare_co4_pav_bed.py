#!/usr/bin/env python3
"""Prepare genomic intervals for candidate Co-4-specific genes."""

import csv
import json
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
BASE = ROOT / "results/comparative_genomics/GENESPACE/run_5_with_TAIR10"
OUT = BASE / "process_pangenome/Co4_specific_functional_analysis/PAV_validation"
OUT.mkdir(parents=True, exist_ok=True)

candidates_path = (
    BASE / "process_pangenome/tables/Co4_only_from_Co_4_reference.tsv"
)
gff_path = (
    ROOT / "results/functional_annotation/integrated/"
    "Co4_final_functional_annotation.gff3"
)

with candidates_path.open() as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

genes = {row["gene_id"] for row in candidates}
found = {}

with gff_path.open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9 or p[2] != "gene":
            continue

        attrs = dict(
            field.split("=", 1)
            for field in p[8].split(";")
            if "=" in field
        )

        gene_id = attrs.get("ID")
        if gene_id in genes:
            found[gene_id] = (
                p[0], int(p[3]) - 1, int(p[4])
            )

missing = genes - set(found)
if missing:
    raise RuntimeError(f"Missing gene models: {sorted(missing)[:10]}")

bed = OUT / "Co4_585_candidate_genes.bed"

with bed.open("w") as f:
    for gene_id in sorted(genes):
        chrom, start, end = found[gene_id]
        f.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

summary = {
    "candidate_genes": len(genes),
    "matched_gene_intervals": len(found),
    "bed_file": str(bed),
    "coordinate_system": "BED: 0-based, half-open",
    "sequence_type": "Genomic gene intervals including introns",
}

(OUT / "preparation.json").write_text(
    json.dumps(summary, indent=2)
)

print(json.dumps(summary, indent=2))
