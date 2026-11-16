#!/usr/bin/env python3
"""Prioritize candidate Co-4-specific genes after genomic screening."""

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
OUT = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/"
    "Co4_specific_functional_analysis"
)
PAV = OUT / "PAV_validation"

def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

pav = read_tsv(PAV / "Co4_PAV_combined.tsv")
functional = {
    r["gene_id"]: r
    for r in read_tsv(
        OUT / "Co4_specific_functional_annotations.tsv"
    )
}
protein = {
    r["gene_id"]: r
    for r in read_tsv(
        OUT / "Co4_candidate_TAIR10_protein_screen.tsv"
    )
}

selected = []
counts = Counter()

for row in pav:
    if row["PAV_screen_class"] != "no_alignment_in_four_references":
        continue

    gene = row["gene_id"]
    func = functional[gene]
    prot = protein[gene]

    aed = float(func["best_AED"]) if func["best_AED"] else None
    category = func["keyword_categories"]
    gene_length = int(func["gene_length_bp"])

    low_aed = aed is not None and aed <= 0.25
    te_related = "transposon_related" in category.split(";")
    long_enough = gene_length >= 300

    # This is a prioritization filter, not proof of gene absence.
    priority = (
        "higher_priority"
        if low_aed and not te_related and long_enough
        else "requires_review"
    )

    result = {
        "pgID": row["pgID"],
        "gene_id": gene,
        "priority": priority,
        "gene_length_bp": gene_length,
        "best_AED": func["best_AED"],
        "protein_screen": prot["screening_class"],
        "TAIR10_protein_hit": prot["TAIR10_hit"],
        "Pfam": func["Pfam"],
        "InterPro": func["InterPro"],
        "functional_categories": category,
        "no_minimap2_alignment_in_four_genomes": "yes",
    }

    selected.append(result)
    counts[priority] += 1
    counts["AED_le_0.25"] += int(low_aed)
    counts["transposon_related"] += int(te_related)
    counts["length_ge_300bp"] += int(long_enough)
    counts["no_TAIR10_protein_besthit"] += int(
        prot["screening_class"] == "no_TAIR10_besthit"
    )

columns = list(selected[0]) if selected else [
    "pgID", "gene_id", "priority"
]

write_tsv(
    PAV / "Co4_126_no_alignment_candidates.tsv",
    selected, columns
)
write_tsv(
    PAV / "Co4_higher_priority_PAV_candidates.tsv",
    [r for r in selected if r["priority"] == "higher_priority"],
    columns
)

report = {
    "no_alignment_gene_records": len(selected),
    "no_alignment_pangene_groups": len(
        {r["pgID"] for r in selected}
    ),
    "priority_statistics": dict(counts),
    "warning": (
        "No minimap2 asm5 alignment is not proof of absence. "
        "Candidates require sensitive sequence searches, "
        "synteny checks, and TE filtering."
    )
}

(PAV / "Co4_PAV_candidate_prioritization.json").write_text(
    json.dumps(report, indent=2)
)

print(json.dumps(report, indent=2))
print("Saved candidate tables to:", PAV)
