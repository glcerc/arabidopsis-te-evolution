#!/usr/bin/env python3
"""Prepare 119 remaining candidates for sensitive genomic screening."""

import csv
import json
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
PAV = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/"
    "Co4_specific_functional_analysis/PAV_validation"
)
PAV.mkdir(parents=True, exist_ok=True)

def read_ids(path):
    with path.open() as handle:
        return {
            row["gene_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }

all_candidates = read_ids(
    PAV / "Co4_126_no_alignment_candidates.tsv"
)
validated_seven = read_ids(
    PAV / "Co4_higher_priority_PAV_candidates.tsv"
)

remaining = all_candidates - validated_seven

if len(all_candidates) != 126 or len(validated_seven) != 7:
    raise RuntimeError("Unexpected input candidate counts")

if len(remaining) != 119:
    raise RuntimeError(f"Expected 119 candidates, found {len(remaining)}")

bed_in = PAV / "Co4_585_candidate_genes.bed"
bed_out = PAV / "Co4_119_remaining_candidates.bed"

found = set()

with bed_in.open() as src, bed_out.open("w") as dst:
    for line in src:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 4 and fields[3] in remaining:
            dst.write(line)
            found.add(fields[3])

if found != remaining:
    raise RuntimeError(
        f"Missing BED records: {sorted(remaining - found)[:10]}"
    )

summary = {
    "original_no_alignment_candidates": len(all_candidates),
    "previously_screened_by_blastn": len(validated_seven),
    "remaining_candidates": len(remaining),
    "bed_file": str(bed_out)
}

(PAV / "Co4_119_preparation.json").write_text(
    json.dumps(summary, indent=2)
)

print(json.dumps(summary, indent=2))
