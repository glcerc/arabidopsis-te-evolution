#!/usr/bin/env python3

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
SOURCE = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/pangenes/TAIR10_pangenes.txt.gz"
)
OUT = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/tables"
)
OUT.mkdir(parents=True, exist_ok=True)

groups = defaultdict(set)
co4_records = defaultdict(list)
all_genomes = set()
flags = defaultdict(int)

with gzip.open(SOURCE, "rt") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    print("Columns:", reader.fieldnames)

    for row in reader:
        pgid = row["pgID"]
        genome = row["genome"]
        all_genomes.add(genome)
        groups[pgid].add(genome)

        if genome == "Co_4":
            co4_records[pgid].append(row)
            flags[row["flag"]] += 1

specific = {
    pgid for pgid, genomes in groups.items()
    if genomes == {"Co_4"}
}

output_rows = []

for pgid in specific:
    for row in co4_records[pgid]:
        output_rows.append({
            "pgID": pgid,
            "gene_id": row["id"],
            "flag": row["flag"],
            "chr": row["chr"],
            "start": row["start"],
            "end": row["end"]
        })

output_file = OUT / "Co4_specific_pangene_candidates.tsv"

with output_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "pgID", "gene_id", "flag",
            "chr", "start", "end"
        ],
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(sorted(
        output_rows,
        key=lambda x: (int(x["pgID"]), x["gene_id"])
    ))

summary = {
    "source": str(SOURCE),
    "genomes_observed": sorted(all_genomes),
    "total_pangene_groups": len(groups),
    "Co4_only_pangene_groups": len(specific),
    "Co4_candidate_gene_records": len(output_rows),
    "Co4_flags_all_records": dict(flags),
    "limitation": (
        "These are candidate Co-4-specific pan-gene groups "
        "within the evaluated GENESPACE dataset. "
        "Absence of detected orthology does not establish "
        "biological absence in other accessions."
    )
}

(OUT / "Co4_specific_pangene_summary.json").write_text(
    json.dumps(summary, indent=2)
)

print(json.dumps(summary, indent=2))
print("Saved:", output_file)
