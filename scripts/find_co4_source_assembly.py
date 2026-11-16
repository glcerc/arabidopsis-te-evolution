#!/usr/bin/env python3
"""Locate the original Co-4 Flye assembly used for annotation."""

import os
from pathlib import Path

ROOTS = [
    Path("/data/users/gercan/assembly_annotation_course"),
    Path("/data/users/gercan/Gene_annotation"),
    Path("/data/users/gercan/TE_annotation"),
]

OUT = Path(
    "/data/users/gercan/Gene_annotation/"
    "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/"
    "Co4_specific_functional_analysis"
)
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_CONTIGS = 139
EXPECTED_BP = 136966697

suffixes = (".fa", ".fasta", ".fna")
candidates = []

for root in ROOTS:
    if not root.exists():
        continue

    for directory, dirs, files in os.walk(root):
        dirs[:] = [
            d for d in dirs
            if d not in {
                "github_submission", "busco_downloads",
                "run_brassicales_odb10", "maker_datastore",
                ".git", "conda", "envs"
            }
            and not d.endswith("_datastore")
        ]

        for name in files:
            if not name.lower().endswith(suffixes):
                continue

            path = Path(directory) / name
            try:
                size = path.stat().st_size
            except OSError:
                continue

            if not 90_000_000 <= size <= 180_000_000:
                continue

            candidates.append(path)

matches = []

for path in sorted(set(candidates)):
    try:
        n = 0
        total = 0
        has_contig_120 = False

        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    n += 1
                    name = line[1:].split()[0]
                    if name == "contig_120":
                        has_contig_120 = True
                else:
                    total += len(line.strip())

        print(
            f"{n} contigs | {total} bp | "
            f"contig_120={has_contig_120} | {path}"
        )

        if (
            n == EXPECTED_CONTIGS
            and total == EXPECTED_BP
            and has_contig_120
        ):
            matches.append(str(path))

    except (OSError, UnicodeError) as exc:
        print("SKIPPED:", path, exc)

target = OUT / "source_assembly_candidates.txt"
target.write_text(
    "\n".join(matches) + ("\n" if matches else "")
)

print("\n=== MATCHING ASSEMBLIES ===")
for match in matches:
    print(match)

print("Number of matches:", len(matches))
print("Saved:", target)
