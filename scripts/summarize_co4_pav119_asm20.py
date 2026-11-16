#!/usr/bin/env python3
"""Summarize asm20 screening of 119 Co-4 candidate genes."""

import csv
import json
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path("/data/users/gercan/Gene_annotation")
PAV = (
    ROOT / "results/comparative_genomics/GENESPACE/"
    "run_5_with_TAIR10/process_pangenome/"
    "Co4_specific_functional_analysis/PAV_validation"
)
ACCS = ["Col-0", "Altai-5", "Are-6", "Etna-2"]

def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

candidates = read_tsv(
    PAV / "Co4_126_no_alignment_candidates.tsv"
)
screened7 = {
    r["gene_id"]
    for r in read_tsv(
        PAV / "Co4_higher_priority_PAV_candidates.tsv"
    )
}
remaining = [
    r for r in candidates if r["gene_id"] not in screened7
]

assert len(remaining) == 119

hits = defaultdict(list)

for acc in ACCS:
    path = PAV / "asm20_119" / f"Co4_119_vs_{acc}.paf"

    with path.open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue

            qlen = int(p[1])
            qstart = int(p[2])
            qend = int(p[3])
            matches = int(p[9])
            block_length = int(p[10])

            if qlen <= 0 or block_length <= 0:
                continue

            hits[(p[0], acc)].append({
                "target": p[5],
                "coverage": 100 * (qend - qstart) / qlen,
                "identity": 100 * matches / block_length,
                "mapq": int(p[11]),
            })

rows = []
counts = Counter()

for candidate in remaining:
    gene = candidate["gene_id"]
    result = {
        "pgID": candidate["pgID"],
        "gene_id": gene,
        "best_AED": candidate["best_AED"],
    }

    aligned_accs = 0
    strong_accs = 0

    for acc in ACCS:
        current = hits.get((gene, acc), [])

        if current:
            aligned_accs += 1
            best = max(
                current,
                key=lambda h: (
                    h["coverage"] * h["identity"],
                    h["mapq"]
                )
            )
            cov = best["coverage"]
            ident = best["identity"]
            target = best["target"]
        else:
            cov = 0.0
            ident = 0.0
            target = ""

        strong = cov >= 90 and ident >= 90
        strong_accs += int(strong)

        result[f"{acc}_hit_count"] = len(current)
        result[f"{acc}_target"] = target
        result[f"{acc}_coverage_pct"] = round(cov, 2)
        result[f"{acc}_identity_pct"] = round(ident, 2)
        result[f"{acc}_strong_match"] = int(strong)

    if strong_accs:
        status = "strong_match_detected"
    elif aligned_accs:
        status = "partial_or_divergent_match"
    else:
        status = "no_asm20_alignment"

    result["aligned_reference_count"] = aligned_accs
    result["strong_reference_count"] = strong_accs
    result["status"] = status

    rows.append(result)
    counts[status] += 1

output = PAV / "Co4_119_asm20_summary.tsv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

unmapped = [
    r for r in rows
    if r["status"] == "no_asm20_alignment"
]

unmapped_path = PAV / "Co4_119_still_unmapped.tsv"

with unmapped_path.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(unmapped)

report = {
    "initial_candidates": 119,
    "classification": dict(counts),
    "still_unmapped_pangene_groups": len({
        r["pgID"] for r in unmapped
    }),
    "notes": [
        "PAF identity uses matching bases / alignment block length.",
        "Coverage refers to the best single query alignment.",
        "No asm20 alignment is not evidence of genomic absence.",
        "Short genes may be missed by minimap2 presets.",
        "Further sensitive BLASTN and synteny checks are needed."
    ]
}

(PAV / "Co4_119_asm20_summary.json").write_text(
    json.dumps(report, indent=2)
)

print(json.dumps(report, indent=2))
print("\nSaved:", output)
print("Saved:", unmapped_path)
