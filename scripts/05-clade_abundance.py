
#!/usr/bin/env python3
"""Estimate genomic annotation abundance per Copia/Gypsy clade."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def load_classifications(paths):
    classifications = {}

    for path in paths:
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                te_id = row["#TE"].split("#")[0]
                superfamily = row["Superfamily"].strip()
                clade = row["Clade"].strip() or "Unclassified"

                if te_id in classifications:
                    raise ValueError(f"Duplicate TE ID: {te_id}")

                classifications[te_id] = (superfamily, clade)

    return classifications


def parse_attributes(text):
    attributes = {}

    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value

    return attributes


def analyze_gff(gff_path, classifications):
    counts = Counter()
    bases = Counter()
    representatives = defaultdict(set)
    unmatched = Counter()
    mismatched = Counter()

    with open(gff_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            attributes = parse_attributes(fields[8])
            genome_class = attributes.get("classification", "")

            if genome_class not in ("LTR/Copia", "LTR/Gypsy"):
                continue

            te_id = attributes.get("Name", "").split("#")[0]

            if te_id not in classifications:
                unmatched[genome_class] += 1
                continue

            superfamily, clade = classifications[te_id]

            if genome_class != f"LTR/{superfamily}":
                mismatched[genome_class] += 1
                continue

            start = int(fields[3])
            end = int(fields[4])
            length = end - start + 1

            key = (superfamily, clade)
            counts[key] += 1
            bases[key] += length
            representatives[key].add(te_id)

    return counts, bases, representatives, unmatched, mismatched


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--copia", required=True, type=Path)
    parser.add_argument("--gypsy", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)

    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    classifications = load_classifications(
        [args.copia, args.gypsy]
    )

    counts, bases, representatives, unmatched, mismatched = analyze_gff(
        args.gff, classifications
    )

    output_file = args.output / "clade_abundance.tsv"

    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "Superfamily",
            "Clade",
            "Library_representatives_with_annotations",
            "Genomic_annotation_count",
            "Summed_annotation_bp",
        ])

        for key in sorted(
            counts,
            key=lambda k: (k[0], -counts[k], k[1])
        ):
            writer.writerow([
                key[0],
                key[1],
                len(representatives[key]),
                counts[key],
                bases[key],
            ])

    print(f"Output: {output_file}")
    print(f"Library representatives classified: {len(classifications)}")

    for family in ("Copia", "Gypsy"):
        total = sum(
            count for (sf, _), count in counts.items()
            if sf == family
        )
        print(f"{family} matched genomic annotations: {total}")
        print(
            f"{family} unmatched genomic annotations: "
            f"{unmatched.get('LTR/' + family, 0)}"
        )
        print(
            f"{family} classification mismatches: "
            f"{mismatched.get('LTR/' + family, 0)}"
        )

    print("\nTop clades:")
    for (family, clade), count in counts.most_common(10):
        print(f"{family:8} {clade:15} {count:6}")


if __name__ == "__main__":
    main()

