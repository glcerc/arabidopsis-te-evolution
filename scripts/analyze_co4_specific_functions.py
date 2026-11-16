#!/usr/bin/env python3
"""Functional and quality assessment of candidate Co-4-specific genes."""

import csv
import json
import re
from collections import defaultdict, Counter
from pathlib import Path

ROOT = Path("/data/users/gercan/Gene_annotation")
FA = ROOT / "results/functional_annotation"
GS = ROOT / "results/comparative_genomics/GENESPACE/run_5_with_TAIR10"
OUT = GS / "process_pangenome/Co4_specific_functional_analysis"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES = (
    GS / "process_pangenome/tables/Co4_only_from_Co_4_reference.tsv"
)
GFF = FA / "integrated/Co4_final_functional_annotation.gff3"

SOURCES = {
    "UniProt": FA / "BLASTP/blastp_uniprot.besthits.tsv",
    "TAIR10": FA / "BLASTP/blastp_TAIR10.besthits.tsv",
    "Pfam": FA / "InterProScan/outputs/Co4_full_Pfam_AED_0.50.tsv",
    "InterPro": FA / "integrated/Co4_InterPro_annotation_by_transcript.tsv",
    "UniProt_functions": (
        FA / "UniProt_functions/Co4_UniProt_mRNA_annotations.tsv"
    ),
}

def attrs(s):
    return dict(
        part.split("=", 1)
        for part in s.split(";")
        if "=" in part
    )

def normalize(s):
    s = s.strip().split()[0].split("|")[0]
    return re.sub(r"(-protein-\d+|-p\d+)$", "", s)

def write_tsv(path, rows, columns):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns, delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

# Read all candidate records.
with CANDIDATES.open() as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

candidate_genes = {r["gene_id"] for r in candidates}

# Extract gene-to-transcript relationships from final GFF3.
genes = {}
transcripts = defaultdict(list)
transcript_attrs = {}
transcript_to_gene = {}

with GFF.open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, att = fields
        a = attrs(att)
        identifier = a.get("ID", "")

        if feature == "gene":
            genes[identifier] = {
                "chr": seqid,
                "start": int(start),
                "end": int(end),
                "gene_length_bp": int(end) - int(start) + 1,
            }

        elif feature in ("mRNA", "transcript"):
            parents = a.get("Parent", "").split(",")
            transcript_attrs[identifier] = a

            for parent in parents:
                transcripts[parent].append(identifier)
                transcript_to_gene[identifier] = parent

# Add all transcript IDs and protein naming variants as lookup keys.
lookup = defaultdict(set)

for gene in candidate_genes:
    lookup[gene].add(gene)
    for tid in transcripts.get(gene, []):
        for key in {
            tid,
            normalize(tid),
            tid + "-p1",
            tid + "-protein-1",
        }:
            lookup[key].add(gene)

# Assign evidence from functional annotation tables.
evidence = defaultdict(lambda: defaultdict(set))
source_stats = {}

for source, path in SOURCES.items():
    if not path.exists():
        source_stats[source] = {"status": "missing"}
        continue

    matches = 0
    lines = 0
    examples = []
    header_skipped = False

    with path.open(errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue

            query = parts[0].strip()

            # Ignore common TSV column headers.
            if query.lower() in {
                "query", "query_id", "protein_id",
                "transcript_id", "mrna_id", "id",
                "sequence_id", "gene_id", "qseqid",
                "protein", "transcript",
            }:
                header_skipped = True
                continue

            lines += 1
            mapped = set()
            for key in {query, normalize(query)}:
                mapped.update(lookup.get(key, set()))

            if mapped:
                matches += 1
                description = " | ".join(parts[1:])[:1500]

                for gene in mapped:
                    evidence[gene][source].add(description)

                if len(examples) < 2:
                    examples.append(query)

    source_stats[source] = {
        "status": "read",
        "records": lines,
        "matched_records": matches,
        "header_skipped": header_skipped,
        "matched_query_examples": examples,
    }

# Putative categories are only keyword-based screening results.
keywords = {
    "defense_related": (
        r"defen[sc]e|resistance|pathogen|immune|"
        r"nucleotide.binding.*leucine.rich|NLR|TIR.NBS"
    ),
    "stress_related": (
        r"stress|drought|salinity|heat.shock|"
        r"cold.response|oxidative|ABA.response"
    ),
    "transposon_related": (
        r"transpos|retrotranspos|reverse.transcriptase|"
        r"integrase|gag.pol"
    ),
}

rows = []
category_counts = Counter()

for candidate in candidates:
    gene = candidate["gene_id"]
    tids = sorted(set(transcripts.get(gene, [])))
    aeds = []

    for tid in tids:
        value = transcript_attrs.get(tid, {}).get("_AED", "")
        try:
            aeds.append(float(value))
        except ValueError:
            pass

    hit = evidence.get(gene, {})
    searchable = " ".join(
        value
        for source_values in hit.values()
        for value in source_values
    )

    categories = [
        name for name, pattern in keywords.items()
        if re.search(pattern, searchable, re.IGNORECASE)
    ]
    for category in categories:
        category_counts[category] += 1

    meta = genes.get(gene, {})

    row = {
        "pgID": candidate["pgID"],
        "gene_id": gene,
        "GFF_match": int(gene in genes),
        "transcript_count": len(tids),
        "transcript_ids": ";".join(tids),
        "chr": meta.get("chr", ""),
        "start": meta.get("start", ""),
        "end": meta.get("end", ""),
        "gene_length_bp": meta.get("gene_length_bp", ""),
        "best_AED": min(aeds) if aeds else "",
        "UniProt": int(bool(hit.get("UniProt"))),
        "TAIR10": int(bool(hit.get("TAIR10"))),
        "Pfam": int(bool(hit.get("Pfam"))),
        "InterPro": int(bool(hit.get("InterPro"))),
        "UniProt_functions": int(bool(hit.get("UniProt_functions"))),
        "keyword_categories": ";".join(categories),
        "annotation_evidence": searchable[:3000],
    }
    rows.append(row)

columns = [
    "pgID", "gene_id", "GFF_match", "transcript_count",
    "transcript_ids", "chr", "start", "end", "gene_length_bp",
    "best_AED", "UniProt", "TAIR10", "Pfam", "InterPro",
    "UniProt_functions", "keyword_categories",
    "annotation_evidence",
]

write_tsv(
    OUT / "Co4_specific_functional_annotations.tsv",
    rows, columns
)

write_tsv(
    OUT / "Co4_specific_unmatched_genes.tsv",
    [r for r in rows if not r["GFF_match"]],
    columns,
)

write_tsv(
    OUT / "Co4_specific_functional_candidates.tsv",
    [r for r in rows if r["keyword_categories"]],
    columns,
)

summary = {
    "candidate_pangene_groups": len(
        {r["pgID"] for r in candidates}
    ),
    "candidate_gene_records": len(rows),
    "unique_candidate_gene_ids": len(candidate_genes),
    "genes_in_final_GFF": sum(r["GFF_match"] for r in rows),
    "genes_with_transcripts": sum(
        r["transcript_count"] > 0 for r in rows
    ),
    "genes_with_AED": sum(
        r["best_AED"] != "" for r in rows
    ),
    "genes_with_UniProt": sum(r["UniProt"] for r in rows),
    "genes_with_TAIR10": sum(r["TAIR10"] for r in rows),
    "genes_with_Pfam": sum(r["Pfam"] for r in rows),
    "genes_with_InterPro": sum(r["InterPro"] for r in rows),
    "genes_with_UniProt_functions": sum(
        r["UniProt_functions"] for r in rows
    ),
    "genes_with_any_external_evidence": sum(
        any(r[s] for s in SOURCES) for r in rows
    ),
    "keyword_category_counts": dict(category_counts),
    "source_read_statistics": source_stats,
    "limitations": [
        "Only candidate accession-specific groups, not proven unique genes.",
        "Functional keyword categories are screening labels, not validated functions.",
        "Source tables must be checked for query-ID and column compatibility.",
        "Missing annotation may reflect an ID mismatch or absent evidence.",
        "Protein isoforms and pangene groups are not interchangeable units.",
    ],
}

(OUT / "summary.json").write_text(
    json.dumps(summary, indent=2)
)

print("\n=== Co-4-specific functional analysis ===")
for key, value in summary.items():
    if key not in ("source_read_statistics", "limitations"):
        print(f"{key}: {value}")

print("\n=== Source matching diagnostics ===")
for source, info in source_stats.items():
    print(source, json.dumps(info))

print("\nResults saved to:", OUT)
