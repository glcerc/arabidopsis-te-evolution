#!/bin/bash
#SBATCH --job-name=co4_pav7
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data/users/gercan/Gene_annotation/logs/co4_pav7_%A_%a.out
#SBATCH --error=/data/users/gercan/Gene_annotation/logs/co4_pav7_%A_%a.err

set -euo pipefail

module load BLAST+/2.15.0-gompi-2021a
module load BEDTools/2.30.0-GCC-10.3.0

ROOT=/data/users/gercan/Gene_annotation
PAV="$ROOT/results/comparative_genomics/GENESPACE/run_5_with_TAIR10/process_pangenome/Co4_specific_functional_analysis/PAV_validation"
REFDIR=/data/courses/genome_annotation_resources/2026/data/Lian_et_al/genome/selected
ASSEMBLY=/data/users/gercan/assembly_annotation_course/assemblies/flye/assembly.fasta

ACCESSIONS=(Col-0 Altai-5 Are-6 Etna-2)
ACC=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

mkdir -p "$PAV/validation_7"
TMP=$(mktemp -d "${SLURM_TMPDIR:-$PAV}/pav7_${ACC}_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

python3 - "$PAV" "$TMP/seven.bed" <<'PY'
import csv
import sys

from pathlib import Path

pav = Path(sys.argv[1])
ids = set()

with (pav / "Co4_higher_priority_PAV_candidates.tsv").open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        ids.add(row["gene_id"])

with (pav / "Co4_585_candidate_genes.bed").open() as src, \
     open(sys.argv[2], "w") as dst:
    for line in src:
        if line.rstrip("\n").split("\t")[3] in ids:
            dst.write(line)
PY

bedtools getfasta \
  -fi "$ASSEMBLY" \
  -bed "$TMP/seven.bed" \
  -nameOnly \
  -fo "$TMP/seven.fasta"

N=$(grep -c '^>' "$TMP/seven.fasta")
echo "Candidate sequences: $N"
test "$N" -eq 7

makeblastdb \
  -in "$REFDIR/${ACC}.fasta" \
  -dbtype nucl \
  -out "$TMP/refdb" >/dev/null

blastn \
  -task blastn \
  -query "$TMP/seven.fasta" \
  -db "$TMP/refdb" \
  -dust yes \
  -evalue 1e-5 \
  -max_target_seqs 10 \
  -max_hsps 10 \
  -num_threads "$SLURM_CPUS_PER_TASK" \
  -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore" \
  -out "$PAV/validation_7/${ACC}_blastn.tsv"

echo "Completed $ACC"
wc -l "$PAV/validation_7/${ACC}_blastn.tsv"
