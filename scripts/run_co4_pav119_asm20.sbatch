#!/bin/bash
#SBATCH --job-name=co4_pav119
#SBATCH --partition=pibu_el8
#SBATCH --array=0-3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/gercan/Gene_annotation/logs/co4_pav119_%A_%a.out
#SBATCH --error=/data/users/gercan/Gene_annotation/logs/co4_pav119_%A_%a.err

set -euo pipefail

module load minimap2/2.20-GCCcore-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

ROOT=/data/users/gercan/Gene_annotation
ASSEMBLY=/data/users/gercan/assembly_annotation_course/assemblies/flye/assembly.fasta
REFDIR=/data/courses/genome_annotation_resources/2026/data/Lian_et_al/genome/selected

PAV="$ROOT/results/comparative_genomics/GENESPACE/run_5_with_TAIR10/process_pangenome/Co4_specific_functional_analysis/PAV_validation"

ACCESSIONS=(Col-0 Altai-5 Are-6 Etna-2)
ACC=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

RESULTS="$PAV/asm20_119"
mkdir -p "$RESULTS"

TMP=$(mktemp -d "${SLURM_TMPDIR:-$RESULTS}/co4_asm20_${ACC}_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

bedtools getfasta \
  -fi "$ASSEMBLY" \
  -bed "$PAV/Co4_119_remaining_candidates.bed" \
  -nameOnly \
  -fo "$TMP/candidates.fasta"

N=$(grep -c '^>' "$TMP/candidates.fasta")

echo "Accession: $ACC"
echo "Candidate sequences: $N"

test "$N" -eq 119

minimap2 \
  -x asm20 \
  -c \
  -t "$SLURM_CPUS_PER_TASK" \
  "$REFDIR/${ACC}.fasta" \
  "$TMP/candidates.fasta" \
  > "$RESULTS/Co4_119_vs_${ACC}.paf"

echo "Completed: $ACC"
wc -l "$RESULTS/Co4_119_vs_${ACC}.paf"
