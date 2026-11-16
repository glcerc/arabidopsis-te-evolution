#!/bin/bash
#SBATCH --job-name=co4_pav_mm2
#SBATCH --array=0-3
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data/users/gercan/Gene_annotation/logs/co4_pav_mm2_%A_%a.out
#SBATCH --error=/data/users/gercan/Gene_annotation/logs/co4_pav_mm2_%A_%a.err

set -euo pipefail

module load minimap2/2.20-GCCcore-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

ROOT=/data/users/gercan/Gene_annotation
REFDIR=/data/courses/genome_annotation_resources/2026/data/Lian_et_al/genome/selected
ASSEMBLY=/data/users/gercan/assembly_annotation_course/assemblies/flye/assembly.fasta
OUT=$ROOT/results/comparative_genomics/GENESPACE/run_5_with_TAIR10/process_pangenome/Co4_specific_functional_analysis/PAV_validation

ACCESSIONS=(Col-0 Altai-5 Are-6 Etna-2)
ACC=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

mkdir -p "$OUT/minimap2_results"

TMPBASE=${SLURM_TMPDIR:-$OUT}
TMP=$(mktemp -d "$TMPBASE/co4mm2_${ACC}_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

bedtools getfasta \
  -fi "$ASSEMBLY" \
  -bed "$OUT/Co4_585_candidate_genes.bed" \
  -nameOnly \
  -fo "$TMP/candidates.fasta"

N=$(grep -c '^>' "$TMP/candidates.fasta")
if [ "$N" -ne 585 ]; then
    echo "ERROR: Expected 585 candidates, found $N"
    exit 1
fi

minimap2 \
  -x asm5 \
  -c \
  -t "$SLURM_CPUS_PER_TASK" \
  "$REFDIR/${ACC}.fasta" \
  "$TMP/candidates.fasta" \
  > "$OUT/minimap2_results/Co4_vs_${ACC}.paf"

echo "Finished: $ACC"
wc -l "$OUT/minimap2_results/Co4_vs_${ACC}.paf"
