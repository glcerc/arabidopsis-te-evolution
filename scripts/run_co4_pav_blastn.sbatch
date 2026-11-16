#!/bin/bash
#SBATCH --job-name=co4_pav
#SBATCH --array=0-3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/data/users/gercan/Gene_annotation/logs/co4_pav_%A_%a.out
#SBATCH --error=/data/users/gercan/Gene_annotation/logs/co4_pav_%A_%a.err

set -euo pipefail

module load BLAST+/2.15.0-gompi-2021a
module load BEDTools/2.30.0-GCC-10.3.0

ROOT=/data/users/gercan/Gene_annotation
ASSEMBLY=/data/users/gercan/assembly_annotation_course/assemblies/flye/assembly.fasta
REFDIR=/data/courses/genome_annotation_resources/2026/data/Lian_et_al/genome/selected
OUT=$ROOT/results/comparative_genomics/GENESPACE/run_5_with_TAIR10/process_pangenome/Co4_specific_functional_analysis/PAV_validation

ACCESSIONS=(Col-0 Altai-5 Are-6 Etna-2)
ACC=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

mkdir -p "$OUT/blast_results"

TMPBASE=${SLURM_TMPDIR:-$OUT}
TMP=$(mktemp -d "$TMPBASE/co4pav_${ACC}_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

echo "Accession: $ACC"
echo "Job ID: $SLURM_JOB_ID"
echo "Temporary directory: $TMP"

cp "$ASSEMBLY" "$TMP/co4_assembly.fasta"
cp "$REFDIR/${ACC}.fasta" "$TMP/reference.fasta"

bedtools getfasta \
  -fi "$TMP/co4_assembly.fasta" \
  -bed "$OUT/Co4_585_candidate_genes.bed" \
  -nameOnly \
  -fo "$TMP/candidates.fasta"

N=$(grep -c '^>' "$TMP/candidates.fasta")
echo "Extracted candidate sequences: $N"

if [ "$N" -ne 585 ]; then
    echo "ERROR: Expected 585 sequences."
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    cp "$TMP/candidates.fasta" "$OUT/Co4_585_gene_regions.fasta"
fi

makeblastdb \
  -in "$TMP/reference.fasta" \
  -dbtype nucl \
  -out "$TMP/reference_db"

blastn \
  -query "$TMP/candidates.fasta" \
  -db "$TMP/reference_db" \
  -task blastn \
  -word_size 11 \
  -evalue 1e-10 \
  -max_target_seqs 20 \
  -max_hsps 20 \
  -num_threads "$SLURM_CPUS_PER_TASK" \
  -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore" \
  -out "$OUT/blast_results/Co4_vs_${ACC}.tsv"

echo "Finished: $ACC"
wc -l "$OUT/blast_results/Co4_vs_${ACC}.tsv"
