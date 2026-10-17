#!/bin/bash
#SBATCH --job-name=TEsorter_intact_LTR
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$REPO_DIR/config/local.env}"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Configuration not found: $CONFIG_FILE" >&2
    exit 1
fi
source "$CONFIG_FILE"
: "${COURSEDIR:?COURSEDIR must be configured}"
: "${WORKDIR:?WORKDIR must be configured}"

# Use the intact-element FASTA
# Edit these two input paths to match your assembly's filename.
RAW="$WORKDIR/output/EDTA"
# Expand the wildcard outside quotes;
shopt -s nullglob
FASTA=("$RAW"/*.mod.EDTA.raw/*.mod.LTR.intact.raw.fa)
GFF=("$RAW"/*.mod.EDTA.raw/*.mod.LTR.intact.raw.gff3)

OUTDIR="$WORKDIR/output/LTR_full_length"
unset SINGULARITY_BIND APPTAINER_BIND SINGULARITY_TMPDIR
mkdir -p "$OUTDIR/tmp"
# Work with copies so ongoing EDTA files are not changed.
cp "${FASTA[0]}" "$OUTDIR/intact_LTR.fa"
cp "${GFF[0]}" "$OUTDIR/intact_LTR.gff3"
cd "$OUTDIR"
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" \
    --bind "$OUTDIR/tmp":/tmp --pwd "$OUTDIR" \
    "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter intact_LTR.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK"
# After this job succeeds, run the R script in the output directory:
# cd "$OUTDIR"
# apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$OUTDIR" \
#   "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
#   Rscript --vanilla "$COURSEDIR/scripts/03-full_length_LTRs_identity_copy.R"

# You can also download the data and script and run it locally in RStudio. 
