#!/bin/bash
#SBATCH --job-name=EDTA_annotation
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --time=6-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Stop on failed commands, unset variables, or failures within pipelines.
set -euo pipefail

# Copy this script to your project and submit from that directory:
# mkdir -p logs
# sbatch 01-run_EDTA.sh
# Shared course software and references; results stay in your own project.
# Load user-specific configuration
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$REPO_DIR/config/local.env}"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Configuration not found: $CONFIG_FILE" >&2
    exit 1
fi
source "$CONFIG_FILE"

: "${COURSEDIR:?COURSEDIR must be configured}"
: "${WORKDIR:?WORKDIR must be configured}"
: "${GENOME:?GENOME must be configured}"

GENOME="$(readlink -f "$GENOME")"
# Arabidopsis coding sequences help remove gene-derived sequences from the TE library.
CDS="$COURSEDIR/data/TAIR10_cds_20110103_representative_gene_model_updated"
EDTA_SIF="$COURSEDIR/containers/EDTA_v2.3.0_TIRLearner4.sif"
OUTDIR="$WORKDIR/output/EDTA"

# Prevent inherited container settings from changing the course environment.
unset SINGULARITY_BIND APPTAINER_BIND SINGULARITY_TMPDIR
# Keep temporary files on disk, outside the container's small temporary filesystem.
mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

# Bind exposes host directories inside the container; --pwd sets its working directory.
# Run the shared EDTA checkout, which supports the container's TIR-Learner 4.
# --step all builds the TE library; --anno 1 also annotates the genome.
# --sensitive 1 includes RepeatModeler families; --wholeelement 1 retains complete LTR elements
apptainer exec -C \
    --bind "$WORKDIR" --bind "$COURSEDIR" \
    --bind "$(dirname "$GENOME")" --bind "$(dirname "$CDS")" \
    --bind "$OUTDIR/tmp":/tmp --pwd "$OUTDIR" "$EDTA_SIF" \
    perl "$COURSEDIR/tools/EDTA_master/EDTA.pl" \
    --genome "$GENOME" \
    --species others \
    --step all \
    --sensitive 1 \
    --cds "$CDS" \
    --anno 1 \
    --wholeelement 1 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --overwrite 0
