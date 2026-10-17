#!/bin/bash
#SBATCH --job-name=LTR_identity
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
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


OUTDIR="$WORKDIR/output/LTR_full_length"

cd "$OUTDIR"

apptainer exec -C \
  --bind "$WORKDIR" \
  --bind "$COURSEDIR" \
  --pwd "$OUTDIR" \
  "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  Rscript --vanilla \
  "$COURSEDIR/scripts/03-full_length_LTRs_identity_copy.R"
