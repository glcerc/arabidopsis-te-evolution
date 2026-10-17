#!/bin/bash
#SBATCH --job-name=TEsorter_library
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
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


OUTDIR="$WORKDIR/output/TEsorter_library"
LIBRARY="$WORKDIR/output/EDTA/assembly.fasta.mod.EDTA.TElib.fa"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Extract Copia and Gypsy library sequences
apptainer exec -C --bind "$WORKDIR" \
  "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  seqkit grep -r -p "Copia" "$LIBRARY" > Copia_sequences.fa

apptainer exec -C --bind "$WORKDIR" \
  "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  seqkit grep -r -p "Gypsy" "$LIBRARY" > Gypsy_sequences.fa

# Verify that both inputs contain sequences
test -s Copia_sequences.fa
test -s Gypsy_sequences.fa

echo "Copia representative sequences:"
grep -c '^>' Copia_sequences.fa

echo "Gypsy representative sequences:"
grep -c '^>' Gypsy_sequences.fa

# Run library-level classification
apptainer exec -C --bind "$WORKDIR" --pwd "$OUTDIR" \
  "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
  TEsorter Copia_sequences.fa -db rexdb-plant -p 4

apptainer exec -C --bind "$WORKDIR" --pwd "$OUTDIR" \
  "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
  TEsorter Gypsy_sequences.fa -db rexdb-plant -p 4

echo "Library-level TEsorter analysis completed."
