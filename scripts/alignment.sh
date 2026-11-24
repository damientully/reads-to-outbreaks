#!/bin/bash
#SBATCH --job-name=vg_alignment
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%j_alignment.out
#SBATCH --error=logs/%j_alignment.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
# -------------------------------------------------------------------

# If the shared folder has no data yet, fall back to a local data/ copy.
[ -d "$RESOURCES/data" ] || RESOURCES="$WORKDIR"

# ---------------------------- ENVIRONMENT --------------------------
module load mafft seqkit
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/alignment logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
GEN=$RESOURCES/data/genomes
OUT=results/alignment
T=${SLURM_CPUS_PER_TASK:-4}

# Relabel: ivar writes 'Consensus_A01_threshold_0.75_quality_20', which is
# unusable on a tree figure. Use your own consensus where you have it, and the
# shared copies for the rest.
: > "$OUT/all.fa"
for S in $(seq -w 1 12); do
  L="A${S}"
  SRC=""
  [ -s "results/consensus/${L}.fa" ] && SRC="results/consensus/${L}.fa"
  [ -z "$SRC" ] && [ -s "$GEN/${L}.fa" ] && SRC="$GEN/${L}.fa"
  [ -z "$SRC" ] && { echo "  no consensus for $L - skipping"; continue; }
  echo ">${L}" >> "$OUT/all.fa"
  grep -v '^>' "$SRC" >> "$OUT/all.fa"
done

# The reference goes in as an outgroup. Without it the tree has no root, and
# 'basal' means nothing.
echo ">REF_Wuhan_Hu_1" >> "$OUT/all.fa"
grep -v '^>' "$REF" >> "$OUT/all.fa"

grep -c '^>' "$OUT/all.fa" | xargs echo "sequences to align:"
mafft --thread "$T" --auto "$OUT/all.fa" > "$OUT/all.aln"
seqkit stats "$OUT/all.aln"
