#!/bin/bash
# Logs are written beside wherever you submit from, as vg_<step>_<jobid>.out/.err.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_alignment
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=vg_alignment_%j.out
#SBATCH --error=vg_alignment_%j.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
# -------------------------------------------------------------------

# Resolve where the reads actually are: the shared bundle, or a local copy.
# Test for a real file, not just a directory - an empty or half-staged shared
# folder must not look populated.
if [ ! -s "$RESOURCES/data/fastq/A01_R1.fastq.gz" ]; then
  if [ -s "$WORKDIR/data/fastq/A01_R1.fastq.gz" ]; then
    RESOURCES="$WORKDIR"
  else
    cat >&2 <<ERR

ERROR: cannot find the sequencing data.

  Looked for:
    $RESOURCES/data/fastq/A01_R1.fastq.gz
    $WORKDIR/data/fastq/A01_R1.fastq.gz

  If you are a participant: tell the facilitator the bundle is not staged.

  If you are the facilitator: stage it with build/stage_bundle.sh, or for a
  local test put the reads in \$WORKDIR/data/fastq/ (note the fastq/
  subdirectory - flat files in data/ will not be found).
ERR
    exit 1
  fi
fi

# ---------------------------- ENVIRONMENT --------------------------
module load mafft seqkit
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/alignment

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
