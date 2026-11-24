#!/bin/bash
#SBATCH --job-name=vg_consensus
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_consensus.out
#SBATCH --error=logs/%j_consensus.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
SAMPLE=A01
MIN_DEPTH=20
MIN_FREQ=0.75
# -------------------------------------------------------------------

# If the shared folder has no data yet, fall back to a local data/ copy.
[ -d "$RESOURCES/data" ] || RESOURCES="$WORKDIR"

# ---------------------------- ENVIRONMENT --------------------------
module load samtools ivar seqkit
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/consensus results/variants logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
BAM=results/mapping/${SAMPLE}.sorted.bam

# consensus
samtools mpileup -aa -A -d 0 -Q 0 --reference "$REF" "$BAM" \
  | ivar consensus -p results/consensus/${SAMPLE} \
      -t "$MIN_FREQ" -m "$MIN_DEPTH" -n N

# variants, including minor variants down to 3%
samtools mpileup -aa -A -d 0 -B -Q 0 --reference "$REF" "$BAM" \
  | ivar variants -p results/variants/${SAMPLE} -r "$REF" -t 0.03 -m "$MIN_DEPTH"

# how much of the genome did we actually observe?
seqkit fx2tab -nlg results/consensus/${SAMPLE}.fa
NS=$(grep -v '^>' results/consensus/${SAMPLE}.fa | tr -cd 'Nn' | wc -c)
LEN=$(grep -v '^>' results/consensus/${SAMPLE}.fa | tr -d '\n' | wc -c)
awk -v n="$NS" -v l="$LEN" 'BEGIN {printf "Ns=%d / %d (%.2f%%)\n", n, l, 100*n/l}' \
  | tee results/consensus/${SAMPLE}.ncount.txt
echo "MIN_DEPTH=$MIN_DEPTH drives that N count. Raise it and the genome gets smaller but safer."
