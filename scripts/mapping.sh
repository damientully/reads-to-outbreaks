#!/bin/bash
#SBATCH --job-name=vg_mapping
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_mapping.out
#SBATCH --error=logs/%j_mapping.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
SAMPLE=A01
# -------------------------------------------------------------------

# If the shared folder has no data yet, fall back to a local data/ copy.
[ -d "$RESOURCES/data" ] || RESOURCES="$WORKDIR"

# ---------------------------- ENVIRONMENT --------------------------
module load minimap2 samtools
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/mapping logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
IN=results/trimming
OUT=results/mapping
T=${SLURM_CPUS_PER_TASK:-4}

minimap2 -t "$T" -ax sr "$REF" "$IN/${SAMPLE}_R1.fastq.gz" "$IN/${SAMPLE}_R2.fastq.gz" \
  | samtools sort -@ "$T" -o "$OUT/${SAMPLE}.sorted.bam"
samtools index "$OUT/${SAMPLE}.sorted.bam"
samtools flagstat "$OUT/${SAMPLE}.sorted.bam" > "$OUT/${SAMPLE}.flagstat"
samtools depth -a -d 0 "$OUT/${SAMPLE}.sorted.bam" > "$OUT/${SAMPLE}.depth.tsv"

awk '{s+=$3; if($3<10) n++} END {printf "mean_depth=%.0f  positions_below_10x=%d\n", s/NR, n+0}' \
  "$OUT/${SAMPLE}.depth.tsv" | tee "$OUT/${SAMPLE}.depth_summary.txt"
