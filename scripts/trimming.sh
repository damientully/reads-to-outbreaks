#!/bin/bash
#SBATCH --job-name=vg_trimming
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_trimming.out
#SBATCH --error=logs/%j_trimming.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
SAMPLE=A01
# -------------------------------------------------------------------

# ---------------------------- ENVIRONMENT --------------------------
module load fastp fastqc multiqc
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/trimming logs

IN=results/host_removal
OUT=results/trimming
T=${SLURM_CPUS_PER_TASK:-4}

fastp \
  -i "$IN/${SAMPLE}_R1.fastq.gz" -I "$IN/${SAMPLE}_R2.fastq.gz" \
  -o "$OUT/${SAMPLE}_R1.fastq.gz" -O "$OUT/${SAMPLE}_R2.fastq.gz" \
  --json "$OUT/${SAMPLE}_fastp.json" --html "$OUT/${SAMPLE}_fastp.html" \
  --thread "$T" --cut_tail --cut_tail_mean_quality 20 --length_required 50

fastqc -t "$T" --outdir "$OUT" "$OUT/${SAMPLE}_R1.fastq.gz" "$OUT/${SAMPLE}_R2.fastq.gz"
echo "done: $OUT"
