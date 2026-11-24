#!/bin/bash
#SBATCH --job-name=vg_qc
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_qc.out
#SBATCH --error=logs/%j_qc.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
SAMPLE=A01                      # A01 .. A12
# -------------------------------------------------------------------

# If the shared folder has no data yet, fall back to a local data/ copy.
[ -d "$RESOURCES/data" ] || RESOURCES="$WORKDIR"

# ---------------------------- ENVIRONMENT --------------------------
module load fastqc multiqc
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/qc logs

R1=$RESOURCES/data/fastq/${SAMPLE}_R1.fastq.gz
R2=$RESOURCES/data/fastq/${SAMPLE}_R2.fastq.gz

# NOTE: fastqc 0.12 dropped '-o'. Use --outdir.
fastqc -t "${SLURM_CPUS_PER_TASK:-4}" --outdir results/qc "$R1" "$R2"
multiqc -o results/qc -n ${SAMPLE}_multiqc results/qc

echo "done: results/qc"
