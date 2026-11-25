#!/bin/bash
# Logs are written beside wherever you submit from, as
#   vg_<step>_<arrayjobid>_<task>.out/.err
# The task number IS the sample number, so _7 is A07.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_trimming
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-12
# 2 CPUs x 12 tasks = 24, which fits under the 32-CPU per-user QOS limit on this
# cluster, so all twelve run in one wave. At 4 CPUs only 8 run and the rest sit
# in PD with QOSMaxCpuPerUserLimit. The tools are I/O bound, not CPU bound.
#SBATCH --time=01:00:00
#SBATCH --output=vg_trimming_%A_%a.out
#SBATCH --error=vg_trimming_%A_%a.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
# This script runs ALL TWELVE samples by default, as a job array (see
# --array=1-12 above). SAMPLE below is a fallback, used only if you submit
# without an array, e.g.  sbatch --array=1 scripts/qc.sh
SAMPLE=A01
# -------------------------------------------------------------------

# ---------------------------- ENVIRONMENT --------------------------
module load fastp fastqc multiqc
# -------------------------------------------------------------------

set -euo pipefail

# An array task processes one sample; SAMPLE above is used only if you submit
# this without --array.
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  SAMPLE=$(printf "A%02d" "$SLURM_ARRAY_TASK_ID")
fi
echo "### $SAMPLE"
cd "$WORKDIR"
mkdir -p results/trimming

IN=results/host_removal
OUT=results/trimming
T=${SLURM_CPUS_PER_TASK:-4}

# ReadItAndKeep writes <prefix>.reads_1.fastq.gz; older versions of this
# pipeline used <prefix>_R1.fastq.gz. Accept either, and stop if neither is there
# rather than handing an empty file to the next step.
IN1=$IN/${SAMPLE}.reads_1.fastq.gz
IN2=$IN/${SAMPLE}.reads_2.fastq.gz
if [ ! -s "$IN1" ]; then
  IN1=$IN/${SAMPLE}_R1.fastq.gz
  IN2=$IN/${SAMPLE}_R2.fastq.gz
fi
if [ ! -s "$IN1" ] || [ ! -s "$IN2" ]; then
  echo "ERROR: no host-depleted reads for $SAMPLE in $IN/" >&2
  echo "       expected ${SAMPLE}.reads_1.fastq.gz - has host_removal.sh finished?" >&2
  ls -la "$IN" >&2 || true
  exit 1
fi

fastp \
  -i "$IN1" -I "$IN2" \
  -o "$OUT/${SAMPLE}_R1.fastq.gz" -O "$OUT/${SAMPLE}_R2.fastq.gz" \
  --json "$OUT/${SAMPLE}_fastp.json" --html "$OUT/${SAMPLE}_fastp.html" \
  --thread "$T" --cut_tail --cut_tail_mean_quality 20 --length_required 50

fastqc -t "$T" --outdir "$OUT" "$OUT/${SAMPLE}_R1.fastq.gz" "$OUT/${SAMPLE}_R2.fastq.gz"
echo "done: $OUT"
