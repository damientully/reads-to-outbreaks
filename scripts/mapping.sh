#!/bin/bash
# Logs are written beside wherever you submit from, as
#   vg_<step>_<arrayjobid>_<task>.out/.err
# The task number IS the sample number, so _7 is A07.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_mapping
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-12
# 2 CPUs x 12 tasks = 24, which fits under the 32-CPU per-user QOS limit on this
# cluster, so all twelve run in one wave. At 4 CPUs only 8 run and the rest sit
# in PD with QOSMaxCpuPerUserLimit. The tools are I/O bound, not CPU bound.
#SBATCH --time=01:00:00
#SBATCH --output=vg_mapping_%A_%a.out
#SBATCH --error=vg_mapping_%A_%a.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
# This script runs ALL TWELVE samples by default, as a job array (see
# --array=1-12 above). SAMPLE below is a fallback, used only if you submit
# without an array, e.g.  sbatch --array=1 scripts/qc.sh
SAMPLE=A01
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
module load minimap2 samtools
# -------------------------------------------------------------------

set -euo pipefail

# An array task processes one sample; SAMPLE above is used only if you submit
# this without --array.
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  SAMPLE=$(printf "A%02d" "$SLURM_ARRAY_TASK_ID")
fi
echo "### $SAMPLE"
cd "$WORKDIR"
mkdir -p results/mapping logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
IN=results/trimming
OUT=results/mapping
T=${SLURM_CPUS_PER_TASK:-4}

for f in "$IN/${SAMPLE}_R1.fastq.gz" "$IN/${SAMPLE}_R2.fastq.gz"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty $f - has trimming.sh finished?" >&2; exit 1; }
done

minimap2 -t "$T" -ax sr "$REF" "$IN/${SAMPLE}_R1.fastq.gz" "$IN/${SAMPLE}_R2.fastq.gz" \
  | samtools sort -@ "$T" -o "$OUT/${SAMPLE}.sorted.bam"
samtools index "$OUT/${SAMPLE}.sorted.bam"
samtools flagstat "$OUT/${SAMPLE}.sorted.bam" > "$OUT/${SAMPLE}.flagstat"
samtools depth -a -d 0 "$OUT/${SAMPLE}.sorted.bam" > "$OUT/${SAMPLE}.depth.tsv"

awk '{s+=$3; if($3<10) n++} END {printf "mean_depth=%.0f  positions_below_10x=%d\n", s/NR, n+0}' \
  "$OUT/${SAMPLE}.depth.tsv" | tee "$OUT/${SAMPLE}.depth_summary.txt"
