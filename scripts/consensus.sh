#!/bin/bash
# Logs are written beside wherever you submit from, as vg_<step>_<jobid>.out/.err.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_consensus
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-12
#SBATCH --time=01:00:00
#SBATCH --output=vg_consensus_%j.out
#SBATCH --error=vg_consensus_%j.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
# This script runs ALL TWELVE samples by default, as a job array (see
# --array=1-12 above). SAMPLE below is a fallback, used only if you submit
# without an array, e.g.  sbatch --array=1 scripts/qc.sh
SAMPLE=A01
MIN_DEPTH=20
MIN_FREQ=0.75
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
module load samtools ivar seqkit
# -------------------------------------------------------------------

set -euo pipefail

# An array task processes one sample; SAMPLE above is used only if you submit
# this without --array.
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  SAMPLE=$(printf "A%02d" "$SLURM_ARRAY_TASK_ID")
fi
echo "### $SAMPLE"
cd "$WORKDIR"
mkdir -p results/consensus results/variants logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
BAM=results/mapping/${SAMPLE}.sorted.bam

# An empty BAM is 312 bytes of header and produces a genome of Ns. Catch it.
if [ ! -s "$BAM" ] || [ "$(samtools view -c "$BAM")" -eq 0 ]; then
  echo "ERROR: $BAM has no alignments - has mapping.sh finished successfully?" >&2
  exit 1
fi

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
