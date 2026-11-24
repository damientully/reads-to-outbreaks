#!/bin/bash
#SBATCH --job-name=vg_host_removal
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_host_removal.out
#SBATCH --error=logs/%j_host_removal.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
SAMPLE=A01
# -------------------------------------------------------------------

# If the shared folder has no data yet, fall back to a local data/ copy.
[ -d "$RESOURCES/data" ] || RESOURCES="$WORKDIR"

# ---------------------------- ENVIRONMENT --------------------------
module load read_keep minimap2 samtools
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/host_removal logs

REF=$RESOURCES/data/reference/MN908947.3.fasta
R1=$RESOURCES/data/fastq/${SAMPLE}_R1.fastq.gz
R2=$RESOURCES/data/fastq/${SAMPLE}_R2.fastq.gz
OUT=results/host_removal
T=${SLURM_CPUS_PER_TASK:-4}

# These are clinical nasopharyngeal samples sequenced without targeted
# enrichment, so most of what came off the machine is patient, not virus.
# ReadItAndKeep keeps only reads matching the SARS-CoV-2 reference.

# The bioconda binary is camelCase; check what the module actually provides.
RIAK=$(command -v readItAndKeep || command -v read_keep || command -v readitandkeep || true)
[ -n "$RIAK" ] || { echo "ERROR: ReadItAndKeep not on PATH after 'module load read_keep'"; exit 1; }

"$RIAK" --tech illumina --ref_fasta "$REF" \
  --reads1 "$R1" --reads2 "$R2" \
  --outprefix "$OUT/${SAMPLE}" 2>&1 | tee "$OUT/${SAMPLE}.riak.log"

# ReadItAndKeep writes <prefix>.reads_1.fastq.gz / .reads_2.fastq.gz.
# Check the names on your version before relying on them downstream.
ls -lh "$OUT/${SAMPLE}".reads_*.fastq.gz

IN_N=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$OUT/${SAMPLE}.riak.log")
KEPT=$(awk -F'\t' '/^Kept reads 1/{print $2}'       "$OUT/${SAMPLE}.riak.log")
awk -v i="$IN_N" -v k="$KEPT" -v s="$SAMPLE" \
  'BEGIN {printf "%s  input_pairs=%s  kept_pairs=%s  retained=%.1f%%\n", s, i, k, 100*k/i}' \
  | tee "$OUT/${SAMPLE}.summary.txt"

cat <<'TXT'

The retained fraction is the number that matters. It varies several-fold between
samples from this same study, because it tracks how much virus was in the swab.
Compare your sample with a neighbour's before moving on.

Two reasons this step exists. Host sequence from a clinical sample is
identifiable patient genetic information, so removing it before data are shared
is a governance requirement, not housekeeping. And everything downstream now
runs on a fraction of the data.

Under the hood this is positive selection on the viral reference. The same thing
can be done with the tools you already know, which is worth seeing once:

  minimap2 -ax sr $REF $R1 $R2 | samtools view -b -F 12 - \
    | samtools sort -n -o paired_mapped.bam
  samtools fastq -1 out_R1.fastq.gz -2 out_R2.fastq.gz -0 /dev/null -s /dev/null \
    paired_mapped.bam

That variant is stricter (it requires BOTH mates to map), so it retains fewer
reads than ReadItAndKeep. Do not mix the two within one dataset.
TXT
