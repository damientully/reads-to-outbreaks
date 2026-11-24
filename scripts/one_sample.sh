#!/bin/bash
# Logs are written beside wherever you submit from, as vg_<step>_<jobid>.out/.err.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_one
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=vg_one_sample_%j.out
#SBATCH --error=vg_one_sample_%j.err
#
# All twelve as a job array instead of one sample:
#   sbatch --array=1-12 scripts/one_sample.sh
#
# One complete sample, start to finish: host removal -> QC -> trim -> map ->
# depth -> consensus -> variants, with timings and a comparison against the
# published metrics. Use it as a smoke test after any change to the cluster,
# the modules, or the bundle.

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
RESOURCES=/mnt/roche/resources/viral-genomics
SAMPLE=A07                      # ignored when run as an array
MIN_DEPTH=20
MIN_FREQ=0.75
# -------------------------------------------------------------------

# ---------------------------- ENVIRONMENT --------------------------
module purge
module load read_keep fastqc fastp minimap2 samtools ivar seqkit
# -------------------------------------------------------------------

set -euo pipefail

# an array job overrides SAMPLE with its index
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  SAMPLE=$(printf "A%02d" "$SLURM_ARRAY_TASK_ID")
fi

cd "$WORKDIR"
T=${SLURM_CPUS_PER_TASK:-4}

# Resolve where the reads are: shared bundle, or a local copy. Test for a real
# file, not a directory - a half-staged folder must not look populated.
if [ ! -s "$RESOURCES/data/fastq/${SAMPLE}_R1.fastq.gz" ]; then
  if [ -s "$WORKDIR/data/fastq/${SAMPLE}_R1.fastq.gz" ]; then
    RESOURCES="$WORKDIR"
  else
    echo "ERROR: no reads for $SAMPLE in $RESOURCES/data/fastq/ or $WORKDIR/data/fastq/" >&2
    exit 1
  fi
fi

REF=$RESOURCES/data/reference/MN908947.3.fasta
[ -s "$REF" ] || { echo "ERROR: reference not found at $REF" >&2; exit 1; }

OUT=results/one_sample/$SAMPLE
mkdir -p "$OUT"/{qc,hostdepleted,trimmed,bam,depth,consensus,variants}

echo "=============================================="
echo " sample     : $SAMPLE"
echo " node       : $(hostname)"
echo " threads    : $T"
echo " reads from : $RESOURCES/data/fastq"
echo " output to  : $WORKDIR/$OUT"
echo "=============================================="

step () { echo; echo "--- $* ---"; }

step "1. host removal (read_keep)"
time readItAndKeep --tech illumina --ref_fasta "$REF" \
  --reads1 "$RESOURCES/data/fastq/${SAMPLE}_R1.fastq.gz" \
  --reads2 "$RESOURCES/data/fastq/${SAMPLE}_R2.fastq.gz" \
  --outprefix "$OUT/hostdepleted/$SAMPLE" 2>&1 | tee "$OUT/${SAMPLE}.riak.log"

IN1=$OUT/hostdepleted/${SAMPLE}.reads_1.fastq.gz
IN2=$OUT/hostdepleted/${SAMPLE}.reads_2.fastq.gz
[ -s "$IN1" ] || { echo "ERROR: expected $IN1 - check the output names of this readItAndKeep version" >&2; exit 1; }

step "2. read QC (fastqc) -- note --outdir, not -o"
time fastqc -t "$T" --outdir "$OUT/qc" "$IN1" "$IN2"

step "3. trim (fastp)"
time fastp -i "$IN1" -I "$IN2" \
  -o "$OUT/trimmed/${SAMPLE}_R1.fastq.gz" -O "$OUT/trimmed/${SAMPLE}_R2.fastq.gz" \
  --json "$OUT/qc/${SAMPLE}_fastp.json" --html "$OUT/qc/${SAMPLE}_fastp.html" \
  --thread "$T" --cut_tail --cut_tail_mean_quality 20 --length_required 50

step "4. map (minimap2 | samtools sort)"
time minimap2 -t "$T" -ax sr "$REF" \
    "$OUT/trimmed/${SAMPLE}_R1.fastq.gz" "$OUT/trimmed/${SAMPLE}_R2.fastq.gz" \
  | samtools sort -@ "$T" -o "$OUT/bam/${SAMPLE}.sorted.bam"
samtools index "$OUT/bam/${SAMPLE}.sorted.bam"
samtools flagstat "$OUT/bam/${SAMPLE}.sorted.bam" | tee "$OUT/qc/${SAMPLE}.flagstat" | head -5

step "5. depth"
samtools depth -a -d 0 "$OUT/bam/${SAMPLE}.sorted.bam" > "$OUT/depth/${SAMPLE}.depth.tsv"

step "6. consensus and variants (ivar)"
time samtools mpileup -aa -A -d 0 -Q 0 --reference "$REF" "$OUT/bam/${SAMPLE}.sorted.bam" \
  | ivar consensus -p "$OUT/consensus/$SAMPLE" -t "$MIN_FREQ" -m "$MIN_DEPTH" -n N
samtools mpileup -aa -A -d 0 -B -Q 0 --reference "$REF" "$OUT/bam/${SAMPLE}.sorted.bam" \
  | ivar variants -p "$OUT/variants/$SAMPLE" -r "$REF" -t 0.03 -m "$MIN_DEPTH"

# ------------------------- observed vs published -------------------------
IN_N=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$OUT/${SAMPLE}.riak.log")
KEPT=$(awk -F'\t' '/^Kept reads 1/{print $2}'       "$OUT/${SAMPLE}.riak.log")
OBS_DEPTH=$(awk '{s+=$3} END {printf "%.0f", s/NR}' "$OUT/depth/${SAMPLE}.depth.tsv")
OBS_NS=$(grep -v '^>' "$OUT/consensus/${SAMPLE}.fa" | tr -cd 'Nn' | wc -c | tr -d ' ')
OBS_MAP=$(awk '/ mapped \(/{gsub(/[(%]/,"");print $5; exit}' "$OUT/qc/${SAMPLE}.flagstat")
NVAR=$(( $(wc -l < "$OUT/variants/${SAMPLE}.tsv") - 1 ))

# published figures, looked up by column name so added columns do not break it
TSV=$RESOURCES/data/samples.tsv
read -r PUB_DEPTH PUB_NS <<< "$(awk -F'\t' -v s="$SAMPLE" '
  NR==1 { for (i=1;i<=NF;i++) { if ($i=="mean_depth") d=i; if ($i=="ambig_ns") n=i } ; next }
  $1==s { print (d?$d:"NA"), (n?$n:"NA") }' "$TSV")"

echo
echo "=============================================="
echo " $SAMPLE : observed vs published"
echo "=============================================="
printf " input pairs        %s\n"        "$IN_N"
printf " kept pairs         %s\n"        "$KEPT"
awk -v k="${KEPT:-0}" -v i="${IN_N:-1}" 'BEGIN{printf " retained           %.1f%%\n", 100*k/i}'
printf " mapped (post)      %s%%\n"      "$OBS_MAP"
printf " mean depth         %-8s (published %s)\n" "$OBS_DEPTH" "${PUB_DEPTH:-NA}"
printf " ambiguous Ns       %-8s (published %s)\n" "$OBS_NS"    "${PUB_NS:-NA}"
printf " variants >=3%%      %s\n"       "$NVAR"
echo
echo " Depth should land within a few percent of published. Ns will differ:"
echo " MIN_DEPTH=$MIN_DEPTH is stricter than the original pipeline, so low-depth"
echo " samples carry more Ns here. That is the point of the threshold, not a bug."
echo "=============================================="
