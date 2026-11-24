#!/usr/bin/env bash
# Per-sample pipeline: host depletion -> QC -> trim -> map -> depth -> consensus -> variants.
# Host depletion now feeds the rest of the pipeline (it did not in the first version).
set -euo pipefail
cd "$(dirname "$0")/.."

REF=data/reference/MN908947.3.fasta
THREADS=${THREADS:-4}
mkdir -p hostdepleted qc trimmed bam depth consensus variants logs

# --- pick a read-QC tool -------------------------------------------------
# A 'fastqc' on PATH is not always FastQC: some environments shadow it with an
# unrelated Python script (symptom: "unrecognized arguments: -o"). Verify, then
# fall back to falco (drop-in, no Java) and finally to fastp's own report.
QC_TOOL="none"
if command -v fastqc >/dev/null 2>&1 && fastqc --version 2>&1 | grep -qi "fastqc"; then
  QC_TOOL="fastqc"
elif command -v falco >/dev/null 2>&1; then
  QC_TOOL="falco"
fi
echo "read-QC tool: $QC_TOOL   ($(command -v fastqc 2>/dev/null || echo 'no fastqc on PATH'))"
if [ "$QC_TOOL" = "none" ]; then
  echo "  -> no working FastQC/falco; relying on fastp HTML reports in qc/"
fi

# Host depletion is the whole point of step 1. If the tool is missing, the most
# likely reason is that the conda environment is not active - so stop, rather
# than quietly producing a bundle with the wrong provenance.
if ! command -v readItAndKeep >/dev/null 2>&1; then
  cat >&2 <<'ERR'

ERROR: readItAndKeep is not on PATH.

  Most likely you have not activated the environment:

      conda activate viral-workshop

  Check with:  command -v readItAndKeep fastqc fastp minimap2 samtools ivar

  To build from raw (undepleted) reads anyway - not recommended, it makes the
  retention figures in the teaching material unreproducible:

      ALLOW_RAW=1 bash scripts/02_build.sh
ERR
  [ "${ALLOW_RAW:-0}" = "1" ] || exit 1
  echo "  ALLOW_RAW=1 set - continuing with raw reads" >&2
fi


for L in $(cut -f1 data/samples.tsv | tail -n +2); do
  echo "=== $L ==="
  R1=raw/${L}_R1.fastq.gz; R2=raw/${L}_R2.fastq.gz
  [ -s "$R1" ] || { echo "  missing $R1 - run 00_fetch.sh"; continue; }

  # 1. host depletion
  if command -v readItAndKeep >/dev/null; then
    # Re-run if the log is missing: the retention counts are teaching material,
    # not a by-product, so depleted FASTQ without its log is incomplete.
    if [ ! -s "hostdepleted/${L}.reads_1.fastq.gz" ] || [ ! -s "logs/${L}_riak.log" ]; then
      readItAndKeep --tech illumina --ref_fasta "$REF" \
        --reads1 "$R1" --reads2 "$R2" --outprefix "hostdepleted/$L" \
        > "logs/${L}_riak.log" 2>&1
    else
      echo "  [skip] host depletion already done"
    fi
    IN1=hostdepleted/${L}.reads_1.fastq.gz
    IN2=hostdepleted/${L}.reads_2.fastq.gz
  else
    echo "  readItAndKeep not found - using raw reads (NOT the intended path)"
    IN1="$R1"; IN2="$R2"
  fi

  # 2. QC on the reads we will actually analyse
  case "$QC_TOOL" in
    fastqc) fastqc -t "$THREADS" --outdir qc/ "$IN1" "$IN2" > "logs/${L}_fastqc.log" 2>&1 \
              || echo "  fastqc FAILED - see logs/${L}_fastqc.log" ;;
    falco)  falco -o "qc/${L}_falco" "$IN1" "$IN2" > "logs/${L}_fastqc.log" 2>&1 \
              || echo "  falco FAILED - see logs/${L}_fastqc.log" ;;
    none)   : ;;
  esac

  # 3. trim
  fastp -i "$IN1" -I "$IN2" \
    -o trimmed/${L}_R1.fastq.gz -O trimmed/${L}_R2.fastq.gz \
    --json qc/${L}_fastp.json --html qc/${L}_fastp.html \
    --thread "$THREADS" --cut_tail --cut_tail_mean_quality 20 --length_required 50 \
    > /dev/null 2> "logs/${L}_fastp.log"

  # 4. map
  minimap2 -t "$THREADS" -ax sr "$REF" trimmed/${L}_R1.fastq.gz trimmed/${L}_R2.fastq.gz \
      2> "logs/${L}_mm2.log" \
    | samtools sort -@ "$THREADS" -o bam/${L}.sorted.bam
  samtools index bam/${L}.sorted.bam
  samtools flagstat bam/${L}.sorted.bam > qc/${L}.flagstat

  # 5. depth, consensus, variants
  samtools depth -a -d 0 bam/${L}.sorted.bam > depth/${L}.depth.tsv
  samtools mpileup -aa -A -d 0 -Q 0 --reference "$REF" bam/${L}.sorted.bam 2>/dev/null \
    | ivar consensus -p consensus/${L} -t 0.75 -m 20 -n N > "logs/${L}_cons.log" 2>&1
  samtools mpileup -aa -A -d 0 -B -Q 0 --reference "$REF" bam/${L}.sorted.bam 2>/dev/null \
    | ivar variants -p variants/${L} -r "$REF" -t 0.03 -m 20 > "logs/${L}_var.log" 2>&1
done

echo
echo "=== observed vs published ==="
printf "%-5s %-13s %7s %9s %9s %7s %7s %9s\n" \
  LABEL SAMPLE KEPT% MAP%_POST PUB_DEPTH OBS_DEP PUB_NS OBS_NS
while IFS=$'\t' read -r label sample sra gb intro clade pango depth ns pct role; do
  [ "$label" = "label" ] && continue
  kept="-"; mapp="-"; obsd="-"; obsn="-"
  f=logs/${label}_riak.log
  if [ -s "$f" ]; then
    i=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$f")
    k=$(awk -F'\t' '/^Kept reads 1/{print $2}' "$f")
    [ -n "${i:-}" ] && [ -n "${k:-}" ] && kept=$(awk -v k="$k" -v i="$i" 'BEGIN{printf "%.1f", k/i*100}')
  fi
  [ -s qc/${label}.flagstat ] && mapp=$(awk '/ mapped \(/{gsub(/[(%]/,"");print $5; exit}' qc/${label}.flagstat)
  [ -s depth/${label}.depth.tsv ] && obsd=$(awk '{s+=$3} END {printf "%.0f", s/NR}' depth/${label}.depth.tsv)
  [ -s consensus/${label}.fa ] && obsn=$(grep -v '>' consensus/${label}.fa | tr -cd 'Nn' | wc -c | tr -d ' ')
  printf "%-5s %-13s %7s %9s %9s %7s %7s %9s\n" \
    "$label" "$sample" "$kept" "$mapp" "$depth" "$obsd" "$ns" "$obsn"
done < data/samples.tsv
echo
echo "KEPT%     = read pairs surviving host depletion (low = lots of non-viral sequence)"
echo "MAP%_POST = mapping rate AFTER depletion, so it should be high for every sample"
echo
echo "provenance: read-QC=$QC_TOOL  host_depletion=$(command -v readItAndKeep >/dev/null 2>&1 && echo readItAndKeep || echo NONE)"
if ! command -v readItAndKeep >/dev/null 2>&1; then
  echo "WARNING: this build used RAW reads. Do not stage it as the teaching bundle."
fi
