#!/usr/bin/env bash
# Download the 12 SRA runs. Run this on macOS where sra-tools is installed.
set -euo pipefail
cd "$(dirname "$0")/.."

command -v prefetch     >/dev/null || { echo "ERROR: prefetch not found (conda install -c bioconda sra-tools)"; exit 1; }
command -v fasterq-dump >/dev/null || { echo "ERROR: fasterq-dump not found"; exit 1; }

mkdir -p raw sra
while IFS=$'\t' read -r label sample sra rest; do
  [ "$label" = "label" ] && continue
  if [ -s "raw/${label}_R1.fastq.gz" ]; then echo "[skip] $label"; continue; fi
  echo "[fetch] $label  $sra"
  prefetch --max-size 50G -O sra "$sra"
  fasterq-dump --split-files --skip-technical -O raw "sra/$sra"
  mv "raw/${sra}_1.fastq" "raw/${label}_R1.fastq"
  mv "raw/${sra}_2.fastq" "raw/${label}_R2.fastq"
  gzip -f "raw/${label}_R1.fastq" "raw/${label}_R2.fastq"
done < data/samples.tsv

echo; echo "=== sizes ==="; ls -lh raw/*.fastq.gz | awk '{print $5, $9}'
