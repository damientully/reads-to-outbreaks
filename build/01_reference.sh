#!/usr/bin/env bash
# Fetch reference genome + the 12 published consensus genomes from GenBank.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/reference data/genomes

if command -v efetch >/dev/null; then
  efetch -db nuccore -id NC_045512.2 -format fasta > data/reference/MN908947.3.fasta
  cut -f4 data/samples.tsv | tail -n +2 | tr '\n' ',' | sed 's/,$//' \
    | xargs -I{} efetch -db nuccore -id {} -format fasta > data/genomes/all_genomes.fa
else
  echo "entrez-direct not found; falling back to direct download"
  curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" \
    > data/reference/MN908947.3.fasta
  ACC=$(cut -f4 data/samples.tsv | tail -n +2 | tr '\n' ',' | sed 's/,$//')
  curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" \
    > data/genomes/all_genomes.fa
fi

command -v samtools >/dev/null && samtools faidx data/reference/MN908947.3.fasta
grep -c '>' data/genomes/all_genomes.fa | xargs echo "genomes fetched:"
echo "Collection dates (for metadata.tsv) are in the GenBank records:"
echo "  efetch -db nuccore -id MT520529.1 -format gb | grep collection_date"
