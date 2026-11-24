#!/usr/bin/env bash
# Alignment + ML tree. Includes the reference as an outgroup so the tree can be
# rooted, and shortens tip labels (iVar writes Consensus_<x>_threshold_...).
set -euo pipefail
cd "$(dirname "$0")/.."

REF=data/reference/MN908947.3.fasta
THREADS=${THREADS:-4}
mkdir -p align tree

# 1. collect consensus genomes with clean labels
: > align/consensus_all.fa
for L in $(cut -f1 data/samples.tsv | tail -n +2); do
  [ -s "consensus/${L}.fa" ] || continue
  echo ">${L}" >> align/consensus_all.fa
  grep -v '^>' "consensus/${L}.fa" >> align/consensus_all.fa
done

# 2. add the reference as outgroup
echo ">REF_Wuhan_Hu_1" >> align/consensus_all.fa
grep -v '^>' "$REF" >> align/consensus_all.fa

grep -c '^>' align/consensus_all.fa | xargs echo "sequences in alignment:"

# 3. align
mafft --thread "$THREADS" --auto align/consensus_all.fa > align/all.aln 2> logs/mafft.log

# 4. tree, rooted on the reference
IQ=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
"$IQ" -s align/all.aln -m GTR+G -B 1000 --prefix tree/all -T "$THREADS" \
  -o REF_Wuhan_Hu_1 -redo > logs/iqtree.log 2>&1

echo "tree: tree/all.treefile (rooted on REF_Wuhan_Hu_1)"
echo
echo "NOTE: ultrafast bootstrap values are only credible at >=95."
cat tree/all.treefile
