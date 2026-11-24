#!/usr/bin/env bash
# Collect the observed-vs-published figures for every sample processed by
# one_sample.sh into a single table. Run on the login node, not via sbatch.
#
#   bash scripts/summarise.sh
#
set -uo pipefail

WORKDIR=${WORKDIR:-$HOME/roche/viral-genomics}
RESOURCES=${RESOURCES:-/mnt/roche/resources/viral-genomics}
[ -s "$RESOURCES/data/samples.tsv" ] || RESOURCES="$WORKDIR"
TSV=$RESOURCES/data/samples.tsv

cd "$WORKDIR"
BASE=results/one_sample
[ -d "$BASE" ] || { echo "no $BASE - run one_sample.sh first"; exit 1; }

printf "%-5s %10s %10s %8s %8s %9s %9s %8s %8s\n" \
  LABEL INPUT KEPT "KEPT%" "MAP%" OBS_DEP PUB_DEP OBS_NS PUB_NS

for S in $(cut -f1 "$TSV" | tail -n +2); do
  d=$BASE/$S
  [ -d "$d" ] || continue

  IN=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$d/${S}.riak.log" 2>/dev/null)
  KEPT=$(awk -F'\t' '/^Kept reads 1/{print $2}'     "$d/${S}.riak.log" 2>/dev/null)
  MAP=$(awk '/ mapped \(/{gsub(/[(%]/,"");print $5; exit}' "$d/qc/${S}.flagstat" 2>/dev/null)
  DEP=$(awk '{s+=$3} END {if(NR) printf "%.0f", s/NR}' "$d/depth/${S}.depth.tsv" 2>/dev/null)
  NS=$(grep -v '^>' "$d/consensus/${S}.fa" 2>/dev/null | tr -cd 'Nn' | wc -c | tr -d ' ')

  read -r PD PN <<< "$(awk -F'\t' -v s="$S" '
    NR==1 { for (i=1;i<=NF;i++) { if ($i=="mean_depth") a=i; if ($i=="ambig_ns") b=i }; next }
    $1==s { print (a?$a:"NA"), (b?$b:"NA") }' "$TSV")"

  PCT=$(awk -v k="${KEPT:-0}" -v i="${IN:-0}" 'BEGIN{ if(i>0) printf "%.1f", 100*k/i; else printf "-" }')

  printf "%-5s %10s %10s %8s %8s %9s %9s %8s %8s\n" \
    "$S" "${IN:--}" "${KEPT:--}" "$PCT" "${MAP:--}" "${DEP:--}" "${PD:-NA}" "${NS:--}" "${PN:-NA}"
done

echo
echo "KEPT%   retention through host depletion; expect 0.5% (A07) to 88.7% (A02)"
echo "MAP%    mapping rate AFTER depletion; expect >98% for every sample"
echo "OBS_DEP should land within a few percent of PUB_DEP"
echo "OBS_NS  will exceed PUB_NS on low-depth samples: MIN_DEPTH=20 is stricter"
echo "        than the original pipeline. Expected, not a fault."
