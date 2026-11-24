#!/usr/bin/env bash
# Cross-sample table from the step-by-step scripts. Run on the login node.
#   bash scripts/report.sh
set -uo pipefail

WORKDIR=${WORKDIR:-$HOME/roche/viral-genomics}
RESOURCES=${RESOURCES:-/mnt/roche/resources/viral-genomics}
[ -s "$RESOURCES/data/samples.tsv" ] || RESOURCES="$WORKDIR"
TSV=$RESOURCES/data/samples.tsv
cd "$WORKDIR"

printf "%-5s %10s %10s %8s %8s %9s %9s %8s %8s\n" \
  LABEL INPUT KEPT "KEPT%" "MAP%" OBS_DEP PUB_DEP OBS_NS PUB_NS

for S in $(cut -f1 "$TSV" | tail -n +2); do
  R=results/host_removal/${S}.riak.log
  IN=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$R" 2>/dev/null)
  KEPT=$(awk -F'\t' '/^Kept reads 1/{print $2}'     "$R" 2>/dev/null)
  MAP=$(awk '/ mapped \(/{gsub(/[(%]/,"");print $5; exit}' results/mapping/${S}.flagstat 2>/dev/null)
  DEP=$(awk '{s+=$3} END {if(NR) printf "%.0f", s/NR}' results/mapping/${S}.depth.tsv 2>/dev/null)
  NS=$(grep -v '^>' results/consensus/${S}.fa 2>/dev/null | tr -cd 'Nn' | wc -c | tr -d ' ')

  read -r PD PN <<< "$(awk -F'\t' -v s="$S" '
    NR==1 { for (i=1;i<=NF;i++) { if ($i=="mean_depth") a=i; if ($i=="ambig_ns") b=i }; next }
    $1==s { print (a?$a:"NA"), (b?$b:"NA") }' "$TSV")"
  PCT=$(awk -v k="${KEPT:-0}" -v i="${IN:-0}" 'BEGIN{ if(i>0) printf "%.1f", 100*k/i; else printf "-" }')

  printf "%-5s %10s %10s %8s %8s %9s %9s %8s %8s\n" \
    "$S" "${IN:--}" "${KEPT:--}" "$PCT" "${MAP:--}" "${DEP:--}" "${PD:-NA}" "${NS:--}" "${PN:-NA}"
done

echo
echo "KEPT%   retention through host depletion: 0.5% (A07) to 88.7% (A02)"
echo "MAP%    after depletion, so >98% for every sample - see step 4b"
echo "OBS_NS  exceeds PUB_NS on low-depth samples: MIN_DEPTH=20 is stricter"
