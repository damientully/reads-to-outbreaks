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
  NS=""; [ -s results/consensus/${S}.fa ] && \
    NS=$(grep -v '^>' results/consensus/${S}.fa | tr -cd 'Nn' | wc -c | tr -d ' ')

  read -r PD PN <<< "$(awk -F'\t' -v s="$S" '
    NR==1 { for (i=1;i<=NF;i++) { if ($i=="mean_depth") a=i; if ($i=="ambig_ns") b=i }; next }
    $1==s { print (a?$a:"NA"), (b?$b:"NA") }' "$TSV")"
  PCT=$(awk -v k="${KEPT:-0}" -v i="${IN:-0}" 'BEGIN{ if(i>0) printf "%.1f", 100*k/i; else printf "-" }')

  printf "%-5s %10s %10s %8s %8s %9s %9s %8s %8s\n" \
    "$S" "${IN:--}" "${KEPT:--}" "$PCT" "${MAP:--}" "${DEP:--}" "${PD:-NA}" "${NS:--}" "${PN:-NA}"
done

echo

# These notes used to be hardcoded, which meant the footer kept asserting
# "0.5% to 88.7%" and ">98% for every sample" whatever the run actually
# produced. Compute them, and say plainly which steps have not run yet.
awk -F'\t' 'NR>1{print $1}' "$TSV" > /tmp/.rep_labels.$$
MISS_MAP=0; MISS_CONS=0; N=0
MINP=; MAXP=; MINS=; MAXS=; MINMAP=
while read -r S; do
  N=$((N+1))
  R=results/host_removal/${S}.riak.log
  IN=$(awk -F'\t' '/^Input reads file 1/{print $2}' "$R" 2>/dev/null)
  KEPT=$(awk -F'\t' '/^Kept reads 1/{print $2}' "$R" 2>/dev/null)
  if [ -n "${IN:-}" ] && [ "${IN:-0}" -gt 0 ] 2>/dev/null; then
    P=$(awk -v k="${KEPT:-0}" -v i="$IN" 'BEGIN{printf "%.1f", 100*k/i}')
    if [ -z "$MINP" ] || awk -v a="$P" -v b="$MINP" 'BEGIN{exit !(a<b)}'; then MINP=$P; MINS=$S; fi
    if [ -z "$MAXP" ] || awk -v a="$P" -v b="$MAXP" 'BEGIN{exit !(a>b)}'; then MAXP=$P; MAXS=$S; fi
  fi
  M=$(awk '/ mapped \(/{gsub(/[(%]/,"");print $5; exit}' results/mapping/${S}.flagstat 2>/dev/null)
  if [ -z "${M:-}" ]; then MISS_MAP=$((MISS_MAP+1));
  elif [ -z "$MINMAP" ] || awk -v a="$M" -v b="$MINMAP" 'BEGIN{exit !(a<b)}'; then MINMAP=$M; fi
  [ -s results/consensus/${S}.fa ] || MISS_CONS=$((MISS_CONS+1))
done < /tmp/.rep_labels.$$
rm -f /tmp/.rep_labels.$$

[ -n "$MINP" ] && printf "KEPT%%   retention through host depletion: %s%% (%s) to %s%% (%s)\n" \
  "$MINP" "$MINS" "$MAXP" "$MAXS"
if [ "$MISS_MAP" -eq "$N" ]; then
  echo "MAP%    no mapping output found - mapping.sh has not been run"
else
  [ -n "$MINMAP" ] && printf "MAP%%    after depletion, so high for every sample (lowest %s%%) - see step 4b\n" "$MINMAP"
  [ "$MISS_MAP" -gt 0 ] && printf "        %s of %s samples have no mapping output yet\n" "$MISS_MAP" "$N"
fi
if [ "$MISS_CONS" -eq "$N" ]; then
  echo "OBS_NS  blank for every sample - consensus.sh has not been run yet"
else
  echo "OBS_NS  exceeds PUB_NS on low-depth samples: MIN_DEPTH=20 is stricter"
  [ "$MISS_CONS" -gt 0 ] && printf "        %s of %s samples have no consensus yet\n" "$MISS_CONS" "$N"
fi

# Which reads did mapping actually consume? A stale script that maps the
# host-depleted reads directly produces a table that looks entirely normal.
if [ "$MISS_MAP" -lt "$N" ] && [ ! -d results/trimming ]; then
  echo
  echo "NOTE    results/mapping exists but results/trimming does not, so these"
  echo "        depths were not computed on trimmed reads. Check that"
  echo "        scripts/mapping.sh reads from results/trimming, then rerun"
  echo "        trimming.sh and mapping.sh."
fi
