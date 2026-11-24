#!/usr/bin/env bash
# Assemble the shared teaching bundle from the local build output, in exactly
# the layout the SLURM scripts expect, then print the command to push it.
#
# Run on the machine where you did the build (needs the output of 00-03).
#
#   bash stage_bundle.sh
#
set -euo pipefail

# ---------------------------- EDIT HERE ----------------------------
BUILD=${BUILD:-$HOME/working/viral-workshop}            # local build output
STAGE=${STAGE:-$HOME/working/bundle-staging/viral-genomics}
CLUSTER=${CLUSTER:-umicroche-login}
REMOTE=${REMOTE:-/mnt/roche/resources/viral-genomics}
SCRIPTS=${SCRIPTS:-$HOME/working/reads-to-outbreaks/scripts}
# -------------------------------------------------------------------

say() { printf '\n== %s\n' "$*"; }
miss=0
need() { [ -e "$1" ] || { echo "  MISSING: $1"; miss=$((miss+1)); }; }

rm -rf "$STAGE"
mkdir -p "$STAGE"/{data/fastq,data/reference,data/genomes,scripts} \
         "$STAGE"/precomputed/{00_hostdepleted,01_qc,02_mapping,03_consensus,04_variants,05_align,06_tree}

LABELS=$(cut -f1 "$BUILD/data/samples.tsv" | tail -n +2)

say "raw reads -> data/fastq"
# Deliberately the RAW reads, not the host-depleted ones: step 3 of the workshop
# runs host removal, and the retained fraction is the whole point of that step.
# These are public SRA runs, so no new disclosure. If you ever swap in your own
# clinical data, deplete BEFORE staging.
for L in $LABELS; do
  for R in R1 R2; do
    src="$BUILD/raw/${L}_${R}.fastq.gz"
    need "$src" && cp -n "$src" "$STAGE/data/fastq/" 2>/dev/null || true
  done
done

say "reference -> data/reference"
need "$BUILD/data/reference/MN908947.3.fasta"
cp "$BUILD/data/reference/MN908947.3.fasta"* "$STAGE/data/reference/" 2>/dev/null || true

say "consensus genomes -> data/genomes  (relabelled)"
# alignment.sh falls back to $RESOURCES/data/genomes/<LABEL>.fa, so the headers
# must be the bare label, not ivar's Consensus_A01_threshold_0.75_quality_20.
: > "$STAGE/data/genomes/all_genomes.fa"
for L in $LABELS; do
  src="$BUILD/consensus/${L}.fa"
  if [ -s "$src" ]; then
    { echo ">${L}"; grep -v '^>' "$src"; } > "$STAGE/data/genomes/${L}.fa"
    cat "$STAGE/data/genomes/${L}.fa" >> "$STAGE/data/genomes/all_genomes.fa"
  else
    echo "  MISSING: $src"; miss=$((miss+1))
  fi
done

say "sample table + public reference files -> data/"
# The participant copy must NOT carry the introduction column: that is the answer
# to Practical 4. The full table with assignments stays in the facilitator run
# sheet, not in the shared folder participants can read.
need "$BUILD/data/samples.tsv"
if [ -s "$BUILD/data/samples.tsv" ]; then
  awk -F'\t' 'BEGIN{OFS="\t"; split("introduction clade pango role",d," "); for(k in d) drop[d[k]]=1}
    NR==1 { for(i=1;i<=NF;i++) if($i in drop) skip[i]=1 }
    { out=""; for(i=1;i<=NF;i++) if(!(i in skip)) out=(out==""?$i:out OFS $i); print out }' \
    "$BUILD/data/samples.tsv" > "$STAGE/data/samples.tsv"
  echo "  stripped from participant samples.tsv: introduction, clade, pango, role"
fi
# The CDC TSV also carries the assignments, so it stays out of the bundle too.
need "$BUILD/data/sarscov2-SNF-A.dnd" && cp "$BUILD/data/sarscov2-SNF-A.dnd" "$STAGE/data/" 2>/dev/null || true
echo "  NOTE: sars-cov-2-SNF-A.tsv withheld - it contains the introduction assignments"

say "precomputed output"
cp "$BUILD"/logs/*_riak.log      "$STAGE/precomputed/00_hostdepleted/" 2>/dev/null || echo "  no riak logs"
cp "$BUILD"/qc/*.html            "$STAGE/precomputed/01_qc/"           2>/dev/null || true
cp "$BUILD"/qc/*.flagstat        "$STAGE/precomputed/02_mapping/"      2>/dev/null || true
cp "$BUILD"/bam/*.bam "$BUILD"/bam/*.bai "$STAGE/precomputed/02_mapping/" 2>/dev/null || echo "  no BAMs (large; optional)"
cp "$BUILD"/consensus/*.fa       "$STAGE/precomputed/03_consensus/"     2>/dev/null || true
cp "$BUILD"/depth/*.depth.tsv    "$STAGE/precomputed/03_consensus/"     2>/dev/null || true
cp "$BUILD"/variants/*.tsv       "$STAGE/precomputed/04_variants/"      2>/dev/null || true
cp "$BUILD"/align/all.aln        "$STAGE/precomputed/05_align/"         2>/dev/null || true
cp "$BUILD"/tree/all.treefile "$BUILD"/tree/all.iqtree "$STAGE/precomputed/06_tree/" 2>/dev/null || true

say "SLURM scripts -> scripts/"
if [ -d "$SCRIPTS" ]; then cp "$SCRIPTS"/*.sh "$STAGE/scripts/"; else echo "  MISSING: $SCRIPTS"; miss=$((miss+1)); fi

say "manifest"
( cd "$STAGE" && find . -type f | sort > MANIFEST.txt )
( cd "$STAGE" && shasum -a 256 $(find . -type f ! -name CHECKSUMS.txt ! -name MANIFEST.txt) > CHECKSUMS.txt 2>/dev/null ) || true

say "summary"
printf "fastq files:    %s  (expect 24)\n" "$(ls -1 "$STAGE"/data/fastq/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')"
printf "genome files:   %s  (expect 12 + all_genomes.fa)\n" "$(ls -1 "$STAGE"/data/genomes/A*.fa 2>/dev/null | wc -l | tr -d ' ')"
printf "scripts:        %s\n" "$(ls -1 "$STAGE"/scripts/*.sh 2>/dev/null | wc -l | tr -d ' ')"
printf "total size:     %s\n" "$(du -sh "$STAGE" | cut -f1)"
printf "missing items:  %s\n" "$miss"

cat <<TXT

Push it. macOS 15 ships Apple's openrsync, which segfaults against older remote
rsync ("unexpected end of file", "child exited with status 11"). Check which you
have with: rsync --version

  # GNU rsync 3.x (install with: conda install -c conda-forge rsync)
  rsync -avh --progress "$STAGE/" ${CLUSTER}:${REMOTE}/

  # or, with no working rsync - one shot, no resume, but fast
  cd "$(dirname "$STAGE")" && tar cf - "$(basename "$STAGE")" \
    | ssh ${CLUSTER} 'cd "$(dirname "${REMOTE}")" && tar xvf -'

Then on the cluster, make it read-only to participants:

  chmod -R a+rX,go-w ${REMOTE}

And check one path the scripts actually use:

  ls ${REMOTE}/data/fastq/A01_R1.fastq.gz ${REMOTE}/data/genomes/A01.fa
TXT

[ "$miss" -eq 0 ] || { echo; echo "WARNING: $miss item(s) missing - fix before pushing."; exit 1; }
