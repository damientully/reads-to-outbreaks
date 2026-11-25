#!/usr/bin/env bash
# One MultiQC report across ALL samples. Run on the login node AFTER the array
# finishes - not with sbatch, and not inside an array.
#
#   bash scripts/multiqc.sh
#
# Why this is separate: qc.sh is a job array over twelve samples. A multiqc
# inside it would scan results/qc while the other eleven tasks were still
# writing there, so every report would cover a different, arbitrary subset.
set -uo pipefail

WORKDIR=${WORKDIR:-$HOME/roche/viral-genomics}
cd "$WORKDIR" || { echo "ERROR: no $WORKDIR" >&2; exit 1; }

command -v multiqc >/dev/null || module load multiqc 2>/dev/null || true
command -v multiqc >/dev/null || { echo "ERROR: multiqc not on PATH. Try: module load multiqc" >&2; exit 1; }

made=0
for STEP in qc trimming; do
  DIR=results/$STEP
  [ -d "$DIR" ] || continue
  N=$(ls -1 "$DIR"/*_fastqc.zip 2>/dev/null | wc -l | tr -d ' ')
  if [ "$N" -eq 0 ]; then
    echo "skipping $DIR: no FastQC output yet"
    continue
  fi
  EXPECT=24   # 12 samples x R1+R2
  [ "$N" -lt "$EXPECT" ] && \
    echo "NOTE: $DIR has $N of $EXPECT FastQC files - has the array finished? (squeue -u \$USER)"
  multiqc --force -o "$DIR" -n all_multiqc "$DIR" >/dev/null 2>&1 \
    && { echo "wrote $DIR/all_multiqc.html  ($N FastQC files)"; made=$((made+1)); } \
    || echo "multiqc failed on $DIR"
done

[ "$made" -gt 0 ] || { echo "nothing to report on - run scripts/qc.sh first" >&2; exit 1; }
cat <<'TXT'

One report, every sample. The General Statistics table at the top is the point:
read counts and quality side by side for all twelve, which is what tells you
these samples are not equivalent. Copy it over to look at it:

  scp -O <user>@umicroche-login:$WORKDIR/results/qc/all_multiqc.html .
TXT
