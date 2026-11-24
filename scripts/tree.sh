#!/bin/bash
# Logs are written beside wherever you submit from, as vg_<step>_<jobid>.out/.err.
# They deliberately do NOT go in a subdirectory: SLURM opens the output file at
# job launch, and a missing directory makes the job fail or hang silently on
# this cluster.
#SBATCH --job-name=vg_tree
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=vg_tree_%j.out
#SBATCH --error=vg_tree_%j.err

# ---------------------------- EDIT HERE ----------------------------
WORKDIR=$HOME/roche/viral-genomics
OUTGROUP=REF_Wuhan_Hu_1
# -------------------------------------------------------------------

# ---------------------------- ENVIRONMENT --------------------------
module load iqtree2
# -------------------------------------------------------------------

set -euo pipefail
cd "$WORKDIR"
mkdir -p results/tree logs

ALN=results/alignment/all.aln
OUT=results/tree

IQ=$(command -v iqtree2 || command -v iqtree3 || command -v iqtree || true)

if [ -n "$IQ" ]; then
  # 1000 ultrafast bootstrap replicates, rooted on the reference.
  "$IQ" -s "$ALN" -m GTR+G -B 1000 --prefix "$OUT/all" \
    -T "${SLURM_CPUS_PER_TASK:-4}" -o "$OUTGROUP" -redo
  TREE=$OUT/all.treefile
  cat <<'TXT'

Support values are ULTRAFAST BOOTSTRAP, on a 0-100 scale.
Treat >=95 as strong. Below that the node is unresolved, not weakly resolved -
and a figure that draws an unresolved node as a clean split misrepresents it.
TXT
else
  FT=$(command -v FastTree || command -v fasttree)
  "$FT" -gtr -nt -gamma "$ALN" > "$OUT/all.nwk" 2> logs/fasttree.log
  TREE=$OUT/all.nwk
  cat <<'TXT'

FALLBACK: FastTree was used because iqtree2 was not found.
Support values are SH-LIKE LOCAL SUPPORTS on a 0-1 scale. These are NOT
bootstrap replicates and are not comparable to 0-100 bootstrap values.
Treat >=0.95 as strong, 0.70-0.95 as suggestive, below 0.70 as unresolved.
TXT
fi

echo
echo "tree: $TREE"
cat "$TREE"
