#!/usr/bin/env bash
# Run this ON THE CLUSTER (login node is fine) to see what is available.
# Usage: bash check_hpc_modules.sh
# Reports, per tool: whether a module matching the name exists, and whether
# the executable is on PATH after attempting to load it.

ESSENTIAL="samtools minimap2 ivar fastp fastqc mafft read_keep iqtree2"
TREE="iqtree2 iqtree iqtree3"
OPTIONAL="seqkit multiqc"

have_module_cmd() { type module >/dev/null 2>&1; }

probe () {
  local name="$1" group="$2" mod="-" exe="-"
  if have_module_cmd; then
    if module avail "$name" 2>&1 | grep -qi "$name"; then mod="yes"; else mod="no"; fi
    module load "$name" >/dev/null 2>&1 || true
  fi
  for cand in "$name" "${name//-/}" "${name,,}"; do
    if command -v "$cand" >/dev/null 2>&1; then exe="$(command -v "$cand")"; break; fi
  done
  printf "%-12s %-18s module:%-5s exe:%s\n" "$group" "$name" "$mod" "$exe"
}

echo "=== module system ==="
if have_module_cmd; then module --version 2>&1 | head -2; else echo "no 'module' command found"; fi
echo
echo "=== ESSENTIAL (session breaks without these) ==="
for t in $ESSENTIAL; do probe "$t" ESSENTIAL; done
echo
echo "=== TREE (need ONE of these) ==="
for t in $TREE; do probe "$t" TREE; done
echo
echo "=== OPTIONAL (degrade gracefully) ==="
for t in $OPTIONAL; do probe "$t" OPTIONAL; done
echo
echo "=== apptainer/singularity fallback ==="
for t in apptainer singularity; do printf "%-12s %s\n" "$t" "$(command -v $t || echo '-')"; done
echo
echo "=== SLURM ==="
for t in sinfo squeue salloc sbatch sacct; do printf "%-10s %s\n" "$t" "$(command -v $t || echo '-')"; done
echo
echo "Partitions (if sinfo present):"
command -v sinfo >/dev/null && sinfo -o "%P %a %l %D %c %m" 2>/dev/null | head -10 || echo "  n/a"
