# Building the teaching bundle

Facilitator scripts. They turn public accessions into the bundle staged on the
cluster. Run on a machine with internet and sra-tools; the cluster has neither.

```
conda env create -f environment.yml && conda activate viral-workshop
bash 00_fetch.sh        # the 12 SRA runs
bash 01_reference.sh    # reference + published consensus genomes from GenBank
bash 02_build.sh        # host depletion, QC, mapping, consensus, variants
bash 03_tree.sh         # alignment + tree, rooted on the reference
bash stage_bundle.sh    # assemble and verify, then rsync to the cluster
```

`02_build.sh` refuses to run without `readItAndKeep` on PATH, and prints
observed depth and N counts beside the published values so you can confirm the
samples behave as documented.

`check_hpc_modules.sh` runs on the cluster and reports which tools resolve.

Design decisions, sample selection, provenance and the findings from each test
run are recorded in the setup specification, not here.
