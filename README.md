# From Reads to Outbreaks

A two-hour hands-on workshop in viral pathogen genomics: twelve SARS-CoV-2
genomes from a real skilled nursing facility outbreak, taken from raw sequencing
reads to an epidemiological conclusion on an HPC cluster under SLURM.

**[Read the workshop walkthrough](docs/index.md)**

## What is here

```
docs/index.md        the participant walkthrough, steps 1-6
scripts/             SLURM scripts, one per step
data/samples.tsv     the twelve samples: SRA runs, introductions, published QC metrics
data/                supporting public files (dataset table, reference tree)
```

Sequencing data is not in this repository. The facilitator stages it into a
shared folder on the cluster; see `docs/index.md`.

## The dataset

Twelve samples selected from the Boston skilled nursing facility outbreak
(SNF-A) in the [CDC SARS-CoV-2 benchmark
datasets](https://github.com/CDCgov/datasets-sars-cov-2), which assign each
sample to one of three documented introductions. Published QC metrics come from
[Bhatt et al., PeerJ 2022](https://peerj.com/articles/13821/); the underlying
genomes are from [Lemieux et al., Science
2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7857412).

The selection is deliberate rather than representative: it spans all three
introductions, both clades, three Pango lineages, and a 200-fold range in
sequencing depth, including two genomes that should not be used for fine-grained
inference.

## Tools

Everything runs from cluster modules: fastqc, multiqc, fastp, read_keep,
minimap2, samtools, ivar, mafft, iqtree2 and seqkit.

## Credits

Dataset and benchmark metrics are the work of the cited authors. Workshop
structure follows the conventions of the
[UMIC Roche metagenomics training](https://github.com/AlfredUg/metagenomics-umic-roche).
