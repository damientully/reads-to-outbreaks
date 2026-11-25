# Pathogen Genomics

**From Reads to Outbreaks** — practical viral pathogen genomics using NGS data.
Twelve SARS-CoV-2 genomes from a real skilled nursing facility outbreak, taken
from raw sequencing reads to an epidemiological conclusion.

By the end you will have built a consensus genome from raw reads, found the
variation in it, placed twelve genomes on a tree, and decided whether an outbreak
was one transmission chain or several introductions.

Everything runs on the cluster through SLURM. **Every step has a precomputed
answer** in `precomputed/`, so if a job fails or runs long, copy the output and
keep moving. That is expected, not a last resort.

Questions marked **▸** are for you to answer as you go.

## What this dataset is

These are real genomes from real patients. In spring 2020, an outbreak tore
through a 142-bed skilled nursing facility in Massachusetts. The mean resident
age was 83. The facility has three parts — short-term care, long-term care and a
memory unit — and cases appeared across them over several weeks. The outbreak
was sequenced, analysed and published, and it caused significant mortality in a
population that had almost no defence against it.

Twelve of those genomes are what you will work with. They were produced by
unbiased metagenomic sequencing on Illumina at 101 bp, which is why your first
step is removing human reads rather than trimming primers.

The twelve were not chosen to be representative. They were chosen to make
specific things visible:

- Both circulating clades, and three different Pango lineages
- A 200-fold range in sequencing depth, from 30× to over 7,000×
- A 170-fold range in how much of each run was actually viral
- Two genomes that should not be used for fine-grained inference, and one that
  looks deceptively good until you count its comparable positions
- At least one genotype shared by samples from different units

**The answer exists.** The published dataset records which introduction each
sample belongs to, and your facilitator has it. You are reconstructing a
conclusion that is already known, which means you can find out afterwards how
close you came — and, more usefully, whether the evidence you had actually
justified the confidence you felt.

### Provenance

| | |
| --- | --- |
| Outbreak and genomes | [Lemieux et al., *Science* 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7857412) |
| Benchmark dataset and introduction assignments | [CDC SARS-CoV-2 benchmark datasets](https://github.com/CDCgov/datasets-sars-cov-2) |
| Published QC metrics | [PeerJ 2022](https://peerj.com/articles/13821/) |
| Reference genome | Wuhan-Hu-1, MN908947.3 |

Every accession is in `data/samples.tsv`, so anything you find today can be
traced back to the public record.

## Layout

```
/mnt/roche/resources/viral-genomics/      shared, read-only
  data/fastq/            A01..A12 reads
  data/reference/        MN908947.3.fasta
  data/genomes/          consensus genomes (fallback for steps 5-6)
  data/samples.tsv       sample table: accessions and published metrics
  precomputed/           output of every step
  scripts/               the scripts below

~/roche/viral-genomics/                   yours
  scripts/               your copies
  vg_<step>_<jobid>.out  job logs, written where you submit from
  results/               one subfolder per step
```

## The seven steps at a glance

Each is a job array over all twelve samples. Submit from `$WORKDIR`, and wait for
each array to clear `squeue -u $USER` before the next — they chain through
`results/`.

```
cd $WORKDIR
sbatch scripts/qc.sh            # 2.  read quality control
sbatch scripts/host_removal.sh  # 3.  remove human reads
sbatch scripts/trimming.sh      # 4a. trim adapters and low-quality bases
sbatch scripts/mapping.sh       # 4b. map to the reference
sbatch scripts/consensus.sh     # 4c. consensus genome and variants
sbatch scripts/alignment.sh     # 5a. align all twelve plus the reference
sbatch scripts/tree.sh          # 5b. maximum likelihood tree

bash scripts/report.sh          # cross-sample table, any time after 4c
```

The sections below explain what each step does and what to look for. Steps 4a to
4c are sub-steps of section 4, and 5a to 5b of section 5.

## 1. Set up your workspace

```
export WORKDIR=$HOME/roche/viral-genomics
export RESOURCES=/mnt/roche/resources/viral-genomics

mkdir -p $WORKDIR/scripts
cp $RESOURCES/scripts/*.sh $WORKDIR/scripts/
cd $WORKDIR
```

Every tool comes from modules, including `read_keep/0.3.0` and `iqtree2/2.4.0`.

Each step runs as a **job array over all twelve samples** — one task per sample,
so twelve run side by side and the whole step finishes in about the time one
sample takes. Comparing samples is the point of several steps, and you cannot do
that one sample at a time.

Every script has an `EDIT HERE` block at the top. **You do not normally change
anything in it** — `sbatch scripts/qc.sh` already runs all twelve, because
`--array=1-12` is set inside the script.

The `SAMPLE=A01` line in that block is only a fallback for running one sample on
its own. To do that, override the array on the command line:

```
sbatch --array=1  scripts/qc.sh     # A01 only
sbatch --array=7  scripts/qc.sh     # A07 only
sbatch --array=1-3 scripts/qc.sh    # A01 to A03
```

The scripts create their own `results/` subdirectories. Job logs are written
beside you as `vg_<step>_<jobid>.out` and `.err`, so **always submit from
`$WORKDIR`** — every path in the scripts is relative to the directory you submit
from.

Submit with `sbatch`, watch your jobs with `squeue -u $USER`, and read failures
with:

```
sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed
```

The cluster runs SLURM 16.05, which predates some options you may have seen
elsewhere. `squeue --me` does not exist here; use `squeue -u $USER`.

### The samples

`A01`–`A12` come from one outbreak in a 142-bed skilled nursing facility,
sequenced by unbiased metagenomic sequencing on Illumina at 101 bp. The facility
has short-term care, long-term care and a memory unit. `data/samples.tsv` maps
each label to its real sample ID and SRA run, so anything you find can be traced
back to the public record.

## 2. Quality control

```
cd $WORKDIR          # always submit from here
sbatch scripts/qc.sh
```

> This walkthrough is also on the cluster at
> `$RESOURCES/docs/index.md`, so you can read it without leaving the terminal:
> `less $RESOURCES/docs/index.md`.

> **Always submit from `$WORKDIR`, not from inside `scripts/`.** Every path in
> the scripts is relative to the directory you submit from, including where the
> results land. Logs appear beside you as `vg_qc_<jobid>.out` and
> `vg_qc_<jobid>.err`.

The HTML report cannot be opened on a compute node. Either copy it to your own
machine or read the text summary in place:

```
scp <user>@umicroche-login:$WORKDIR/results/qc/A01_R1_fastqc.html .
# or
unzip -p results/qc/A01_R1_fastqc.zip A01_R1_fastqc/summary.txt
```

**What you should see:** high quality across most of the read, declining toward
the 3' end. Alarming duplication levels. A failed "per base sequence content"
plot.

**What it means:** FastQC is judging your data against assumptions that do not
hold here. A tool's warning is relative to a model of what normal looks like, and
you need to know whether that model applies before you act on it.

**Note:** FastQC 0.12 dropped the short `-o` flag. The scripts use `--outdir`. If
you type it yourself and get `unrecognized arguments: -o`, that is why.

**▸ Q1.** Open the report and find these four sections in the left-hand panel:
*Per base sequence quality*, *Per base sequence content*, *Sequence Duplication
Levels* and *Adapter Content*. Some will be marked with a red cross. Which of
those failures would actually concern you in this data, and which are expected
to look bad here?

(FastQC calls its report sections "modules", which is unrelated to the cluster's
`module load`. This walkthrough says "sections" to keep them apart.)

## 3. Host removal

```
sbatch scripts/host_removal.sh
```

These are clinical nasopharyngeal samples sequenced without targeted enrichment,
so most of what came off the machine is patient, not virus. ReadItAndKeep keeps
only reads matching the SARS-CoV-2 reference.

**What you should see:** in `results/host_removal/A01.summary.txt`, a kept-read
count far below the input. For `A01`, 487,837 read pairs survive out of
1,442,994 — about 34%.

Now look at all twelve at once:

```
for f in results/host_removal/*.riak.log; do
  s=$(basename $f .riak.log)
  i=$(awk -F'\t' '/^Input reads file 1/{print $2}' $f)
  k=$(awk -F'\t' '/^Kept reads 1/{print $2}' $f)
  awk -v s=$s -v i=$i -v k=$k 'BEGIN{printf "%-5s %10s %10s %6.1f%%\n", s, i, k, 100*k/i}'
done
```


| Sample | Retained | | Sample | Retained |
| --- | --- | --- | --- | --- |
| A01 | 33.8% | | A07 | **0.5%** |
| A02 | **88.7%** | | A08 | 4.0% |
| A03 | 2.3% | | A09 | 1.2% |
| A04 | 4.3% | | A10 | 12.8% |
| A05 | 4.0% | | A11 | 8.4% |
| A06 | 24.8% | | A12 | 9.8% |

From 0.5% to 88.7% — a 170-fold range, in one study, one laboratory, one
protocol. For `A07`, 4,545 useful read pairs came out of 854,605 sequenced.

**▸ Q2.** Take `A07`: 4,545 kept pairs. Multiply by 2 (both mates) and by 101 bp,
then divide by the genome length, 29,903 bp. Compare your answer with the
published mean depth in `data/samples.tsv`. What does that tell you about where
sequencing depth comes from?

**▸ Q3.** Now try the same arithmetic for `A02`, which kept 1,820,703 pairs. It
does not agree with the published figure, and the difference is large. Give two
explanations that would produce a published depth *lower* than the reads imply,
and say how you would test each one.

**What it means.** Two separate things.

Host sequence from a clinical sample is identifiable patient genetic information.
Removing it before data are shared is a governance requirement, not housekeeping,
and nobody downstream can do it for you.

And the fraction of a sequencing run that was useful is a property of the sample —
how much virus was in that swab — which you cannot know until you measure it. Two
samples on the same run can differ by two orders of magnitude in how much data
they actually gave you.

**▸ Q4.** Your laboratory wants to share raw reads from a clinical metagenomic
run with a collaborator abroad. What has to happen first, and what would you want
documented?

## 4. From reads to a consensus genome

### 4a. Trimming

```
sbatch scripts/trimming.sh
```

Note what we are *not* doing: trimming primer sequences. These are metagenomic
data, so there are no primers. If they were amplicon data there would be an
`ivar trim` step here, and skipping it can produce confidently wrong genotypes
where a variant sits under a primer binding site.

### 4b. Mapping

```
sbatch scripts/mapping.sh
```

**What you should see:** a very high mapping rate, and in
`results/mapping/A01.depth_summary.txt`, the mean depth and the number of
positions below 10×.

**What it means:** that mapping rate looks reassuring and is the wrong thing to
be reassured by. You are mapping the reads that survived step 3, which were
selected for matching SARS-CoV-2. Of course they map. For `A01` the *raw* reads
map at 37.65%.

A quality metric can be flattered by a filtering decision made upstream of it. If
you report a mapping rate, report what it was computed on.

All twelve at once, with the published figures alongside:

```
bash scripts/report.sh
```

Depth should land within a few percent of the published value for most samples.
One does not, and finding it is Q3.

### 4c. Consensus and variants

```
sbatch scripts/consensus.sh
```

**What it means:** two numbers in the `EDIT HERE` block decide what your genome
says. `MIN_DEPTH=20` means any position with fewer than 20 reads becomes an N
rather than a guess. `MIN_FREQ=0.75` means a base is only called when 75% of
reads agree. Lower the depth cutoff for a more complete genome that is less
reliable; raise it for a genome full of Ns that you can trust. There is no
setting that avoids the trade-off, only settings that hide it.

`results/variants/A01.tsv` holds variants down to 3% frequency — far below the
75% needed to change the consensus. Everything between those thresholds was in
the sample and is absent from the FASTA you are about to put on a tree.

```
cut -f2,3,4,11,12,13 results/variants/A01.tsv | column -t | head -30
```

**▸ Q5.** Your consensus has an N at a position. Name three different reasons
that could have happened. Which can you distinguish from the FASTA alone?

**▸ Q6.** You see a variant at 40% frequency. List the possible explanations.
How would you tell them apart?

**▸ Q7.** `A01` and `A07` were sequenced identically. One yields a near-complete
genome and one does not. Where in your analysis would that difference first have
become visible?

### Compare the variant tables

You have a variant table per sample. Compare two from different parts of the
tree:

```
cut -f2,3,4,11 results/variants/A01.tsv | sort -k1,1n > /tmp/a01.txt
cut -f2,3,4,11 results/variants/A05.tsv | sort -k1,1n > /tmp/a05.txt
comm -3 /tmp/a01.txt /tmp/a05.txt | head -20
```

**What it means:** the positions where two samples differ are the entire basis of
everything you will do in steps 5 and 6. There are not many of them. Before you
build a tree, it is worth seeing how thin the evidence actually is — a handful of
positions out of thirty thousand, some of which sit in regions where one of the
two samples has no coverage at all.

## 5. Phylogenetics

### 5a. Alignment

```
sbatch scripts/alignment.sh
```

This collects your consensus genomes, falling back to the shared copies for
samples you have not processed, relabels them, and appends the reference as an
outgroup.

**▸ Q8.** Open the alignment and find a column that is N in some sequences.
Should that column contribute to the tree?

### 5b. Tree

```
sbatch scripts/tree.sh
```

IQ-TREE 2 with 1,000 ultrafast bootstrap replicates, rooted on the reference.
Support values are on a **0–100 scale**. Treat 95 and above as strong; below
that the node is unresolved, not weakly resolved.

Copy the tree to your own machine and open it in FigTree, or paste the Newick
into [icytree.org](https://icytree.org):

```
scp <user>@umicroche-login:$WORKDIR/results/tree/all.treefile .
```

Then, in order:

1. **Confirm it is rooted** on `REF_Wuhan_Hu_1`. An unrooted tree displayed as
   though rooted invites you to read ancestry that is not there.
2. **Read the scale bar.** Branch lengths mean something only relative to it.
3. **Turn on support values.** Then see which groupings survive them.
4. **Find the polytomies** — nodes where three or more lineages emerge together.
   These are where the data could not resolve the order of events.

**What it means:** with twelve closely related genomes you will have short
branches, several identical sequences, and low support on some internal nodes.
None of that is a failure. It is what limited genetic diversity looks like, and
it is the normal situation in an outbreak investigation.

The mistake to avoid is treating an unresolved node as evidence of a close
relationship. A polytomy says the data cannot order these events; it does not say
they happened simultaneously.

Ask the facilitator to show you the same alignment analysed **without** the
reference. The groupings change. Rooting is not a cosmetic choice.

**▸ Q9.** Two genomes in your tree are identical. Name at least three distinct
scenarios consistent with that — then check how many positions were actually
comparable between them before you commit to any of the three.

## 6. Genomic epidemiology

No new commands. You have a tree, and the metadata:

```
column -t $RESOURCES/data/samples.tsv
```

That table gives you, per genome: label, sample ID, SRA run, collection date,
unit, and the published depth and ambiguity counts. It does not tell you which
introduction each sample belongs to — working that out is the exercise.

### The situation

Twelve residents and staff at the facility tested positive over several weeks,
across three units. Management wants to know whether infection control failed
inside the building or whether the virus was repeatedly introduced from outside,
because the two answers demand different responses.

### Your task

Annotate your tree with dates and units — by hand on a printout is entirely
acceptable and often faster than configuring a viewer. Then answer:

1. How many distinct introductions does the genomic evidence support?
2. Which cases, if any, are consistent with transmission within the facility?
3. Which pairs cannot be ordered, and why?
4. What would you tell the infection control team, and what would you decline to
   tell them?

### Before you answer: what a tree cannot tell you

| Limitation | Why it bites here |
| --- | --- |
| **Clustering is not transmission** | Identical genomes are consistent with direct transmission, with a shared unsampled source, and with two independent infections from the same circulating lineage. |
| **Direction is not recoverable** | A tree does not say who infected whom. Collection date reflects when someone was tested, not when they were infected. |
| **Sampling determines what you see** | You have twelve genomes. If half the cases were never sequenced, an apparent chain may be two chains joined by a missing intermediate. |
| **Resolution has a floor** | SARS-CoV-2 accumulates roughly a couple of mutations a month. Over a few weeks many transmission pairs are genetically indistinguishable. |
| **Sequencing error mimics signal** | A single-base difference may be a miscall at a low-coverage position. Check the depth before building an argument on one SNP. |

### What a defensible conclusion sounds like

Not "these three cases form a transmission chain", but something closer to: two
groups of three genomes each are separated from the rest by several
substitutions and supported at 99-100% bootstrap, which is consistent with
separate introductions; a third group is visible but supported in the low
seventies, below the threshold at which the number means anything; within
groups, genomes differ by 0 to 1 substitutions, which is consistent with
transmission in the facility without establishing it; and one genome appears
identical to two others only because a third of its positions are N, so that
identity is manufactured by missing data rather than observed.

That is less satisfying than the answer management asked for. It is the one the
data supports, and saying so clearly — including which parts are solid and which
are artefacts of coverage — is the actual skill this session teaches.

## 7. Closing discussion

These are the questions the session exists for. Each has a defensible answer and
a tempting wrong one.

**1. Your consensus genome is 94% complete. Is it fit for purpose?**

**2. Two patients have identical viral genomes. Did one infect the other?**

**3. A variant appears at 8% frequency in your sample. Do you report it?**

**4. Your phylogeny has a node with 45% bootstrap support. How do you present
it?**

**5. Infection control asks you for a yes or no: was this one outbreak?**

If you take one sentence away from today, take this one: phylogenetic clustering
does not establish transmission. Sampling, epidemiological context, timing,
recombination and sequencing quality all bear on interpretation, and the tree is
one input among them.

## Appendix: everything in one script

`scripts/one_sample.sh` runs the whole chain for a single sample — host removal,
QC, trim, map, depth, consensus, variants — with timings and a comparison
against the published metrics:

```
cd $WORKDIR
sbatch scripts/one_sample.sh            # A07 by default, edit SAMPLE to change
```

All twelve as a job array:

```
sbatch --array=1-12 scripts/one_sample.sh
```

The step-by-step scripts are the better way to learn what each stage does, which
is why the walkthrough uses them. This one is for when you want the result, for
re-running on your own data later, and as a smoke test after anything changes on
the cluster.

## Troubleshooting

| Symptom | Cause | Fix |
| --- | --- | --- |
| `command not found` | Modules not loaded, or lost after reconnecting | Re-run the `module load` line; `module list` to confirm |
| `fastqc: error: unrecognized arguments: -o` | FastQC 0.12+ dropped the short flag | Use `--outdir` |
| `salloc` hangs | Queue wait | Tell a facilitator; use the reservation or work from `precomputed/` |
| Job killed, no clear error | Out of memory or past the time limit | `sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed` |
| Permission denied writing output | Still in the read-only shared folder | `cd $WORKDIR` |
| `[E::faidx] Failed to open` | `$REF` unset in a new shell | Re-export it |
| Mapping rate suspiciously high | You are mapping host-depleted reads | Expected — see 4b |
| `ivar` column numbers wrong | Column order differs by version | `head -1 results/variants/A01.tsv` and adjust |
| Anything still running after 3 minutes | The cluster, not you | Use `precomputed/` |

## Tools

All from `/software/modulefiles` unless noted.

| Need | Module |
| --- | --- |
| Read QC | `fastqc/0.12.1`, `multiqc/1.35` |
| Trimming | `fastp/1.3.7` |
| Host removal | `read_keep/0.3.0` (binary: `readItAndKeep`) |
| Mapping | `minimap2/2.31`, `samtools/1.22.1` |
| Consensus, variants | `ivar/1.4.4` |
| Alignment | `mafft/7.525` |
| Tree | `iqtree2/2.4.0` |
| FASTA handling | `seqkit/2.14.0` |

There is no csvtk — table formatting uses `cut` and `column -t`. There is no
sra-tools, so reads cannot be fetched on the cluster; the facilitator stages the
bundle in advance. Clade and lineage assignment is not part of this session: the
structure you need comes from the variants and the tree.

## Taking this further

- The ARTIC field bioinformatics pipeline runs the amplicon equivalent of step 4
  end to end, and is what a surveillance laboratory would deploy.
- Nextstrain and UShER place genomes in a global rather than local context, which
  answers a different question from the one asked here.
- For within-host diversity and drug resistance — HIV, influenza, hepatitis B —
  the consensus-versus-population distinction in step 4c becomes the central
  analytical problem rather than a caveat.
