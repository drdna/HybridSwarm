# Hybrid-swarm simulations for testing admixture tools

These scripts simulate hybrid swarms of a **haploid, heterothallic fungus**, write the
input files for ADMIXTURE, STRUCTURE and sNMF, and compare each tool's ancestry estimates
with the true ancestry of every simulated strain.

The question they are built for: *is it harder to recover true admixture proportions in
young swarms (few rounds of intermating) than in older ones?*

---

## Contents

| File | What it does | Runs on |
|---|---|---|
| `haploswarm.py` | Simulates swarms from a model file (sources, backcrosses, outcrosses, survival, clonality, several populations) or a one-line shortcut, with several chromosomes, 4-spore tetrads and mating types; writes all tool inputs and the truth | Mac or cluster |
| `models/` | Example model files | – |
| `genome_example.tsv` | Chromosome lengths, map lengths and SNP counts for `haploswarm.py`. **Edit for your species.** | – |
| `split_generations.py` | Splits a simulated data set into one input set per generation (each with the same reference isolates) | Mac or cluster |
| `run_admixture_by_gen.sh` | SLURM array: PLINK2 + unsupervised and supervised ADMIXTURE for every generation | Cluster |
| `control_admixture.py` | Positive control: ADMIXTURE-conforming data from the same parents, matched to each generation's swarm | Mac |
| `reciprocal_crossovers.py` | Counts reciprocal crossovers (meiotic-sib signatures) per generation, in simulated ancestry or a real VCF | Mac or cluster |
| `plot_admixture.R` | Plots one set of ADMIXTURE runs against the truth | Mac |
| `compare_generations.R` | Runs `plot_admixture.R` for every generation and summarises accuracy across generations and scenarios | Mac |
| `make_admix_inputs.py`, `recseg_to_tracts.py`, `check_tetrads.py`, `haploid_cross/` | Earlier route using AdmixSim2 (one chromosome only). See [Older AdmixSim2 route](#older-admixsim2-route). | – |

**Recommended design.** Analyse the swarm sample on its own, **unsupervised**, one generation
at a time, as in real studies where the source lineages are unknown and are inferred as the K
clusters. The simulation's truth is used only afterwards, to score the estimates. If you want
references at all (for example to test supervised ADMIXTURE), they must be the **two parents
themselves** (`--refs 0 --include-parents`), which have zero SNP differences from the swarm's
sources. Do not use `--refs N` (other isolates drawn from each lineage's allele frequencies):
they are not the swarm's sources, and they add SNPs that vary only among the references, which
made the swarm look like a separate population in earlier runs.

**Requirements**
- Python 3 with `numpy`
- R with `ggplot2`, `dplyr`, `tidyr`
- On the cluster: `plink2` and `admixture`

---

## Quick start

```bash
# 1. Simulate (Mac): a swarm founded by one tetrad, sampled every generation
#    (the two parents are written too, as the only possible references)
python haploswarm.py --genome genome_example.tsv --founder-tetrads 1 \
    --sizes 200,1000,5000x9 --output-gens all --sample 200 --refs 0 --include-parents --seed 11 \
    --out bench3/tetrad

# 2. One swarm-only input set per generation (Mac)
python split_generations.py bench3/tetrad --no-refs --outdir bench3/tetrad_swarmonly

# 3. Copy to the cluster and run ADMIXTURE on every generation (12-task array)
rsync -av bench3/tetrad_swarmonly run_admixture_by_gen.sh username@cluster.address:path/to/
sbatch run_admixture_by_gen.sh tetrad_swarmonly 5 10 allmaf

# 4. Copy the results back (Mac)
rsync -av username@cluster.address:path/to/tetrad_swarmonly/ bench3/tetrad_swarmonly/

# 5. Plot and summarise (Mac)
Rscript compare_generations.R --set tetrad=bench3/tetrad_swarmonly --set perithecium=bench3/perithecium_swarmonly
```

---

## 1. Simulating swarms: `haploswarm.py`

### The model

A run is described by a **model file** (`--model`), written in the spirit of AdmixSim2's
`.mod` files. For the simplest case, a single swarm founded by parents A × B, the shortcut
options (`--sizes`, `--founder-tetrads`, `--clonal`, `--survive`, ...) build the equivalent
model for you.

```
# comments start with #
source A  MAT=1  lineage=A                # fixed parental genomes; can mate in any generation
source B  MAT=2  lineage=B
source C  MAT=2  lineage=C   private=50   # optional: 50 SNPs/Mb that only this source carries

*0      swarm  size=4     mates=A,B  ferts=1  meioses=1  mdist=fixed   # one meiosis founds the swarm
*1      swarm  size=50    mates=swarm                                  # intermating among the 4 sibs
*2-6    swarm  size=500   mates=swarm:0.9,A:0.1  survive=0.5  clonal=0.3
*7-11   swarm  size=1000  mates=swarm:0.95,C:0.05 survive=0.5
```

**`source` lines** define fixed parental genomes. Each has a mating type and a lineage (the
ancestry label used in the truth tables; default: its name). Two sources of the same lineage
are genetically identical unless given `private=` SNPs. Sources never change and can mate in
any generation, so they model backcrosses (mating with an original parent) and outcrosses
(mating with a new parent).

**Generation lines** (`*G` or `*G1-G2`, then a population name) define a population in
those generations. Generations must be consecutive. Several populations can exist at once
(see `models/two_swarms_merge.mod`). Each generation of a population with `size=N` is made of:

| Part | How many | Where from |
|---|---|---|
| Survivors | `survive` × N | The same individuals as in this population's previous generation, chosen at random without replacement. Fungi don't die after mating, so a strain can persist, and keep mating, for many generations. |
| Clones | `clonal` × N | Copies of random strains of this population's previous generation (with replacement) |
| New spores | the rest | Fertilizations, then meiosis (below) |

**Fertilizations.** Each pairs one MAT1-1 and one MAT1-2 partner. Each partner is drawn from
the populations and sources listed in `mates=`, in proportion to their weights (among those
that have that mating type), then at random within that population. A strain can mate many
times. For example, `mates=swarm:0.9,A:0.1` makes about 10% of the MAT1-1 partners the
original parent A, a backcross. `mates=swarm:0.95,C:0.05` makes about 5% of the MAT1-2
partners the new parent C, an outcross. Each fertilization gives Poisson(`meioses`)
tetrads (at least 1), or exactly `meioses` with `mdist=fixed`. Fertilizations continue
until there are enough spores (`ferts=auto`, the default) or exactly `ferts=F` times; the
needed spores are then drawn at random and the rest discarded.

| Key | Meaning | Default |
|---|---|---|
| `size` | Strains in the population at this generation | required |
| `mates` | `pop[:weight],...` populations (at the previous generation) or sources that supply partners | the population itself |
| `survive` | Share that are surviving strains of the previous generation | 0 |
| `persist` | Alternative to `survive`: share of the previous generation's strains that survive into this one (e.g. `persist=0.5` = half of them live on) | 0 |
| `clonal` | Share that are clonal copies of the previous generation | 0 |
| `ferts` | Fertilizations, or `auto` | `auto` |
| `meioses` | Tetrads per fertilization (Poisson mean) | 20 |
| `mdist` | `poisson`, or `fixed` for exactly `meioses` | `poisson` |

If a generation's mates have lost one mating type, its sexual part is produced clonally
instead (with a warning).

Example models are in `models/`:
- `one_meiosis_backcross_outcross.mod`: the model above.
- `two_swarms_merge.mod`: swarms A × B and A × C founded separately, kept apart for four
  generations, then merging.
- `single_tetrad_swarm.mod`: one A × B fertilization; the 4 sibs of its tetrad intermate;
  every later fertilization gives exactly 10 meioses, and half of each generation survives
  into the next (`persist=0.5`).

**Mating type.** Unlinked by default: it then segregates as if it sat right next to a
centromere. Use `--mat CHROM:POS` to place it on a chromosome.

**SNPs.** Genomes are stored as ancestry tracts labelled by source, so memory stays small.
SNPs are added afterwards, so the SNP model never changes the strains or their ancestry. With
`--parent-snps-per-mb`, exactly that many differences per Mb between the lineages are placed
at random. Otherwise a lineage allele-frequency model is used (Balding–Nichols, `--fst`, over
`n_snps` sites per chromosome). At each SNP, a strain carries the allele of the source its
tract came from.

### How crossovers are assembled

For each chromosome in each meiosis:

1. The number of crossovers in the bivalent is Poisson(2L), where L is the map length in
   Morgans. Their positions are uniform along the chromosome (no crossover interference).
2. Each crossover joins one chromatid from each homolog, chosen at random (no chromatid
   interference). In the code, the four product chromatids are traced from one end of the
   chromosome to the other. At each crossover, one product currently on a homolog-1 strand
   and one on a homolog-2 strand exchange everything beyond that point. Tracing from the
   telomere or outward from the centromere gives the same four recombinant chromatids; only
   the labels differ.
3. The centromere decides segregation. The two products whose centromere comes from the
   same homolog go together at the first division, to spores 1+2 or 3+4. The second division
   then sends one to each spore. Chromosomes assort independently.

Checks on 20,000 meioses with exactly two crossovers gave two-strand, three-strand and
four-strand doubles in proportions 0.247 : 0.498 : 0.255 (expected 1 : 2 : 1). Spores 1+2
always shared one centromere and 3+4 the other. Each crossover therefore leaves a reciprocal
pair of recombinant spores in its tetrad.

### Genome file

Tab- or space-separated, one row per chromosome. `#` lines are ignored.

```
chrom  length_bp  map_cM  n_snps  [centromere_bp]
1      8000000    200     5000
```

`n_snps` is the number of candidate SNP sites. Sites that aren't polymorphic among the
output samples are dropped. The centromere defaults to the middle of the chromosome.

### Options

| Option | Meaning | Default |
|---|---|---|
| `--genome` | Genome file (required) | – |
| `--model` | Model file (sources and generations). Replaces the shortcut options below. | – |
| `--sizes` | Shortcut: strains per generation after founding, e.g. `200,1000,5000x9` (`NxR` = N repeated R times) | – |
| `--founder-tetrads` | Meioses of the founding fertilization | 1 |
| `--founders` | Spores kept as founders | all |
| `--meioses` | Mean meioses per fertilization in later generations | 20 |
| `--ferts` | Fertilizations per generation | just enough spores for the next generation |
| `--clonal` | Shortcut: share of each generation that is clonal copies, 0–1. One value, or one per generation, e.g. `0,0,1x5,0x4` | 0 |
| `--survive` | Shortcut: share of each generation that is surviving strains of the previous one, 0–1 | 0 |
| `--mat CHROM:POS` | Position of the mating-type locus | unlinked |
| `--parent-snps-per-mb` | Place exactly this many SNPs per Mb between parent A and parent B, at random positions. Overrides `n_snps` and `--fst`; any references become copies of the parents. | off |
| `--fst` | Differentiation between lineages A and B (used only without `--parent-snps-per-mb`) | 0.6 |
| `--refs` | Other isolates per lineage (not the sources; not recommended) | 0 |
| `--include-parents` | Also output every source as a reference sample (`parent<NAME>`) | off |
| `--output-gens` | Generations to sample, e.g. `0,1,3,11`, or `all` (0 = founders) (required) | – |
| `--sample` | Strains sampled per output generation and population (one value, or one per generation) | 200 |
| `--thin` | Keep at most one SNP per this many bp | all |
| `--linkage` | Add map distances to the STRUCTURE file and set `LINKAGE 1` | off |
| `--seed` | Random seed | – |
| `--out` | Output prefix (required) | – |

Examples:

```bash
# founded by one tetrad (4 spores)
python haploswarm.py --genome genome_example.tsv --founder-tetrads 1  --sizes 200,1000,5000x9 --output-gens all --sample 200 --include-parents --seed 11 --out bench3/tetrad
# founded by one perithecium (20 meioses, 80 spores)
python haploswarm.py --genome genome_example.tsv --founder-tetrads 20 --sizes 200,1000,5000x9 --output-gens all --sample 200 --include-parents --seed 11 --out bench3/perithecium
# a model file: one meiosis, backcrossing to A, then an outcross to C
python haploswarm.py --genome genome_example.tsv --model models/one_meiosis_backcross_outcross.mod --parent-snps-per-mb 300 --include-parents --output-gens all --sample 50 --seed 5 --out bench4/bc_oc
```

Each takes about 2 minutes and under 300 MB of memory.

A swarm founded by **one meiosis**, then grown mainly clonally with a little intermating:

```bash
python haploswarm.py --genome genome_example.tsv --founder-tetrads 1 \
    --sizes 20,100,500,1000x5 --clonal 0,0.9x7 --meioses 5 \
    --output-gens all --sample 50 --refs 0 --include-parents --seed 3 --out bench4/one_meiosis
```

Generation 0 is the 4 spores of the single meiosis. With the default `--ferts`, each
generation's sexual part comes from just enough fertilizations (often only one or two),
so Ne stays very low. Set `--ferts` and `--meioses` to control this.

**Setting the divergence between the parents.** By default, the parents' SNPs come from the
two-lineage model: 35,000 candidate sites, of which about 36% differ between the parents
(about 315 SNPs per Mb in the example genome). To set the density directly, use
`--parent-snps-per-mb`, e.g. `--parent-snps-per-mb 2000` gives 81,000 SNPs in the 40.5 Mb
example genome. The SNP model is applied after the swarm is simulated, so the same `--seed`
gives the same strains and ancestry at any density, which lets you test the effect of marker
density alone. Run time and file sizes grow with the number of SNPs (2,000 per Mb took about
5 minutes for 2,200 strains).

### Output (`<out>.*`)

| File | Contents |
|---|---|
| `dip.vcf.gz` | Genotypes coded as homozygous diploid (`0/0`, `1/1`), for PLINK → ADMIXTURE |
| `vcf.gz` | Haploid genotypes |
| `geno` | sNMF / LEA input: one row per SNP, one character per sample (0/1) |
| `str`, `mainparams`, `extraparams` | STRUCTURE input (`PLOIDY 1`). References have `POPFLAG 1`, so setting `USEPOPINFO 1` gives a supervised run. |
| `pop` | ADMIXTURE `--supervised` labels: `A`, `B`, or `-` |
| `samples.tsv` | Sample order. This is the same in every file, and in the PLINK `.fam`. |
| `truth.tsv` | True ancestry per strain: `A_frac_bp` (share of the genome) and `A_frac_snp` (share of SNPs) for the first lineage, plus `MAT` and `n_tracts`. With more than two lineages, also `<lineage>_frac_bp` for each. |
| `anc.tsv.gz` | True ancestry at every SNP (1 = first lineage, 0 = other), so the truth can be recomputed for any SNP subset |
| `anc_lineage.tsv.gz` | With more than two lineages: the lineage at every SNP |
| `tracts.tsv` | Ancestry tracts of every sampled strain, with the lineage and source of each |
| `strains.tsv` | For every sampled strain: generation born, mating type, fertilization (family), tetrad, parents, permanent ID (`uid`) and origin (`spore` or `clone`). A survivor has a later sample name than its `uid`. |
| `generations.tsv` | Per generation and population: size, survivors, clones, fertilizations, meioses, mating types, mean share of each lineage, tracts per strain |
| `log` | The command used |

Sample names are `g<gen>_s<n>` for a population called `swarm`, `<pop>_g<gen>_s<n>` for other populations, and `parent<NAME>` for sources.

---

## 2. One input set per generation: `split_generations.py`

```bash
python split_generations.py bench2/tetrad                     # all generations -> bench2/tetrad_bygen/
python split_generations.py bench2/tetrad --gens 1,3,6,11      # some generations
python split_generations.py bench2/tetrad --max-strains 40     # subsample each generation to 40 strains
python split_generations.py bench2/tetrad --no-refs            # swarm strains only
python split_generations.py bench4/merge --pop merged          # one population of a multi-population model
```

Each `g<gg>/` folder holds `<name>_g<gg>.*`: `dip.vcf.gz`, `vcf.gz`, `geno`, `str`
(+ `mainparams`, `extraparams`), `pop`, `samples.tsv`, `truth.tsv`, `anc.tsv.gz`, and
`refs.txt` (reference IDs, for `plink2 --keep`). SNPs that are monomorphic within a
generation are dropped.

The split `.str` files keep every SNP (about 26,000), which is too many for STRUCTURE.
Thin them first, or use the SNPs ADMIXTURE kept.

---

## 3. ADMIXTURE on the cluster: `run_admixture_by_gen.sh`

```bash
sbatch run_admixture_by_gen.sh <folder> [KMAX=5] [RUNS=10] [PRUNE=allmaf]
```

- **SLURM array:** each task handles one `g*` folder. Set `#SBATCH --array=0-(number of generations - 1)`.
  Without SLURM, it loops over all folders.
- **Software:** if your cluster needs it, edit the `module load` line, or wrap the commands
  in your `singularity run` calls.
- **`PRUNE`** sets how SNPs are chosen:
  - `allmaf` (default): MAF ≥ 0.05 across all samples of the generation (swarm + references), no LD pruning.
  - `all`: the same, plus `--indep-pairwise 50 10 0.1` on all samples (a standard pipeline).
  - `refsmaf` / `refs`: MAF (and for `refs`, LD pruning) computed on the references only.
  - `none`: every SNP polymorphic in the generation.

  **Which SNPs are kept changes the answer.** Selecting on the references only keeps SNPs that
  vary among the references but are fixed in the swarm, because both parents carried the same
  allele. In generation 1 of the tetrad set, 60% of the 3,405 `refs` SNPs were fixed in the swarm.
  ADMIXTURE then put every swarm strain in its own cluster (Q ≈ 0 for A in every generation),
  while the pooled analysis, where the swarm dominated the MAF filter and none of the SNPs were
  fixed in the swarm, gave intermediate estimates. Part of this is real (a swarm founded by two
  genomes lacks its lineages' diversity), but the simulated references are far more diverse
  than real clonal lineages, which exaggerates it.

Output in each `g<gg>/`:

```
<tag>_ADMIX.bed/.bim/.fam                  the SNPs ADMIXTURE used
unsup/<tag>_ADMIX_K<K>_run<r>.Q/.P/.log    unsupervised, K = 1..KMAX, RUNS seeds each
sup/<tag>_ADMIX_K2_run<r>.Q/.P/.log        supervised with the A / B references
```

ADMIXTURE writes the `.Q` and `.P` files only when a run finishes. A run that exits with an
error is logged (`K=… run … FAILED` in the SLURM output, the error in its `.log`) and skipped,
so one bad seed doesn't end the generation's task.

**Slow runs at high K.** When K is larger than the real structure (here, more than 2 sources),
the extra clusters have nothing real to fit. The likelihood surface is flat, and some seeds
crawl for hours before the gain per iteration drops below ADMIXTURE's 1e-4 stopping threshold.
This is worse with many identical samples (clonal sources) and in generations with few
distinct genotypes: generation 0 of a single-tetrad swarm has only 6. A missing panel at high K
usually means that run hadn't finished. To save time, set KMAX to 3 or 4, or split K across
array tasks. Check that the best few seeds per K agree in `_runs.tsv`.

### Running ADMIXTURE your own way

Any naming works for the plotting scripts, as long as:
- each run's `.Q` file starts with a common prefix (e.g. `tetrad.dip_ADMIX_K2_run1.Q`)
- its screen output is saved next to it with the same name ending `.log` or `.out`. The
  CV error and log-likelihood only appear on screen.
- you keep the `.fam`/`.bim` that ADMIXTURE actually ran on (after pruning)

If you add extra samples (e.g. a reference genome such as CMW25467), the `.pop` file needs
a `-` line in that sample's position. The plotting scripts show such samples as "other".

---

## 4. Plotting one set of runs: `plot_admixture.R`

```bash
Rscript plot_admixture.R --prefix tetrad.dip_ADMIX_K --fam tetrad.dip_ADMIX.fam \
    --samples tetrad.samples.tsv --truth tetrad.truth.tsv --anc tetrad.anc.tsv.gz \
    --label "single-tetrad founders"
```

| Option | Meaning |
|---|---|
| `--prefix` | Start of the `.Q` file names (required). K is read from the number of columns. |
| `--dir` | Folder with the `.Q` and log files (default `.`) |
| `--fam` | `.fam` of the ADMIXTURE input, **after pruning** (sets the sample order) |
| `--bim` | `.bim` of the ADMIXTURE input (default: the `.fam` name with `.bim`) |
| `--samples`, `--truth` | From `haploswarm.py` or `split_generations.py` (required) |
| `--anc` | Per-SNP truth. Recommended: it lets the SNP-based truth use only the SNPs ADMIXTURE saw. |
| `--basis` | `bp` (default): truth = share of the **genome** inherited from A. `snp`: share of the SNPs used. |
| `--out` | Output prefix (default `<dir>/<prefix>_plots`) |
| `--label` | Text added to plot titles |
| `--kmin`, `--kmax` | Only use runs with K in this range (default: every K found) |

With several runs per K, the run with the highest log-likelihood is plotted. The script
compares the `.bim` with the `.P` files (one row per SNP used). If you passed the
pre-pruning `.fam`, it switches to the matching `.bim` and tells you.

**Output**

`<out>.pdf`, one plot per page:
1. Cross-validation error for every run
2. Log-likelihood for every run
3. Admixture bar plots for every K, with the true ancestry on top
4. Cluster membership (Q ≥ 0.9 shown as "pure")
5. Estimate vs truth at K = 2, per generation
6. Accuracy by generation: bias, RMSE, r, and the true spread among strains
7. Mean ancestry by generation: true vs estimated
8. Where the error comes from: the markers vs the tool
9. Error vs truth

`<out>_comparison.tsv` has, per strain: `A_frac_bp`, `A_frac_snp`, `truth`, `est_A` (K = 2),
the spread across runs, and the errors. `<out>_runs.tsv` has the log-likelihood and CV error
of every run.

**Cluster colours.** ADMIXTURE numbers clusters arbitrarily in every run, so the script
reorders them: C1 is the cluster where the A references score highest (or, without
references, the one best correlated with the truth), and C2 is the cluster where the B
references score highest. Any further clusters follow, by overall size. With references,
blue = A and red = B in every generation. Clusters C3 and above can't be matched between
separate analyses, so their colours aren't comparable between generations.

**Reading the K = 2 results.** If the B references *also*
score mostly in C1, the two clusters are not A and B, and the script prints a **WARNING**
and adds it to the plot titles. In that case `est_A` is not an estimate of A ancestry.

---

## 5. Comparing generations and scenarios: `compare_generations.R`

```bash
Rscript compare_generations.R --set tetrad=bench2/tetrad_bygen --set perithecium=bench2/perithecium_bygen
```

For every `g*/unsup` and `g*/sup` folder with `.Q` files, this runs `plot_admixture.R`. It
skips folders already plotted unless you add `--redo yes`. Then it writes
`generations_summary.pdf` and `generations_summary.tsv` (change the name with `--out`).
`--kmin` and `--kmax` (e.g. `--kmin 2 --kmax 4`) limit every plot to that range of K.
After updating `plot_admixture.R`, add `--redo yes`: the per-generation tables are reused
otherwise, and keep the old cluster order.

The summary plots, each by generation and scenario, for supervised vs unsupervised:
- bias, RMSE, correlation with the truth, and the share of strains called pure
- RMSE against the true spread among strains, and against the error from the markers alone
- true vs estimated mean ancestry
- whether the K = 2 clusters separate the A and B references
- the K chosen by cross-validation

`plot_admixture.R` must be in the same folder as `compare_generations.R`.

---

## Positive control: `control_admixture.py`

Does ADMIXTURE recover ancestry when the data follow its own model? From a `haploswarm.py`
run made with `--include-parents`, this writes, for every generation:

| Set | Strains | Sources |
|---|---|---|
| `swarm_sources/gNN` | the real swarm strains of that generation | `--sources` clonal copies each of parent A and parent B (identical to the parents) |
| `model/gNN` | one strain per swarm strain, with the **same true A fraction q**, but every SNP drawn independently (A with probability q): no linkage, no relatedness | the same clonal copies |

Both sets use every SNP where the parents differ. The only difference between them is the
swarm's linkage and relatedness, which ADMIXTURE's model leaves out.

```bash
python control_admixture.py bench5/single_tetrad --sources 20 --outdir bench5/control
sbatch run_admixture_by_gen.sh control/swarm_sources 5 10 none     # on the cluster
sbatch run_admixture_by_gen.sh control/model 5 10 none
Rscript compare_generations.R --set swarm=bench5/control/swarm_sources --set model=bench5/control/model
```

The `.pop` files label the clonal sources `A` and `B`, so the cluster script also runs
supervised ADMIXTURE. Checks on generation 5: in the model set, the ancestry correlation
between neighbouring SNPs is 0.000 and between strains 0.000, compared with 0.996 and 0.23
in the swarm. Realised A fractions match their targets to within 0.003.

**Results (bench5, single-tetrad swarm, all SNPs, unsupervised, K = 2–3).**
- **Model set:** K = 2 matches the truth closely in every generation. At K = 3 the third
  cluster gets only noise-level slivers; A and B are recovered as at K = 2.
- **Swarm set:** K = 2 is noisy in generations 1–3 and biased later. At K = 3, from
  generation 2 onward, the swarm strains are modelled as A (blue) plus a **third cluster
  (green) that the B parents don't share**. Almost none of the B parents' cluster (red)
  appears in the swarm. ADMIXTURE infers an unsampled third ancestor that doesn't exist.
- **Why:** founder drift and linkage. From generation 3 onward, 1,898 SNPs carry only the
  A allele in every swarm strain, and 208 only the B allele. At those sites parent B carries
  alleles that no swarm strain has, so the swarm's own allele-frequency profile fits better
  than parent B. Shared tracts and relatedness count each haplotype block many times over,
  as if its SNPs were independent. The green fraction exceeds the true B fraction, so the
  ghost also absorbs A ancestry. It describes the swarm itself, not drifted B. This is the
  bottleneck pitfall described by Lawson, van Dorp & Falush (2018, *Nat. Commun.*), here
  caused by a single founding meiosis.
- **Still to check:** whether CV error favours K = 3 for the swarm but not the model set,
  and whether the ghost persists when tracts are broken (shuffled SNP order, or strict LD
  pruning) but fixed blocks are kept.

## Reciprocal crossovers: `reciprocal_crossovers.py`

The two reciprocal products of a crossover switch ancestry at the same place in opposite
directions (A→B in one, B→A in the other). Finding such pairs in a sample indicates meiotic
siblings, or descendants of them.

```bash
# simulated: true ancestry per SNP, one row per generation
python reciprocal_crossovers.py --anc bench3/tetrad.anc.tsv.gz --samples bench3/tetrad.samples.tsv --out bench3/tetrad_recip
# real (or simulated) isolates: ancestry from matching two parental / proxy samples
python reciprocal_crossovers.py --vcf isolates.vcf.gz --parent-a P1 --parent-b P2 --out recip
```

Breakpoints are located to the interval between adjacent informative SNPs. Strains with
identical ancestry (clones) are counted once unless you add `--keep-clones`. The summary
reports, per generation: strains, distinct genotypes, breakpoints, **reciprocal events**
(SNP intervals where one strain switches A→B and another B→A), **expected by chance** (the
distinct A→B and B→A intervals placed independently), and intervals shared in the same
direction (shared descent). `<out>.events.tsv` lists every reciprocal event with the strains
involved.

Examples (about 315 SNPs per Mb):

| Scenario | Generation | Reciprocal events | Expected by chance |
|---|---|---|---|
| tetrad founders, all sexual, 200 → 5,000 strains | 0 (the 4 founders) | 15 (all 15 crossovers) | 0 |
| | 1 | 700 | 55 |
| | 11 | 257 | 124 (5 of the founding meiosis's 15 crossovers still detected) |
| one meiosis, then 90% clonal, 20 → 1,000 strains | 1 | 51 | 0.6 |
| | 8 | 57 | 2.5 |

The chance expectation falls with SNP density, so denser SNPs (`--parent-snps-per-mb`)
separate true reciprocal pairs from coincidences better.

## Notes and caveats

- **Truth.** The default truth is the share of each strain's genome inherited from A (by bp).
  The SNP-based truth shows how much error the marker set alone causes. With heavy LD pruning
  (e.g. ~145 SNPs) this marker error is about 0.03 per strain, which is a large part of the
  true spread among strains in older swarms.
- **Analyse each generation on its own, swarm only.** Pooling all generations (about 2,200
  strains) mixes populations that descend from one another, and adding reference isolates
  that are not the parents changes which SNPs are kept and what the clusters mean.
- **Matching clusters to A without references.** ADMIXTURE's cluster numbers are arbitrary.
  With no references, the plotting scripts call "A" the cluster that correlates best with the
  true A fraction. This is the most favourable labelling possible; a real study has no truth
  to do it with.
- **Family structure.** With the default `--ferts`, a generation of 5,000 strains comes from
  only about 60 fertilizations. Increase `--ferts` or lower `--meioses` to change this.
- **Recombination** is uniform along each chromosome, with no interference. Selfing and
  clonal reproduction are not modelled.
- **Replication.** Each scenario has been run with one seed so far. Run several seeds
  (`--seed`) before drawing conclusions.

---

## Older AdmixSim2 route

Before `haploswarm.py`, the haploid cross was modelled with AdmixSim2, which handles one
chromosome only. The patched copy `haploid_cross/AdmixSim2_tetrad.py` adds `--tetrad`
(4-spore tetrads, pairing one MAT1-1 and one MAT1-2 spore per zygote) and `--mat POS`.

- `make_admix_inputs.py <AdmixSim2 output prefix> --ref cross.hap:A:1 --ref cross.hap:B:2 --out x`
  writes the same tool inputs and truth tables from an AdmixSim2 run.
- `recseg_to_tracts.py` turns a `.recseg` file into merged A/B tracts.
- `check_tetrads.py` checks tetrad segregation.

For new work, use `haploswarm.py`.
