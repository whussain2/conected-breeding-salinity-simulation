# Connected Breeding vs closed-pool genomic recurrent selection — simulation code

Companion code for:

> Hussain, W., Anumalla, M., Catolos, M., Ramos, J., Sta. Cruz, M. T., Zhang, X.,
> Sreenivasulu, N., and Bhosale, S. **An Integrated Framework to Develop and Deliver
> Salt-Tolerant Rice Varieties for Coastal Ecologies.** *Plant Communications*.

This repository contains everything needed to reproduce **Figure 4C, D**,
**Supplemental Figure S2**, **Supplemental Table S5** (simulated columns) and
**Supplemental Table S6**, together with the raw per-replicate output from which every
reported value is derived.

---

## What the simulation does and does not model

**Read this before interpreting any output.**

Grain yield is modelled as a single, purely additive, highly polygenic trait
(~1,000 QTL). **Salinity tolerance is not modelled.** Specifically there is:

* no salinity-tolerance locus or tolerance sub-trait;
* no genotype-by-salinity interaction;
* no stage-specific tolerance (no seedling / vegetative / reproductive distinction);
* no differential survival, mortality or phenotype truncation of susceptible
  genotypes — which is precisely the phenomenon the Transition from Trait to
  Environment (TTE) design is intended to eliminate, and therefore the reason this
  simulation **cannot test TTE**;
* no environmental covariate in the baseline scenario (genotype-by-environment
  interaction is added only in the sensitivity scenarios).

What the simulation *does* test is the population-genetic mechanism on which
Connected Breeding depends: whether adding a managed, mean-neutral inflow of donor
diversity to recurrent genomic selection sustains genetic gain and additive genetic
variance better than recurrent genomic selection confined to a closed elite pool.

The baseline comparison is **not resource-neutral** and favours Connected Breeding by
construction. `R/02_scenarios.R` defines the equal-budget, development-lag and
cycle-time-penalty scenarios that remove that advantage; see Table S6.

---

## Requirements

* R ≥ 4.3
* [AlphaSimR](https://cran.r-project.org/package=AlphaSimR) ≥ 2.1.0

```r
install.packages("AlphaSimR")
```

No other packages are needed; plotting uses base graphics.

---

## Reproducing the results

```bash
# 1. main grid: baseline + all sensitivity + resource-matched scenarios
#    args: <nRep> <firstSeed> <outfile>
Rscript R/03_run.R 10 9000 out/raw_A.csv
Rscript R/03_run.R 10 9010 out/raw_B.csv     # a second block of replicates

# 2. the additional established strategies used in Table S5
Rscript R/04_run_strategies.R 10 9000 out/raw_strategies.csv

# 3. analysis, Figure S2 and Table S6
Rscript R/05_analyse.R

# 4. per-strategy summary for Table S5
Rscript R/06_summarise_strategies.R
```

The two `03_run.R` calls are independent and may be run in parallel; they use
disjoint seed blocks and are pooled by `05_analyse.R`.

### Reproducibility

Every run is seeded per replicate. **AlphaSimR parallelises recombination using
per-thread random-number streams, so bit-for-bit reproducibility requires
single-threaded execution.** The scripts set `SP$nThreads <- 1L` and
`OMP_NUM_THREADS=1`; do not override this if you need to reproduce published values
exactly.

---

## Files

| File | Purpose |
|---|---|
| `R/01_functions.R` | Founder simulation, trait scaling, donor sampling, and the breeding-programme engines (`GSRS`, `CB`, `PSRS`, `MAGIC`) |
| `R/02_scenarios.R` | Baseline, sensitivity and resource-matched scenario definitions |
| `R/03_run.R` | Driver for the scenario grid; writes one row per cycle per scenario per replicate |
| `R/04_run_strategies.R` | Phenotypic recurrent selection and the one-time multi-parent resource |
| `R/05_analyse.R` | Paired contrasts, Table S6, Figure S2, and the numbers quoted in the main text |
| `R/06_summarise_strategies.R` | Per-strategy summary for Table S5 |
| `out/` | Raw per-replicate output, summary tables, and figures |

## Key parameters

All defaults are in `defaultParams()` in `R/01_functions.R` and correspond to the
parameter table in Supplemental File S1:

| Parameter | Value |
|---|---|
| Genome | 12 chromosomes × 1.2 Morgans (~1,440 cM) |
| QTL (grain yield) | ~1,000 additive QTL (83 per chromosome) |
| Marker panel | 1,200 SNPs (100 per chromosome) |
| Elite base mean | 6.0 t ha⁻¹ |
| Base additive variance | 0.5 (t ha⁻¹)², i.e. 1 base σ_A = 0.707 t ha⁻¹ |
| Elite − donor differential | 2.5 t ha⁻¹ (set by construction; see `makeFounders`) |
| Heritability (training phenotypes) | 0.30 |
| Elite / donor founders | 60 / 60 |
| Crosses × progeny per cycle | 30 × 40 (1,200 candidates) |
| Parents recycled per cycle | 20 |
| GS training set | 200 lines per cycle, sliding window of 1,000 records |
| Elite-equivalent bridge lines per cycle | 2 |
| Cycle length / horizon | 3 years × 20 cycles = 60 years |

### A note on the elite base variance

The trait is created with variance 1 in the *unselected* founder population and scaled
so that one base additive standard deviation equals 0.707 t ha⁻¹ (base additive
variance 0.5). The elite pool is then produced by a burn-in of recurrent selection,
which necessarily depletes part of that variance; the realised elite additive variance
at the start of the comparison is therefore lower than 0.5 and is recorded per
replicate (`eliteVarObs`). Variance retention is reported as a percentage of the
*elite* base, which is the quantity relevant to a breeding programme's remaining
selection potential.

## Licence

Code released under the MIT Licence. Please cite the paper above if you use it.
