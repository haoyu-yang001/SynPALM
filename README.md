# SynPALM: Synthetic Phenotype Assisted Linear Mixed Models

`SynPALM` is a robust and computationally scalable statistical framework for
proteome-wide GWAS in the presence of partially observed measurements.

## Introduction

The UK Biobank Pharma Proteomics Project (UKB-PPP) generated plasma proteomic
data for 54,219 participants. Proteomic measurements remain unavailable for
approximately 90% of the 500,000 UK Biobank participants, sharply limiting power
for genetic discovery. Conventional imputation can yield spurious associations
when the prediction model is misspecified. SynPALM offers a robust alternative.

SynPALM jointly analyses partially observed proteomic measurements and complete
synthetic proteomic data (predicted by machine learning) while accounting for
cryptic relatedness and population structure using mixed models. It is designed
to be:

- **Robust** — controls false positives even when the proteomic prediction model
  is misspecified.
- **Powerful** — gains statistical power as prediction accuracy improves.
- **Scalable** — handles UK Biobank–scale cohorts.

## Key features

- **Information recovery** — recovers association signals by modelling the joint
  distribution of the target and surrogate traits.
- **Mixed model integration** — accounts for genetic relatedness through linear
  mixed models.
- **Computational efficiency** — block-wise matrix inversion and Cholesky
  decomposition, exploiting the block structure of a sparse GRM.
- **Ablation support** — built-in comparators (observed-only, oracle, and
  independent-sample variants) for like-for-like comparison.

---

## System requirements

### Operating systems tested

| Platform | OS version | R version |
|---|---|---|
| macOS, Apple silicon | 15.x (Darwin 24.6.0) | 4.3.2 |
| Linux, Harvard FASRC cluster | Rocky Linux 8.10 | 4.3.1 |

Windows has not been tested. See the note on parallelisation below.

### Software dependencies

R (>= 4.3.0) and the following packages:

| Package | Version tested | Notes |
|---|---|---|
| Matrix | 1.6.1.1 | ships with R as a recommended package |
| dplyr | 1.1.4 | from CRAN |
| methods, parallel, stats, utils | — | ship with R |

No compilation is required; SynPALM contains only R code.

### Non-standard hardware

None. SynPALM runs on a standard desktop or laptop CPU, and no GPU is required.

Block-wise inversion is parallelised over two cores with
`parallel::mclapply`, which relies on process forking and therefore works on
macOS and Linux but not on Windows. On Windows the package will load, but the
inversion routines will not run.

Memory scales with cohort size and with the density of the genetic relatedness
matrix. The bundled demo (20,000 individuals, 500 variants) peaks at
approximately **583 MB**. The UK Biobank analysis in the manuscript
(N = 398,800, of whom 29,578 had observed protein measurements) was run on the Harvard FASRC cluster with 20 GB of memory per SLURM array task.

---

## Installation guide

```r
# install.packages("devtools")
devtools::install_github("haoyu-yang001/SynPALM")
```

**Typical install time on a normal desktop computer:** approximately 5 seconds.
There is nothing to compile. If `dplyr` is not already present, CRAN fetches it
and its dependencies first, which typically adds one to three minutes depending
on the connection. `Matrix` ships with R and is never downloaded.

Measured on an Apple silicon Mac, R 4.3.2, with dependencies already installed.

---

## Demo

A fully self-contained demonstration is bundled with the package. It simulates a
cohort of 20,000 individuals in 5,000 four-member families, discards 90% of the
target phenotypes, and tests 500 variants of which the first is causal. No
external data are required.

### Instructions to run

```r
source(system.file("examples", "quickstart.R", package = "SynPALM"))
```

Simulation settings are collected at the top of that file (`N_FAM`, `FAM_SIZE`,
`MISS_RATE`, `N_SNP`, `BETA_G`, `R_TS`, `TAU`, `SIGMA`) and can be edited
freely. The random seed is fixed so that the output below is reproducible.

### Expected output

```
GRM: 20000 x 20000  | class: dgCMatrix  | nonzero fraction: 2e-04
Cohort: 20000 individuals | 1999 with an observed target phenotype (10%), missingness 90%
Observed-vs-synthetic correlation: 0.582

--- structure of mydf ---
List of 4
 $ X_all:'data.frame':	20000 obs. of  7 variables:
 $ S    : num [1:20000] 1.1057 -0.9911 -0.1809 1.091 -0.0466 ...
 $ Y_obs: num [1:20000] 0.871 NA NA NA -1.291 ...
 $ GRM  :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots

Relatedness blocks found: 5000 (expected 5000)

Fitting the null model (step 1)...
Step 1 elapsed (s): 0.2

Running the SynPALM score test (step 2)...
Score test elapsed (s): 2.1  (500 variants, 4.2 ms per variant)

--- components returned by score_test_SynSurrG_multiply ---
[1] "T_score_SynSurrG"             "negative_log10_pval_SynSurrG"
[3] "hat_beta_SynSurrG"            "var_hat_beta_SynSurrG"

--- Causal variant (true standardised effect = 0.08) ---
   variant neglog10P hat_beta se_beta
 rs_demo_1     6.828   0.1475  0.0281

--- Calibration on 499 null variants ---
  median -log10(P)          : 0.305   (expected 0.301)
  genomic inflation lambda  : 1.024   (expected 1.000)
  proportion P < 0.05       : 0.0501   (expected 0.05)
  proportion P < 0.01       : 0.004   (expected 0.01)
  mean hat_beta             : -0.00159   (expected 0)

--- Timing ---
Whole demo elapsed (s): 3.5
Peak R memory (MB)    : 583
```

The four calibration figures are the substantive check. Under the null,
p-values should be uniform, so the median of `-log10(P)` should sit at 0.301 and
the genomic inflation factor at 1.0.

Note that `hat_beta` is a per-allele effect on the raw 0/1/2 genotype scale,
whereas `BETA_G` in the simulation is specified on the standardised genotype
scale. The two therefore differ by a factor of `1 / sd(G)`, which is
approximately 1.5 for the allele frequencies used here.

Elapsed times and peak memory will vary with hardware.

### Expected run time

Approximately **3.5 seconds** on an Apple silicon Mac (R 4.3.2), of which 0.2 s
is the null model fit and 2.1 s the 500 score tests.

---

## Instructions for use — running SynPALM on your own data

Every entry point takes a single list, conventionally called `mydf`, with four
elements:

| Element | Type | Description |
|---|---|---|
| `X_all` | data frame, N × p | Covariates for **all** individuals: age, sex, genotyping array, ancestry PCs, recruitment centre, and so on. **No intercept column** — the package adds one. |
| `S` | numeric, length N | Synthetic (predicted) phenotype, complete for all individuals. Inverse-normal transformed. |
| `Y_obs` | numeric, length N | Measured target phenotype, `NA` for individuals without a measurement. Inverse-normal transformed over the measured subset. |
| `GRM` | sparse matrix, N × N | Genetic relatedness matrix over all individuals. |

Rows of `X_all`, `S`, `Y_obs` and `GRM` must refer to the same individuals in
the same order. Align them by an explicit identifier such as `eid`, never by row
position.

### Two input requirements that are easy to miss

**The GRM must be a general sparse matrix, not a symmetric-class one.**
Internally `summary()` is used to read the (i, j, x) triplet. For a
symmetric-class matrix such as `dsCMatrix` that returns only the lower triangle,
which would yield incorrect relatedness blocks without raising an error.
Convert explicitly:

```r
GRM <- methods::as(methods::as(GRM, "CsparseMatrix"), "generalMatrix")
```

**Individuals must be ordered so that relatedness blocks are contiguous.** Block
detection sweeps the matrix once and cannot recover blocks whose members are
scattered across the ordering. Check the result:

```r
blocks <- find_blocks_vectorized(GRM)
length(blocks)                       # number of independent blocks
max(vapply(blocks, length, 1L))      # size of the largest block
```

### Worked pipeline

```r
library(SynPALM)
library(Matrix)

## 1. one individual per relatedness block, for the comparators that assume
##    independence
blocks <- find_blocks_vectorized(GRM)
independent_indices <- vapply(blocks,
                              function(b) b[sample.int(length(b), 1)],
                              integer(1))
obs_protein_index <- which(!is.na(mydf$Y_obs))

## 2. step 1 -- variance components and null model quantities.
##    Depends only on the phenotype, so it is fitted once per protein and
##    reused across every genotype chunk.
SynSurrG_step1 <- SynSurrG_ablation_estimate(mydf)

## 3. genotypes, read in chunks to bound memory
Gmat <- fix_constant_columns(
  Gmat, intersect(obs_protein_index, independent_indices))

## 4. step 2 -- score test across the chunk
res <- score_test_SynSurrG_multiply(g_matrix   = Gmat,
                                    step1_pars = SynSurrG_step1)

## 5. results, one element per variant
res$negative_log10_pval_SynSurrG
res$hat_beta_SynSurrG
res$var_hat_beta_SynSurrG
```

In the manuscript analysis variants were processed in chunks of 200 and the
per-chunk results combined with `merge_results()`.

### Comparator analyses

The same two-step pattern applies to the methods SynPALM is compared against.
Note which ones require `independent_indices`:

| Analysis | Step 1 | Score test |
|---|---|---|
| SynPALM, mixed model | `SynSurrG_ablation_estimate(mydf)` | `score_test_SynSurrG_multiply(g, pars)` |
| Observed only, mixed model | `ObsG_ablation_estimate(mydf)` | `score_test_ObsG_multiply(g, pars)` |
| SynPALM, independent samples | `SynSurr_ablation_estimate(mydf, idx)` | `score_test_SynSurr_multiply(g, pars, idx)` |
| Observed only, independent samples | `Obs_ablation_estimate(mydf, idx)` | `score_test_Obs_multiply(g, pars, idx)` |
| Oracle, independent samples | `Oracle_ablation_estimate(mydf, idx)` | `score_test_Oracle_multiply(g, pars, idx)` |

### Expected wall time at scale

In the manuscript analysis, covering 591,558 directly genotyped variants in
398,800 individuals, variance component estimation (step 1) took 3.14 minutes on
a single CPU core, and each subsequent per-variant association test took 0.22
seconds. Because step 1 is computed once per protein and the SNP-invariant matrix
quantities are shared across variants, parallelising the per-variant tests across
100 threads reduced the total runtime to approximately 21 minutes per protein.

The per-variant cost of 4.2 ms in the demo reflects its much smaller cohort.

---

## Reproducing the manuscript results

The UK Biobank analysis reported in the manuscript cannot be re-run outside an
approved compute environment. It requires:

- individual-level UK Biobank phenotype, proteomic and genotype data, available
  only under approved access (application 52008);
- a pre-computed sparse GRM for the full cohort;
- an HPC cluster. The analysis was run on the Harvard FASRC cluster as an array
  of SLURM jobs, one per protein and genotype chunk.

None of these can be redistributed, so the analysis is documented here rather
than packaged as a runnable example.

`reproduce/ukb_ppp_analysis.R` is the driver script used to produce the
proteome-wide results. It is provided so that reviewers can inspect the exact
call sequence, input formats, covariate set, chunking scheme and output
structure. All file paths are collected in a configuration block at the top of
the script; running it would require substituting paths valid within the
reader's own approved environment.

For a self-contained, runnable check of the method and its calibration, see the
[Demo](#demo) section above, which requires no external data.

---

## Repository contents

| Path | Contents |
|---|---|
| `R/` | Package source. `synpalm_functions.R` holds the analysis functions; `SynPALM-package.R` holds imports and package-level documentation. |
| `inst/examples/quickstart.R` | The self-contained demo described above. |
| `reproduce/` | Driver script for the UK Biobank analysis, for inspection. |
| `man/` | Generated function documentation. |
| `tools/` | Development helpers, not part of the installed package. |

---

## License

MIT License. See [LICENSE](LICENSE) and [LICENSE.md](LICENSE.md).

## Citation

Yang, H., Wang, R., Song, S., and Lin, X. Synthetic Phenotype Assisted Linear
Mixed Models Improve Proteome-Wide Genetic Discovery in Incomplete Biobank Data.
Under review at *Nature Communications* (manuscript NCOMMS-26-058511-T).

A Zenodo DOI for the frozen release used during peer review will be added here.

Summary statistics are browsable at <https://syn-palm.genohub.org/>.
NA
release used during peer review.**
