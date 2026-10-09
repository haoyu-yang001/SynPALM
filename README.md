# SynPALM: Synthetic Phenotype Assisted Linear Mixed Models

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21782545.svg)](https://doi.org/10.5281/zenodo.21782545)

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
| ranger | 0.18.0 | from CRAN; random forests for the synthetic phenotype |
| BEDMatrix | 2.0.4 | optional, from CRAN; reading PLINK `.bed` files |
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

A self-contained demonstration is bundled with the package. It simulates a
cohort of 20,000 individuals in 5,000 four-member families, measures the protein
in only 10% of them, offers the random forest 200 continuous and 30 binary
candidate surrogates (of which only a handful carry information), and tests 500
variants of which the first is causal. No external data are required.

### Instructions to run

```r
source(system.file("examples", "quickstart.R", package = "SynPALM"))
```

Simulation settings are collected at the top of that file and can be edited
freely. The random seed is fixed so that the output below is reproducible.

### Expected output

```
Cohort: 20000 individuals | 1999 with a measured protein | 200 continuous + 30 binary candidate surrogates

--- Synthetic phenotype (5-fold cross-fit) ---
  correlation with the measured protein, per fold and overall:
 fold n_train n_labelled n_predicted   rho rho_spearman    r2
    1    1595        404        4000 0.898        0.893 0.806
    2    1584        415        4000 0.916        0.912 0.839
    3    1603        396        4000 0.902        0.902 0.813
    4    1635        364        4000 0.881        0.874 0.775
    5    1579        420        4000 0.899        0.893 0.808
  all      NA       1999       20000 0.900        0.896 0.810
  informative surrogates selected in all 5 folds   : 10 of 10

--- First rows of fit$results ---
    variant    af n_missing SynSurrG_beta SynSurrG_se SynSurrG_p
1 rs_demo_1 0.286         0        0.1025      0.0166   7.01e-10
2 rs_demo_2 0.316         0       -0.0173      0.0163   2.87e-01
3 rs_demo_3 0.236         0        0.0137      0.0179   4.45e-01
  SynSurrG_log10p ObsG_beta ObsG_se   ObsG_p ObsG_log10p rho_oof
1           9.155    0.1379  0.0321 1.78e-05       4.749     0.9
2           0.542   -0.0206  0.0317 5.16e-01       0.288     0.9
3           0.351    0.0602  0.0347 8.30e-02       1.081     0.9

--- Causal variant (true standardised effect = 0.08) ---
               method log10p   beta     se
   SynPALM (SynSurrG)  9.155 0.1025 0.0166
 observed only (ObsG)  4.749 0.1379 0.0321

--- Calibration on 499 null variants ---
  SynSurrG lambda = 1.047   P<0.05: 0.052   P<0.01: 0.012
  ObsG     lambda = 1.010   P<0.05: 0.052   P<0.01: 0.014

--- Variance components (SynSurrG) ---
  tau_T2   tau_TS   tau_S2 sigma_T2 sigma_TS sigma_S2 
  0.5549   0.4333   0.3228   0.3333   0.3694   0.5908 
  valid without projection: FALSE

--- Timing ---
synpalm_gwas() elapsed (s): 19.2
```

SynPALM recovers the causal variant far more strongly than the analysis of
measured individuals alone, while both stay calibrated on the null variants.
`beta` is a per-allele effect on the 0/1/2 genotype scale, whereas `BETA_G` in
the simulation is on the standardised scale, so they differ by `1 / sd(G)`.

### Expected run time

About 20 seconds on one core of a Linux cluster node (R 4.4.1), most of it in
the five random forests.

---

## Instructions for use — running SynPALM on your own data

### Inputs

| Argument | Type | Description |
|---|---|---|
| `protein` | named numeric vector | The measured protein. Names are sample IDs. It may list only the measured individuals; anyone absent or `NA` is unmeasured. |
| `rf_features` | data frame, rownames = IDs | Candidate surrogates for the random forest (labs, vitals, questionnaire items, ...). As many columns as you like: within each fold the `n_top` (default 100) most correlated with the protein are kept. |
| `rf_binary` | data frame, rownames = IDs | Optional 0/1 candidate surrogates (e.g. diagnosis categories), screened by a Wilcoxon test within each fold. |
| `rf_fixed` | character | Columns of `rf_features` always given to the forest (e.g. age, sex), not screened. |
| `covariates` | data frame, rownames = IDs | Adjustment covariates of the mixed model: age, sex, PCs, batch, ... Factors are expanded automatically. No missing values. |
| `grm` | sparse matrix, dimnames = IDs | Sparse genetic relatedness matrix, e.g. from FastSparseGRM. |
| `genotype` | `.bed` path or matrix | A PLINK `.bed` file (read with BEDMatrix; position and alleles come from the `.bim`), or any matrix-like object with sample IDs as rownames and variant IDs as colnames. |

Inputs are aligned by sample ID, never by position. The analysis set is the
individuals present in `grm`, `covariates`, `rf_features` and `rf_binary`;
everyone in it must also be in the genotype data. `protein` does not restrict
it: the unmeasured individuals are exactly what the synthetic phenotype adds. The individuals are reordered internally so that every
relatedness cluster is contiguous, which the block-wise algorithms need.

### One call

```r
library(SynPALM)

fit <- synpalm_gwas(protein     = protein,
                    covariates  = covariates,
                    rf_features = surrogates,
                    rf_binary   = diagnoses,
                    rf_fixed    = c("age", "sex"),
                    grm         = grm,
                    genotype    = "genotypes/chr22.bed")

head(fit$results)          # one row per variant
fit$accuracy               # synthetic vs measured protein, per fold and overall
```

`fit$results` has, per variant, `chr`, `pos`, `effect_allele` (the counted
allele), `other_allele`, `af`, `n_missing`, for each method `_beta`, `_se`,
`_p` and `_log10p`, and `rho_oof`. Missing genotypes are mean-imputed
(`2 * af`) before testing; `n_missing` counts them.

`fit$accuracy` measures how well the synthetic phenotype reproduces the
protein. The forest is 5-fold cross-fitted, so every measured individual is
predicted by a forest that never saw them or their relatives, and the
correlation is an honest out-of-sample one. The table has one row per fold
and a row `fold = "all"` pooling every measured individual (that pooled
Pearson correlation is `rho_oof`): `n_train`, `n_labelled`, `n_predicted`,
`rho` (Pearson), `rho_spearman` and `r2`. By default the methods are `SynSurrG` (SynPALM) and `ObsG`
(the same mixed model on measured individuals only); add `"SynSurr"` and
`"Obs"` to `methods` for the versions on one individual per relatedness
cluster.

### Step by step, for a genome-wide scan in parallel jobs

The null model depends only on the phenotype, so it is fitted once per protein
and reused by every genotype chunk:

```r
folds <- synpalm_folds(grm, K = 5)                       # relatedness-aware folds
pred  <- synpalm_predict(protein, surrogates, folds,     # cross-fitted random forest
                         rf_binary = diagnoses, rf_fixed = c("age", "sex"))
null  <- synpalm_null(pred$protein, pred$protein_hat,    # variance components
                      covariates, grm, folds)
saveRDS(null, "null_model.rds")

## then, in one job per chromosome:
null <- readRDS("null_model.rds")
res  <- synpalm_scan(null, sprintf("genotypes/chr%s.bed", chr))
```

`pred$accuracy` holds the cross-fitted correlations; save it with the null
model. `inst/examples/biobank_template.R` is a ready-to-edit script for exactly
this two-stage workflow (`fit` once per protein, writing
`<protein>_rf_accuracy.tsv`; `scan` once per chromosome).

### What happens inside

1. **Folds.** Whole relatedness clusters are assigned to five folds, so no
   individual's synthetic phenotype depends on a relative's measurement.
2. **Synthetic phenotype.** In each fold, surrogates are screened and a random
   forest (`ranger`, 300 trees) is trained on measured individuals outside the
   fold, then predicts everyone in the fold. Its correlation with the measured
   protein is reported per fold and overall.
3. **Null model.** The measured protein (measured individuals) and the
   synthetic phenotype (everyone) are inverse-normal transformed; variance
   components are estimated by Haseman-Elston regression on the sparse GRM.
4. **Score tests.** Variants are tested in chunks of 200.

By default the raw variance-component estimates are used, as in the
manuscript analysis (`vc_constraint = "none"`); `fit$null$pars$SynSurrG$vc`
records whether they formed a valid covariance. See `?SynPALM-package` for
`vc_constraint = "project"` and the other numerical-stability options.

### Lower-level interface

The functions behind the pipeline can be called directly. They take a single
list, `mydf`, with `X_all` (covariates, no intercept column), `S` (synthetic
phenotype, inverse-normal transformed), `Y_obs` (measured phenotype,
inverse-normal transformed, `NA` where unmeasured) and `GRM`, all in the same
order, with relatedness clusters contiguous. The GRM must be a general, not a
symmetric-class, sparse matrix.


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
| `R/` | Package source. `pipeline.R` holds the end-to-end interface (`synpalm_gwas()` and its steps); `synpalm_functions.R` the model-fitting and score-test functions; `SynPALM-package.R` imports and package-level documentation. |
| `inst/examples/quickstart.R` | The self-contained demo described above. |
| `inst/examples/biobank_template.R` | Template for a real biobank: fit once per protein, scan once per chromosome. |
| `tests/` | Unit tests on simulated data. |
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

Archived release used during peer review: Zenodo, DOI [10.5281/zenodo.21782545](https://doi.org/10.5281/zenodo.21782545) (tag `v0.1.1`).

Summary statistics are browsable at <https://syn-palm.genohub.org/>.
