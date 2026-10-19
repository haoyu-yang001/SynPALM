# Code for the revised manuscript

This folder contains the code used for every simulation, analysis, table and figure of the
revised manuscript and Supplementary Information. The SynPALM package in the repository root
implements the method; these scripts are the study code that calls it (through
`common/funcs_v2.R`, the analysis functions the package was built from, and the vectorised score
test `common/funcs_fastscan.R`).

**No individual-level data are included.** The UK Biobank analyses need approved access to the
UK Biobank phenotype, proteomic and genotype data and a pre-computed sparse GRM; those scripts
are provided so that every step can be inspected. All simulations generate their own data and
can be rerun anywhere; every random draw is seeded and the seed is stored with the output.

Site-specific paths are replaced by placeholders:

| Placeholder | Meaning |
|---|---|
| `<PROJECT_DIR>` | project root with `data/`, `genotype_data/`, `results/`, `simulation_code/`, `code/` |
| `<GRM_DIR>` | directory with the pre-computed UK Biobank sparse GRM |
| `<UKB_DIR>` | UK Biobank data directory (UKB-PPP summary statistics used for the comparison with Sun et al.) |
| `<LAB_DIR>/protein_map` | table of protein gene positions (GRCh38) and Olink panels |
| `<INTEL_ONEAPI_DIR>` | Intel oneAPI installation (MKL BLAS, loaded with `LD_PRELOAD`) |

Jobs were run with SLURM; `submit_*.sh` files are the array submissions and record the resources
used. All scripts skip outputs that already exist, so interrupted arrays can be resubmitted.

## Folders

| Folder | Content |
|---|---|
| `common/` | `funcs_v2.R` (model, method-of-moments variance components, score tests, simulation data-generating functions), `funcs_fastscan.R` (vectorised Syn-PALM score test, identical to the per-SNP version up to 1e-13), `Fun.R` (POP-GWAS) |
| `simulations/` | simulations of the revision (`*_rev1.R`), their submission scripts, analysis and plotting scripts (common figure style: `theme_rev1.R`); `sim_dgp_rev1.R` is the general data-generating model of the new simulations; `power_bias_rho075_rev1.R` and `true_geno_power_rho075_rev1.R` run every $\rho$ (argument `rho`) |
| `ukb_ppp/pipeline/` | UKB-PPP analysis: sample and predictor preparation, five-fold cross-fitting folds by relatedness cluster, step 1 (cross-fitted random-forest synthetic phenotype), step 2 (genome-wide scan), split-sample and permutation analyses, and their SLURM scripts |
| `ukb_ppp/post/` | per-protein post-processing, cross-protein summaries, external validation, fine-mapping (SuSiE), SNPsea input, *trans* hubs, stratified analyses, predictor tables |
| `ukb_ppp/figs/` | figure panels of Figs 4 and 5 and supplementary UKB-PPP figures |
| `height_ablation/` | height ablation analysis with cluster-based cross-fitting |

## Figures and tables

| Result | Script(s) |
|---|---|
| Fig. 2A (bias / robustness) | `simulations/power_bias_rho075_rev1.R`, `merge_power_bias_rev1.R`, `plot_for_point_robustness_rev1.R` |
| Fig. 2B, Table 1, Supp. Figs S2–S4 (type I error) | `simulations/typeIerror_rev1.R`, `analyze_t1e_rev1.R`, `plot_fig2b_rev1.R`, `table1_tex_rev1.py` |
| Fig. 2C, Fig. 3, Supp. Figs S5–S7 (power) | `simulations/power_bias_rho075_rev1.R`, `power_analysis_kinsetting_plots_rev1.R` |
| Supp. Section S3.2, Table S3, Fig. S8 (SNPs affecting only the synthetic phenotype) | `simulations/alphaG_rev1.R`, `analyze_alphaG_vc_rev1.R` |
| Table S1, Supp. Section S3.3, Fig. S9 (genetic vs residual sharing) | `simulations/rep_sims_rev1.R` (S36), `analyze_rep_rev1.R`, `plot_sharing_rev1.R` |
| Supp. Section S3.4, Table S4 (proportion of related individuals) | `simulations/rep_sims_rev1.R` (S4), `analyze_rep_rev1.R` |
| Supp. Section S3.5, Table S5 (missing at random) | `simulations/rep_sims_rev1.R` (S5), `analyze_rep_rev1.R` |
| Supp. Section S3.6, Table S6 (relatives split between labelled and unlabelled sets; POP-GWAS) | `simulations/typeIerror_split_rev1.R`, `analyze_t1e_rev1.R` (`T1E_DIR`) |
| Supp. Section S5, Table S10, Figs S11–S12 (UK Biobank genotypes) | `simulations/true_geno_typeIerror_calgeno_rev1.R`, `true_geno_power_rho075_rev1.R`, `analyze_true_geno_rev1.R` |
| Supp. Section S6, Figs S14–S15 (variance components, MoM vs REML vs BOLT-REML) | `simulations/variance_components_rev1.R`, `analyze_alphaG_vc_rev1.R` |
| Height ablation (main-text table and Supp. Section S4) | `height_ablation/*.R` |
| UKB-PPP headline numbers, Table S13 (gain by prediction accuracy) | `ukb_ppp/post/aggregate_all.R`, `headline_numbers.R`, `category_and_rho_strata.R` |
| External validation (GWAS Catalog, Sun et al.), additional associations | `ukb_ppp/post/external_validation.R` |
| Fig. 4 | `ukb_ppp/figs/fig4_panels_cf.R`, `fig4c_top20_cf.R`, `fig4f_snpsea.R` (SNPsea: `ukb_ppp/post/snpsea_prep.R`, `ukb_ppp/pipeline/run_snpsea.sh`) |
| Fig. 5 | `ukb_ppp/figs/fig4_panels_cf.R` (panel A), `ukb_ppp/post/finemap_susie.R` + `ukb_ppp/figs/fig5b_circos_susie.R` (panel B), `ukb_ppp/figs/manhattan_protein_cf.R` (panel C) |
| *Trans* hubs | `ukb_ppp/post/trans_hubs.R` |
| Split-sample replication, Table S12 | `ukb_ppp/pipeline/split_validation_crossfit.R`, `ukb_ppp/post/summarize_split.R` |
| Permutation QQ plot, Fig. S18 | `ukb_ppp/pipeline/step2_protein_fast_jointperm.R`, `ukb_ppp/post/qq_permuted.R` |
| Prediction accuracy, Fig. S16, Table S11 | `ukb_ppp/figs/prediction_r2_figure.R`, `ukb_ppp/post/aggregate_all.R` |
| Predictors (Supplementary Data) | `ukb_ppp/post/predictor_tables.R` |
| LILRB5 conditional analysis, Supp. Section S7.3 | `ukb_ppp/pipeline/submit_step2_cond.sh` (`SYNPALM_COND_SNPS` in `step2_protein_fast.R`) |

## Order of the UKB-PPP analysis

1. `pipeline/data_selection.R` (sample and predictors), `pipeline/make_folds_once.R` (relatedness-cluster folds).
2. `pipeline/submit_step1.sh` → `step1_protein_crossfit.R`: cross-fitted synthetic phenotype per protein.
3. `pipeline/submit_step2.sh` → `step2_protein_fast.R`: genome-wide scan per protein (Syn-PALM, Obs-GWAS, SynSurr, Obs).
4. `pipeline/post_parallel_local.sh` → `post/post_protein.R`: significant variants, loci and summary statistics per protein.
5. `pipeline/final_all.sh`: cross-protein summaries, external validation, fine-mapping, figure panels.
