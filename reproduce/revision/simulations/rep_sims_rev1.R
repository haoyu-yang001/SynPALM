## Replicate simulations of the revision on the general DGP of sim_dgp_rev1.R:
##   S36 : genetic vs residual sharing at fixed total rho (Reviewer 2 comment 8),
##         including pi = 0 (fully observed; Supp Table S1, Reviewer 1 comment 3,
##         Reviewer 3 minor 1)
##   S4  : proportion of related individuals, whole-family vs individual-level
##         labelling, and no related pairs at all (Reviewer 2 comment 3)
##   S5  : missingness at random given the covariate x ("age") or given S
##         (Reviewer 2 comments 1 and 9); also null SNPs for type I error
## Per replicate one phenotype pair is drawn with n_snp SNPs that each explain
## her of the variance of Y (beta_G) and of S (alpha_G = beta_G, as in the
## submitted power simulations), plus n_null null SNPs; all methods are fitted
## once and score every SNP. MAF 0.25; tauT = 0.7, tauS = 0.4, sigmaT = 0.7,
## sigmaS = 0.5 (SDs).
##   Rscript rep_sims_rev1.R <study> <setting> <array> [n_rep=25]
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R")
                   source("simulation_code/rev1/sim_dgp_rev1.R") })
options(synsurrg.verbose = FALSE)
args     <- commandArgs(trailingOnly = TRUE)
study    <- args[1]
setting  <- as.integer(args[2])
array_id <- as.integer(args[3])
n_rep    <- if (length(args) >= 4) as.integer(args[4]) else 25L
maf <- 0.25; n_obs0 <- 3e4

## ---- settings -----------------------------------------------------------------
S <- switch(study,
  S36 = rbind(
    ## total rho 0.5: rho_g from 0 (residual sharing only) to 0.99
    expand.grid(rho_tot = 0.5,  rho_g = c(0, 0.25, 0.5, 0.75, 0.99), miss = c(0, 0.5, 0.9)),
    ## total rho 0.25: rho_g from 0 to 0.566 (genetic sharing only, rho_e = 0)
    expand.grid(rho_tot = 0.25, rho_g = c(0, 0.15, 0.3, 0.45, 0.566), miss = c(0, 0.5, 0.9)),
    ## submitted 3/5 : 2/5 split, fully observed (new Supp Table S1)
    data.frame(rho_tot = c(0, 0.25, 0.5, 0.75), rho_g = NA, miss = 0)),
  S4 = expand.grid(prop_related = c(0, 0.1, 0.25, 0.5, 1), label = c("pairs", "random"),
                   miss = 0.5, rho_tot = 0.5, stringsAsFactors = FALSE),
  S5 = expand.grid(label = c("random", "mar_x", "mar_s"), miss = c(0.5, 0.9), rho_tot = 0.5,
                   stringsAsFactors = FALSE))
st <- S[setting, , drop = FALSE]
g <- function(k, d) if (!is.null(st[[k]]) && !is.na(st[[k]])) st[[k]] else d
miss  <- g("miss", 0.5); prop <- g("prop_related", 1); lab_type <- g("label", "pairs")
vc <- if (is.na(g("rho_g", NA))) vc_from_total(g("rho_tot", 0.5)) else {
  re <- (g("rho_tot", 0.5) * sqrt(0.98 * 0.41) - g("rho_g", NA) * 0.7 * 0.4) / (0.7 * 0.5)
  vc_from_rho(g("rho_g", NA), re) }
pve_x  <- if (study == "S5") 0.10 else 0
her    <- 5e-4
n_snp  <- 20L
n_null <- if (study == "S5") 200L else 0L
n_all  <- round(n_obs0 / (1 - miss))
coh <- make_cohort(n_all, prop_related = prop, r = 0.5)

OUT <- sprintf("simulation_code/results_rev1/%s/", study)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("%ssetting_%02d_array_%03d.rds", OUT, setting, array_id)
if (file.exists(out_file)) quit(save = "no")
seed <- 20261016L + 1000000L * match(study, c("S36", "S4", "S5")) + 1000L * setting + array_id
set.seed(seed)
cat(study, "setting", setting, ":", paste(names(st), unlist(st), collapse = ", "),
    sprintf("| rho_g %.3f rho_e %.3f rho %.3f | n %d | seed %d\n", vc$rho_g, vc$rho_e, vc$rho, n_all, seed))

V_T <- vc$Gm[1, 1] + vc$Em[1, 1]
beta_G <- sqrt(her * V_T / (2 * maf * (1 - maf)))
beta_X <- sqrt(pve_x * V_T)
rows <- list()
for (k in seq_len(n_rep)) {
  G  <- gene_drop_cohort(coh, n_snp + n_null, maf)
  x  <- rnorm(n_all)
  e  <- draw_noise(coh, vc)
  b  <- c(rep(beta_G, n_snp), rep(0, n_null))
  Y  <- as.numeric(G %*% b) + beta_X * x + e[, 1]
  Sv <- as.numeric(G %*% b) + beta_X * x + e[, 2]
  lab <- label_set(coh, miss, lab_type, x = x, s = Sv)
  mydf <- make_mydf(coh, Y, Sv, x, lab)
  ind <- one_per_family(coh)
  syn <- add_fastscan_cache(SynSurrG_typeIerror_step1(mydf))
  ssr <- add_fastscan_cache(SynSurr_typeIerror_step1(mydf, ind))
  r <- list(
    SynSurrG = score_test_SynSurrG_fast(G, syn),
    ObsG     = score_test_ObsG_multiply(G, ObsG_typeIerror_step1(mydf)),
    OracleG  = score_test_OracleG_multiply(G, OracleG_typeIerror_step1(mydf)),
    SynSurr  = score_test_SynSurrG_fast(G, ssr, subset = ind, label = "SynSurr"),
    Obs      = score_test_Obs_multiply(G, Obs_typeIerror_step1(mydf, ind), ind))
  for (m in names(r)) rows[[length(rows) + 1]] <- data.frame(
    rep = k, method = m, snp = seq_along(b), beta = b,
    hat_beta = r[[m]][[paste0("hat_beta_", m)]], var_hat_beta = r[[m]][[paste0("var_hat_beta_", m)]],
    nlp = r[[m]][[paste0("negative_log10_pval_", m)]], n_lab = length(lab),
    tau_T2 = syn$vc$used[["tau_T2"]], tau_TS = syn$vc$used[["tau_TS"]])
}
out <- do.call(rbind, rows)
attr(out, "settings") <- c(as.list(st), list(study = study, rho = vc$rho, rho_g = vc$rho_g, rho_e = vc$rho_e,
                           Gm = vc$Gm, Em = vc$Em, her = her, beta_G = beta_G, pve_x = pve_x,
                           n_all = n_all, n_pair = coh$n_pair, seed = seed, array = array_id))
saveRDS(out, out_file)
cat("saved", out_file, "\n")
