## Revision copy of simulation_code/plot_for_point_robustness.R: reads simulation_code/results_rev1/ (rho = 0.75 rerun,
## other rho unchanged; see merge_power_bias_rev1.R) and writes to simulation_code/results_rev1/figs/.
library(dplyr)
library(ggplot2)
library(xtable)

#' Simulation
#' 
#' Evaluate performance of the synthetic surrogate estimator, compared
#' with oracle and marginal, as a function of the target-surrogate
#' correlation.
#' 
#' @param n0 Subjects with observed target outcomes.
#' @param reps Simulation replicates.
#' @param rho Target-surrogate correlation.
#' @return Data.frame.
Sim <- function(
    n0 = 1e3,
    reps = 100,
    rho = 0.00
) {
  
  # Proportion of variance explained by genotype.
  PVE_G <- 0.01
  
  # Missingness levels. miss = 0 included by default.
  MISS <- c(0.25, 0.50, 0.75, 0.90)
  
  # The outer loop runs over simulation replicates.
  # The inner loop runs over target-surrogate correlations.
  # Because the performance of the oracle and marginal estimators
  #   does not depend on correlation, these are estimated using
  #   the rho = 0 setting. 
  
  outer_loop <- lapply(seq_len(reps), function(i) {
    
    # Generate data with miss = 0.
    data <- DGP(n0 = n0, miss = 0, rho = rho, pve_g = PVE_G)
    
    # Oracle estimator.
    oracle_est <- OracleEst(data)
    
    # Observed data estimator.
    obs_est <- ObsEst(data)
    
    # Synthetic surrogate estimator.
    ss_est <- SynSurrEst(data)
    
    # Collect oracle and observed.
    ests <- rbind(oracle_est, obs_est, ss_est)
    
    out <- data.frame(
      miss = 0,
      rho = rho,
      method = c("oracle", "obs", "ss"),
      est = ests[, 1],
      se = ests[, 2],
      row.names = NULL
    )
    
    # Inner loop over non-zero correlation levels.
    ss_est <- lapply(MISS, function(miss) {
      
      data <- DGP(n0 = n0, miss = miss, rho = rho, pve_g = PVE_G)
      ss_est <- SynSurrEst(data)
      out <- c(
        miss = miss,
        rho = rho,
        est = as.numeric(ss_est[1]),
        se = as.numeric(ss_est[2])
      )
      
      return(out)
    })
    ss_est <- do.call(rbind, ss_est)
    ss_est <- data.frame(ss_est) %>%
      tibble::add_column(method = "ss", .after = "rho")
    
    out <- rbind(out, ss_est)
    out$bias <- out$est - sqrt(PVE_G)
    return(out)
  })
  
  results <- do.call(rbind, outer_loop)
  return(results)
}


# -----------------------------------------------------------------------------


#' Data Generating Process
#' 
#' Generate data. "g" is genotype, "x" is a covariate, "y" is the outcome, and
#' "yobs" is the outcome after introduction of missingness, "s" is the
#' surrogate.
#'
#' The number of subjects with observed target outcomes is always n. If
#' missingness is > 0, then subjects with target outcomes missing are added to
#' achieve the desired level of missingness.
#' 
#' @param n0 Subjects with observed target outcomes.
#' @param miss Missingness rate.
#' @param rho Target-surrogate correlation.
#' @param maf Minor allele frequency.
#' @param pve_g Proportion of variation explained by genotype.
#' @param pve_x Proportion of variation explained by covariate.
#' @return Data.frame.
DGP <- function(
    n0,
    miss,
    rho,
    maf = 0.25,
    pve_g = 0.005,
    pve_x = 0.10
) {
  
  # Subjects with observed target outcomes.
  n1 <- round(n0 * miss / (1 - miss))
  n <- n0 + n1
  
  # Genotype and covariate.
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  g <- scale(g)
  x <- stats::rnorm(n = n)
  x <- scale(x)
  
  # Linear predictor.
  beta <- c(sqrt(pve_g), sqrt(pve_x))
  eta <- as.numeric(cbind(g, x) %*% beta)
  
  # Residuals.
  sigma <- 1 - sum(c(pve_g, pve_x))
  eps_s <- stats::rnorm(n, sd = sqrt(sigma))
  eps_y <- rho * eps_s + sqrt(1 - rho^2) * stats::rnorm(n, sd = sqrt(sigma))
  
  # Oracle outcome.
  s <- eta + eps_s
  y <- eta + eps_y
  
  # Outcome with missingness.
  draw <- sample(seq_len(n), size = round(miss * n), replace = FALSE)
  yobs <- y
  yobs[draw] <- NA
  
  out <- data.frame(
    g = g,
    x = x,
    s = s,
    y = y,
    yobs = yobs
  )
  return(out)
}


#' Estimate Genetic Effect
#' 
#' Estimates the genetic effect and standard error.
#' 
#' @param covar Data.frame of covariates. Should not include "g" or "y".
#' @param geno Genotype vector.
#' @param outcome Outcome vector.
#' @return Numeric vector containing the estimate "est" and standard error "se".
EstBetaG <- function(covar, geno, outcome) {
  df <- data.frame(
    y = outcome,
    g = geno,
    x = covar
  ) %>% na.omit()
  
  fit <- lm(y ~ g + x, data = df)
  fit_summary <- summary(fit)$coefficients
  
  out <- fit_summary[
    rownames(fit_summary) == "g", c("Estimate", "Std. Error")]
  names(out) <- c("est", "se")
  return(out)
}


#' Oracle Estimator
#' 
#' @param data Data.frame.
#' @return Numeric vector containing "oracle_est" and "oracle_se".
OracleEst <- function(data) {
  oracle_est <- EstBetaG(
    covar = data$x,
    geno = data$g,
    outcome = data$y
  )
  names(oracle_est) <- paste0("oracle_", names(oracle_est))
  return(oracle_est)
}


#' Observed Data Estimator
#' 
#' @param data Data.frame.
#' @return Numeric vector containing "obs_est" and "obs_se".
ObsEst <- function(data) {
  obs_est <- EstBetaG(
    covar = data$x,
    geno = data$g,
    outcome = data$yobs
  )
  names(obs_est) <- paste0("obs_", names(obs_est))
  return(obs_est)
}


#' SynSurr Estimator
#' 
#' @param data Data.frame.
#' @param imp_param Imputation parameters.
#' @return Numeric vector containing "ss_est" and "ss_se".
SynSurrEst <- function(data) {
  
  # SynSurr estimator.
  fit <- SurrogateRegression::FitBNR(
    t = data$yobs,
    s = data$s, 
    X = data[, c("g", "x")],
  )
  param <- coef(fit, type = "Target")
  
  out <- c(
    ss_est = param$Point[1],
    ss_se = param$SE[1]
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Plot Simulation Results
#'
#' Violin + box plots of the estimates by method: Obs-GWAS (orange; all of its
#' estimates, as before) and Syn-PALM at each missing rate (light to dark blue),
#' one panel per rho. Common look: theme_rev1.R.

source("<PROJECT_DIR>/simulation_code/rev1/theme_rev1.R")

PlotSim <- function(sim, base_size = 12, ncol = NULL) {
  df <- sim %>% dplyr::filter(method != "oracle")
  miss_ss <- sort(unique(df$miss[df$method == "ss"]))
  lev <- c("obs", paste0("ss_", miss_ss))
  df$grp <- factor(ifelse(df$method == "obs", "obs", paste0("ss_", df$miss)), levels = lev)
  fills <- c(obs = method_cols[["Obs-GWAS"]],
             setNames(unname(miss_cols[as.character(miss_ss)]), paste0("ss_", miss_ss)))
  xlabs <- c(obs = "Obs-\nGWAS", setNames(sprintf("%g%%", 100 * miss_ss), paste0("ss_", miss_ss)))
  ## label over the Syn-PALM violins (positions 2 .. n), in the headroom of each panel
  n <- length(lev)
  grp_lab <- data.frame(x = (2 + n) / 2, y = Inf, lab = "Syn-PALM, by missing rate")
  ggplot(df, aes(grp, bias)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
    geom_vline(xintercept = 1.5, colour = "grey70", linewidth = 0.35) +
    geom_violin(aes(fill = grp), trim = FALSE, colour = "grey25", linewidth = 0.3) +
    geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", colour = "grey15", linewidth = 0.3) +
    geom_text(data = grp_lab, aes(x, y, label = lab), inherit.aes = FALSE, vjust = 1.5,
              size = 0.8 * base_size / .pt, colour = "grey20") +
    facet_wrap(~ rho, ncol = ncol, labeller = rho_labeller) +
    scale_fill_manual(values = fills, guide = "none") +
    scale_x_discrete(labels = xlabs) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = NULL, y = "Bias") +
    theme_rev1(base_size, legend = "none") +
    theme(panel.grid.major.x = element_blank())
}

# -----------------------------------------------------------------------------

kinsetting <- "siblings"

GRMr = 0.5

miss_vec <- c(0,0.25,0.5,0.75,0.9)
rho_vec <- c(0,0.75)
snp_heritability_vec <- 0

settings <- expand.grid(miss = miss_vec, rho = rho_vec, pve_g = snp_heritability_vec)

res_tab_list <- list() # 先初始化一个空list

for (rr in 1:nrow(settings)) {
  
  miss <- settings[rr,]$miss
  rho <- settings[rr,]$rho
  pve_g <- settings[rr,]$pve_g
  
  output_path <- paste0("<PROJECT_DIR>/simulation_code/results_rev1/point_estimate/GRMr_",
                        GRMr,"_missing_",miss,"_rho_",rho,"_rhogeno_",kinsetting,"_her_",pve_g, "_negative_log10_pvals_for_power.rds")
  
  aa <- as.data.frame(readRDS(output_path))
  
  res_tab <- data.frame(
    miss = rep(miss, 3 * nrow(aa)),
    rho = rep(rho, 3 * nrow(aa)),
    method = rep(c("oracle", "obs", "ss"), each = nrow(aa)),
    bias = c(aa$OracleG_hat_beta, aa$ObsG_hat_beta, aa$SynSurrG_hat_beta),
    se = c(sqrt(aa$OracleG_var_hat_beta), sqrt(aa$ObsG_var_hat_beta), sqrt(aa$SynSurrG_var_hat_beta))
  )
  
  res_tab_list[[rr]] <- res_tab # 存到list
}

# 合并所有结果
all_res_tab <- do.call(rbind, res_tab_list)
rownames(all_res_tab) <- NULL



# -----------------------------------------------------------------------------

## Robustness_Syn-PALM.pdf: rho = 0 above rho = 0.75; fig2a_robustness_rev1.pdf (Fig 2A of the
## main text): the two side by side
q <- PlotSim(all_res_tab, ncol = 1)
ggsave(plot = q, filename = "<PROJECT_DIR>/simulation_code/results_rev1/figs/Robustness_Syn-PALM.pdf",
       device = cairo_pdf, height = 5.6, width = 6, units = "in")

q_row <- PlotSim(all_res_tab, ncol = 2)
ggsave(plot = q_row, filename = "<PROJECT_DIR>/simulation_code/results_rev1/figs/fig2a_robustness_rev1.pdf",
       device = cairo_pdf, height = 2.45, width = 7, units = "in")
