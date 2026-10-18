## Supplementary true-genotype (UKB) results of the revision: type I error table
## (Supp Table "typeIcontrol_ture") and power figures (Supp Figs
## power_Syn-PALM_ObsG_true_geno, power_ratio_Syn-PALM_vs_{obsG,SynSurr}_true_geno).
## Same computations, table layout and plot code as
## typeIerror_control_true_geno_analysis.R and power_true_geno_analysis.R; runs on
## partial results and says what is missing (incomplete rows are flagged).
##
## Type I error (all four rho): true_geno_typeIerror_calgeno_rev1.R, i.e. every
##   autosomal genotyped SNP (575,236) scored under 17 seeded null-phenotype
##   replicates (40 SNP arrays per replicate) -> ~9.78e6 Syn-PALM tests per
##   (missing rate, rho); results_rev1/true_geno_typeIerror_calgeno/. This replaces
##   the submitted imputed-variant rows (results/true_geno_error_rate_results/chr_*,
##   whose script and genotype files no longer exist); those are written alongside
##   for reference (rho < 0.75 only).
## Power: results_rev1/true_geno_power_rho<100 rho>/ (true_geno_power_rho075_rev1.R; the
##   autosomal SNPs of the with_seed(241150) sample). rho < 0.75 are rerun as well
##   (2026-10-16); until a (missing rate, rho, pve_g) rerun is complete the submitted
##   result is used (results/true_geno_power_results/, 20,000 rows in array order,
##   non-autosomal rows dropped) and the source is reported.
##
##   Rscript simulation_code/rev1/analyze_true_geno_rev1.R
## Outputs (simulation_code/results_rev1/figs/):
##   typeIerror_SynG_true_geno_rev1.csv          (columns as table_results/typeIerror_SynG_true_geno.csv)
##   typeIerror_SynG_true_geno_rev1_rows.tex     (LaTeX rows of the supplementary table, alpha 1e-2..1e-5)
##   typeIerror_SynG_true_geno_rev1_details.csv  (number of tests, replicates x arrays present, per row)
##   typeIerror_SynG_true_geno_submitted_imputed.csv (submitted rho < 0.75 rows, for reference)
##   power_Syn-PALM_ObsG_true_geno_rev1.pdf, power_ratio_Syn-PALM_vs_obsG_true_geno_rev1.pdf,
##   power_ratio_Syn-PALM_vs_SynSurr_true_geno_rev1.pdf, power_true_geno_rev1.csv
##   (then: cd simulation_code/rev1/figs_compose && pdflatex figS_true_geno_ratio_rev1.tex)
setwd("<PROJECT_DIR>")
suppressMessages({ library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(withr) })

R0  <- "simulation_code/results/"
R1  <- "simulation_code/results_rev1/"
FIG <- paste0(R1, "figs/"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CACHE <- paste0(R1, "true_geno_analysis/"); dir.create(CACHE, showWarnings = FALSE)
arraynum <- 1000                         # power: SNP sample split into 1000 arrays
n_rep_t1e <- 17; n_arr_t1e <- 40         # type I error: replicates x SNP arrays
miss_vec <- c(0, 0.25, 0.5, 0.75, 0.9)
rho_vec  <- c(0, 0.25, 0.5, 0.75)
rho_old  <- c(0, 0.25, 0.5)
alpha    <- c(10^-2, 10^-3, 10^-4, 10^-5, 10^-6)
snp_heritability_vec <- seq(0.0001, 0.0007, 0.0001)
missing_report <- character(0)
note <- function(...) { m <- paste0(...); missing_report <<- c(missing_report, m); message("NOTE: ", m) }

## SNP order of the genotype file used by the v2 scripts (BEDMatrix colnames = bim V2)
bim <- fread("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bim", header = FALSE, data.table = FALSE)
is_auto <- setNames(bim$V1 %in% as.character(1:22), bim$V2)
cat(sprintf("Genotyped SNPs: %d, non-autosomal (X/Y/XY/MT): %d (%.2f%%)\n",
            nrow(bim), sum(!is_auto), 100 * mean(!is_auto)))

########################## 1. type I error ###################################
count_t1e <- function(p) c(n = length(p), sapply(alpha, function(a) sum(10^(-p) < a)))
fmt_row <- function(cnt) as.numeric(format(cnt[-1] / cnt[1], scientific = TRUE, digits = 2))

## 1a. revision runs, all four rho (calibrated autosomal SNPs, 17 replicates)
T1E <- Sys.getenv("T1E_DIR", paste0(R1, "true_geno_typeIerror_calgeno/"))   # T1E_DIR: test runs only
grid_ra <- expand.grid(array = seq_len(n_arr_t1e), rep = seq_len(n_rep_t1e))
settings <- expand.grid(miss = miss_vec, rho = rho_vec)
tab <- NULL; det <- NULL
for (rr in seq_len(nrow(settings))) {
  miss <- settings$miss[rr]; rho <- settings$rho[rr]
  fs <- sprintf("%srep_%d_missing_%s_rho_%s_chunkarray_%d_negative_log10_pvals_type1error_rev1.rds",
                T1E, grid_ra$rep, miss, rho, grid_ra$array)
  have <- which(file.exists(fs))
  cc <- 0
  for (k in have) {
    x <- readRDS(fs[k])
    stopifnot(all(is_auto[attr(x, "snps")]))
    cc <- cc + count_t1e(x[, "SynSurrG_negative_log10_pvalues"])
  }
  if (length(have) < length(fs))
    note(sprintf("type I error miss = %s, rho = %s: %d of %d replicate x array files", miss, rho, length(have), length(fs)))
  tab <- rbind(tab, if (!length(have)) rep(NA, length(alpha)) else fmt_row(cc))
  det <- rbind(det, data.frame(miss, rho, n_tests = if (length(have)) cc[[1]] else 0,
                               files = length(have), files_expected = length(fs),
                               complete = length(have) == length(fs)))
}
tab <- cbind(settings, tab); colnames(tab) <- c("miss", "rho", alpha)
write.csv(tab, paste0(FIG, "typeIerror_SynG_true_geno_rev1.csv"), row.names = FALSE)
write.csv(det, paste0(FIG, "typeIerror_SynG_true_geno_rev1_details.csv"), row.names = FALSE)

## LaTeX rows in the layout of the submitted table (alpha = 1e-2 ... 1e-5)
e <- function(x) ifelse(is.na(x), "--", formatC(x, format = "e", digits = 1))
tex <- character(0)
for (rr in seq_len(nrow(tab))) {
  line <- paste0(paste(c(tab$miss[rr], tab$rho[rr], e(unlist(tab[rr, 3:6]))), collapse = "&"), "\\\\")
  if (rr == nrow(tab)) line <- paste0(line, "\\bottomrule")
  else if (tab$rho[rr + 1] != tab$rho[rr]) line <- paste0(line, "\\hline")
  if (!det$complete[rr]) line <- paste0(line, " % partial: ", det$files[rr], " of ", det$files_expected[rr], " files")
  tex <- c(tex, line)
}
writeLines(c("% Rows of Supp Table typeIcontrol_ture (analyze_true_geno_rev1.R): calibrated autosomal",
             sprintf("%% genotyped SNPs x %d null replicates; tests per setting: %s.", n_rep_t1e,
                     paste(unique(format(det$n_tests, big.mark = ",")), collapse = " / ")),
             tex), paste0(FIG, "typeIerror_SynG_true_geno_rev1_rows.tex"))
cat("\nType I error of Syn-PALM (calibrated autosomal SNPs, all rho):\n")
print(cbind(tab, n_tests = det$n_tests))

## 1b. submitted imputed-variant rows (rho < 0.75), for reference; cached counts
f_cache <- paste0(CACHE, "t1e_counts_submitted_chr.rds")
if (file.exists(f_cache)) cnt_chr <- readRDS(f_cache) else {
  cnt_chr <- list()
  for (rho in rho_old) for (miss in miss_vec) {
    fs <- sprintf("%strue_geno_error_rate_results/chr_%d_missing_%s_rho_%s_negative_log10_pvals_type1error.rds",
                  R0, 1:22, miss, rho)
    if (!all(file.exists(fs))) stop("submitted chr_ files missing for miss ", miss, " rho ", rho)
    cnt_chr[[paste(miss, rho)]] <- Reduce(`+`, lapply(fs, function(f) count_t1e(readRDS(f)[, "SynSurrG_negative_log10_pvalues"])))
  }
  saveRDS(cnt_chr, f_cache)
}
sub_set <- expand.grid(miss = miss_vec, rho = rho_old)
sub_tab <- cbind(sub_set, t(sapply(seq_len(nrow(sub_set)), function(i) fmt_row(cnt_chr[[paste(sub_set$miss[i], sub_set$rho[i])]]))),
                 n_tests = sapply(seq_len(nrow(sub_set)), function(i) cnt_chr[[paste(sub_set$miss[i], sub_set$rho[i])]][[1]]))
colnames(sub_tab)[3:7] <- alpha
write.csv(sub_tab, paste0(FIG, "typeIerror_SynG_true_geno_submitted_imputed.csv"), row.names = FALSE)

############################## 2. power ######################################
rs_all <- bim$V2
select_index <- with_seed(241150, sample(1:length(rs_all), 200*100, replace = F))
rs_sel <- rs_all[select_index]
auto_sel <- is_auto[rs_sel]
cat(sprintf("Power SNP sample: %d SNPs, %d non-autosomal (%.2f%%) dropped\n",
            length(rs_sel), sum(!auto_sel), 100 * mean(!auto_sel)))

POW_DIR <- function(rho) sprintf("%strue_geno_power_rho%03d/", R1, round(100 * rho))
settings <- expand.grid(miss = miss_vec, rho = rho_vec, pve_g = snp_heritability_vec)
com_results <- list(); pow_arr <- NULL
for (rr in seq_len(nrow(settings))) {
  miss <- settings$miss[rr]; rho <- settings$rho[rr]; pve_g <- settings$pve_g[rr]
  fs <- sprintf("%smissing_%s_rho_%s_her_%s_chunkarray_%d_negative_log10_pvals_for_power_rev1.rds",
                POW_DIR(rho), miss, rho, pve_g, 1:arraynum)
  have <- which(file.exists(fs))
  pow_arr <- rbind(pow_arr, data.frame(miss, rho, pve_g, arrays = length(have),
                                       source = if (rho == 0.75 || length(have) == arraynum) "rerun" else "submitted"))
  if (rho < 0.75 && length(have) < arraynum) {
    aa <- readRDS(paste0(R0, "true_geno_power_results/missing_", miss, "_rho_", rho, "_her_", pve_g,
                         "_negative_log10_pvals_for_power.rds"))
    if (nrow(aa) != length(rs_sel)) stop("unexpected row count in submitted power file (miss ", miss, ", rho ", rho, ")")
    aa <- aa[auto_sel, , drop = FALSE]
  } else {
    if (!length(have)) next
    aa <- do.call(rbind, lapply(fs[have], readRDS))
  }
  com_results[[rr]] <- data.frame(setting_num = rr, miss = miss, rho = rho, pve_g = pve_g, aa)
}
com_results <- bind_rows(com_results)
print(aggregate(arrays ~ rho + source, pow_arr, function(a) c(settings = length(a), min = min(a))))
for (m in miss_vec) {
  a <- pow_arr$arrays[pow_arr$miss == m & pow_arr$rho == 0.75]
  if (any(a < arraynum))
    note(sprintf("power rho = 0.75, miss = %s: %d-%d of %d arrays per pve_g (7 values)", m, min(a), max(a), arraynum))
}
if (any(pow_arr$source == "submitted")) note("power rho < 0.75: submitted results used for ",
                                              sum(pow_arr$source == "submitted"), " settings (rerun incomplete)")

color_ramp <- colorRampPalette(colors = c("white", "#0073C2FF"))
biv_palette <- color_ramp(n = 9)[c(3, 5, 7, 9)]
names(biv_palette) <- c("0", "0.25", "0.5", "0.75")
missing_names <- c(`0` = "Missing Rate = 0%", `0.25` = "Missing Rate = 25%",
                   `0.5` = "Missing Rate = 50%", `0.75` = "Missing Rate = 75%",
                   `0.9` = "Missing Rate = 90%")
rho_linetypes <- c("0" = "dotted", "0.25" = "dotdash", "0.5" = "dashed", "0.75" = "solid")
alpha_level <- 5e-08   # genome-wide level of the revision (was 1e-8 in the submission)

com_results$power_SynSurrG <- 10^{-com_results$SynSurrG_negative_log10_pval} < alpha_level
com_results$power_ObsG     <- 10^{-com_results$ObsG_negative_log10_pval} < alpha_level
com_results$power_OracleG  <- 10^{-com_results$OracleG_negative_log10_pval} < alpha_level
com_results$power_Obs      <- 10^{-com_results$Obs_negative_log10_pval} < alpha_level
com_results$power_Oracle   <- 10^{-com_results$Oracle_negative_log10_pval} < alpha_level
com_results$power_SynSurr  <- 10^{-com_results$SynSurr_negative_log10_pval} < alpha_level

power_tab <- com_results %>% group_by(miss, rho, pve_g) %>%
  summarise(n_snps = n(), across(starts_with("power_"), mean), .groups = "drop")
write.csv(power_tab, paste0(FIG, "power_true_geno_rev1.csv"), row.names = FALSE)

## ---- figures (common look: theme_rev1.R) ----------------------------------------
source("simulation_code/rev1/theme_rev1.R")
alpha_lab <- "5 %*% 10^-8"
## (i) power of Syn-PALM and Obs-GWAS (Supp Fig Syn-PALM_true_geno_power_vs_oracle)
pw <- com_results %>%
  group_by(miss, rho, pve_g) %>%
  summarise(`Syn-PALM` = mean(power_SynSurrG), `Obs-GWAS` = mean(power_ObsG), .groups = "drop") %>%
  pivot_longer(c(`Syn-PALM`, `Obs-GWAS`), names_to = "method", values_to = "power") %>%
  mutate(method = factor(method, c("Syn-PALM", "Obs-GWAS")), rho = factor(rho),
         panel = factor(miss_lab(miss), miss_lab(miss_vec)))
power_plot <- ggplot(pw, aes(pve_g * 100, power, colour = method, linetype = rho, shape = method,
                             group = interaction(rho, method))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_wrap(~ panel, nrow = 2) +
  scale_colour_manual(values = method_cols, name = "Method") +
  scale_shape_manual(values = c("Syn-PALM" = 16, "Obs-GWAS" = 17), name = "Method") +
  scale_linetype_manual(values = rho_lty, name = expression(rho)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
  scale_x_continuous(breaks = snp_heritability_vec * 100) +
  labs(x = "Per-SNP heritability (%)", y = expression(paste("Power at ", alpha == 5 %*% 10^-8))) +
  theme_rev1(11) +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.24), legend.box = "horizontal") +
  guides(linetype = guide_legend(order = 2, override.aes = list(colour = "grey20", shape = NA)),
         colour = guide_legend(order = 1), shape = guide_legend(order = 1))
ggsave(plot = power_plot, filename = paste0(FIG, "power_Syn-PALM_ObsG_true_geno_rev1.pdf"),
       device = cairo_pdf, height = 5.2, width = 7.6, units = "in")

## (ii) power ratios at SNP heritability 0.05% (Supp Fig Power_ratio_Syn-PALM_true_geno)
ind_pve_g <- 5
ratio_breaks <- function(limits) sort(unique(c(scales::pretty_breaks(n = 5)(limits), 1)))
ratio_data <- function(other) {
  com_results[, c("setting_num", "miss", "rho", "pve_g", "power_SynSurrG", other)] %>%
    group_by(miss, rho, pve_g) %>%
    summarise_at(vars(power_SynSurrG, all_of(other)), mean) %>%
    mutate(ratio = power_SynSurrG / .data[[other]]) %>%
    filter(abs(pve_g - snp_heritability_vec[ind_pve_g]) < 1e-12)
}
sub_txt <- bquote("Per-SNP heritability" ~ .(round(snp_heritability_vec[ind_pve_g] * 100, 2)) * "%," ~ alpha == 5 %*% 10^-8)

data_re <- ratio_data("power_ObsG")
power_ratio_plot_SynSurrG_vs_ObsG <- ggplot(data_re, aes(miss * 100, ratio, colour = factor(rho), group = rho)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_colour_manual(values = rho_cols, name = expression(rho)) +
  scale_x_continuous(breaks = miss_vec * 100) +
  scale_y_continuous(breaks = ratio_breaks) +
  labs(x = "Missing rate (%)", y = "Power ratio (Syn-PALM / Obs-GWAS)", subtitle = sub_txt) +
  theme_rev1(10, legend = "right")
ggsave(plot = power_ratio_plot_SynSurrG_vs_ObsG, filename = paste0(FIG, "power_ratio_Syn-PALM_vs_obsG_true_geno_rev1.pdf"),
       device = cairo_pdf, height = 3.4, width = 4.3, units = "in")

data_re <- ratio_data("power_SynSurr")
power_ratio_plot_SynSurrG_vs_SynSurr <- data_re %>%
  mutate(miss_pct = factor(miss * 100, levels = sort(unique(miss * 100)))) %>%
  ggplot(aes(miss_pct, ratio, fill = factor(rho))) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  scale_fill_manual(values = rho_cols, name = expression(rho)) +
  scale_y_continuous(breaks = ratio_breaks, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Missing rate (%)", y = "Power ratio (Syn-PALM / SynSurr)", subtitle = sub_txt) +
  theme_rev1(10, legend = "right") + theme(panel.grid.major.x = element_blank())
ggsave(plot = power_ratio_plot_SynSurrG_vs_SynSurr, filename = paste0(FIG, "power_ratio_Syn-PALM_vs_SynSurr_true_geno_rev1.pdf"),
       device = cairo_pdf, height = 3.4, width = 4.3, units = "in")

cat("\nPower at alpha = 5e-8, pve_g = 0.05%:\n")
print(as.data.frame(power_tab %>% filter(abs(pve_g - 5e-4) < 1e-12) %>%
                      select(miss, rho, n_snps, power_SynSurrG, power_ObsG, power_SynSurr)))

if (length(missing_report)) {
  cat("\n*** INCOMPLETE results; outputs use what exists. Missing:\n")
  cat(paste0("  ", missing_report), sep = "\n")
} else cat("\nAll type I error and power files present.\n")
cat("Written to", FIG, "\n")
