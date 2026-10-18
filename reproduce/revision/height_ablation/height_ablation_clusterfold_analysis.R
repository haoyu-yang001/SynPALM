## Tables and figure for the height ablation rerun with relatedness-cluster
## cross-fit folds (height_ablation_crossfit_clusterfold.R). Same quantities as
## ablation_v2_analysis.R; the Oracle columns are taken from the 0% file because
## the oracle does not depend on the missing rate.
##   Rscript simulation_code/height_ablation_clusterfold_analysis.R
## Writes to ablation_results/v2geno_main_results_crossfit_clusterfold/:
##   table_counts_5e-8.csv   Oracle / Syn-PALM / Obs-GWAS hits, overlap with the
##                           oracle and additional hits (manuscript height table)
##   table_mean_se.csv       mean SE of beta-hat at the oracle hits (Supp)
##   ratio_plot.pdf          ratio of hits, Syn-PALM vs SynSurr / Obs-GWAS / Obs
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(scales) })
setwd("<PROJECT_DIR>")
D     <- "ablation_results/v2geno_main_results_crossfit_clusterfold/"
RATES <- c(0, 0.25, 0.5, 0.75, 0.9)
ALPHA <- 5e-8
f_res <- function(r) paste0(D, "height_missing_", r, "_negative_log10_pvals_v2geno_crossfit_clusterfold.rds")

## autosomes only, as stated in the paper (the .bed also holds X, XY, Y and MT)
bim  <- data.table::fread("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bim", header = FALSE,
                          col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))
auto <- bim$snp[bim$chr %in% as.character(1:22)]
read_res <- function(r) { x <- readRDS(f_res(r)); as.data.frame(x[rownames(x) %in% auto, , drop = FALSE]) }

oracle <- read_res(0)
cat("autosomal SNPs:", nrow(oracle), "\n")
hit_or <- which(10^(-oracle$OracleG_negative_log10_pvalues) < ALPHA)

counts <- list(); se <- list(); rho <- c()
for (r in RATES) {
  x <- read_res(r)
  stopifnot(identical(rownames(x), rownames(oracle)))
  hit <- function(m) which(10^(-x[[paste0(m, "_negative_log10_pvalues")]]) < ALPHA)
  hS <- hit("SynSurrG"); hO <- hit("ObsG")
  counts[[length(counts) + 1]] <- data.frame(
    missing_rate = r, OracleG = length(hit_or),
    SynSurrG_findings_num = length(hS), SynSurrG_overlapped = length(intersect(hS, hit_or)),
    SynSurrG_overlapped_prob = round(100 * length(intersect(hS, hit_or)) / length(hit_or), 2),
    SynSurrG_add_num = length(setdiff(hS, hit_or)),
    ObsG_findings_num = length(hO), ObsG_overlapped = length(intersect(hO, hit_or)),
    ObsG_overlapped_prob = round(100 * length(intersect(hO, hit_or)) / length(hit_or), 2),
    ObsG_add_num = length(setdiff(hO, hit_or)),
    SynSurr_findings_num = length(hit("SynSurr")), Obs_findings_num = length(hit("Obs")))
  msd <- function(m) mean(sqrt(x[[paste0(m, "_var_hat_beta")]][hit_or]))
  se[[length(se) + 1]] <- data.frame(missing_rate = r,
    OracleG = mean(sqrt(oracle$OracleG_var_hat_beta[hit_or])),
    ObsG = msd("ObsG"), SynSurrG = msd("SynSurrG"), Obs = msd("Obs"), SynSurr = msd("SynSurr"))
  fr <- paste0(D, "height_miss_", r, "_correlation_clusterfold.rds")
  rho <- c(rho, if (file.exists(fr)) unlist(readRDS(fr))[1] else NA)
}
counts <- bind_rows(counts); se <- bind_rows(se)
counts$rho_oof <- rho
write.csv(counts, paste0(D, "table_counts_5e-8.csv"), row.names = FALSE)
write.csv(se, paste0(D, "table_mean_se.csv"), row.names = FALSE)
print(counts)
print(cbind(missing_rate = se$missing_rate, apply(se[, -1], 2, formatC, format = "e", digits = 2)), quote = FALSE)

## ratio plot (same layout as plot_figures/ablation_v2_ratio.pdf)
df <- tibble(
  setting = factor(RATES, levels = RATES, labels = percent_format(accuracy = 1)(RATES)),
  `Syn-PALM vs SynSurr`  = counts$SynSurrG_findings_num / counts$SynSurr_findings_num,
  `Syn-PALM vs Obs-GWAS` = counts$SynSurrG_findings_num / counts$ObsG_findings_num,
  `Syn-PALM vs Obs`      = counts$SynSurrG_findings_num / counts$Obs_findings_num) %>%
  pivot_longer(-setting, names_to = "Methods", values_to = "ratio") %>%
  mutate(Methods = factor(Methods, levels = c("Syn-PALM vs SynSurr", "Syn-PALM vs Obs-GWAS", "Syn-PALM vs Obs")))
g <- ggplot(df, aes(x = setting, y = ratio, fill = Methods)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  xlab("Missing Rate (%)") + ylab("Ratio of GWAS SNPs discovered") +
  scale_x_discrete(labels = c("0", "25", "50", "75", "90")) +
  scale_fill_manual(breaks = levels(df$Methods),
                    values = c("Syn-PALM vs Obs-GWAS" = "#d95f02", "Syn-PALM vs SynSurr" = "#1b9e77",
                               "Syn-PALM vs Obs" = "#0073C2FF")) +
  scale_y_continuous(breaks = function(x) { b <- pretty(x); if (!1 %in% b) b <- sort(c(b, 1)); b },
                     limits = c(0, max(df$ratio, na.rm = TRUE) + 0.2), expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8))
ggsave(paste0(D, "ratio_plot.pdf"), g, device = "pdf", height = 3.5, width = 5, units = "in")
cat("wrote", paste0(D, "table_counts_5e-8.csv"), "and ratio_plot.pdf\n")
