## Source Data for the Nature Communications revision (R1) of the Syn-PALM paper:
## the data underlying the main display items (Figures 2-6, Tables 1-2), one sheet per
## figure panel / table, with a README sheet first.
##
##   module load R/4.4.1-fasrc01
##   Rscript code/post_cf3/source_data/build_source_data.R
##
## Writes  results/source_data_rev1/Source_Data/<sheet>.csv   (one CSV per sheet, README.csv)
##         results/source_data_rev1/Source_Data.xlsx           (via csv_to_xlsx.py, openpyxl)
##         results/source_data_rev1/Source_Data_csv.zip        (the CSVs, as a fallback)
##
## Only summary-level data are written: simulation summaries, per-variant summary statistics,
## per-protein / per-locus / per-group counts and published annotations. No individual-level
## UK Biobank data (no participant identifiers, phenotypes or genotypes) enter any sheet.
## The data of each panel are recomputed from the same inputs and with the same rules as the
## figure scripts named in the README (simulation_code/rev1/*.R, code/post_cf3/figs/*.R).
suppressPackageStartupMessages({ library(data.table) })
BASE <- "<PROJECT_DIR>"
setwd(BASE)
OUT  <- file.path(BASE, "results/source_data_rev1")
CSV  <- file.path(OUT, "Source_Data")
dir.create(CSV, recursive = TRUE, showWarnings = FALSE)
SIM  <- "simulation_code/results_rev1/"
TAB  <- "results/table_results_cf_noPRS_noPC/"
PORTAL <- "https://syn-palm.genohub.org/"

readme <- list()
put <- function(sheet, dt, item, content, cols) {
  dt <- as.data.table(dt)
  stopifnot(nchar(sheet) <= 31, all(names(cols) %in% names(dt)), all(names(dt) %in% names(cols)))
  fwrite(dt, file.path(CSV, paste0(sheet, ".csv")))
  readme[[sheet]] <<- data.table(Sheet = sheet, `Display item` = item, Content = content, Rows = nrow(dt),
                                 Columns = paste(sprintf("%s: %s", names(cols), unname(cols)), collapse = "\n"))
  cat(sprintf("%-16s %8d rows  %s\n", sheet, nrow(dt), item))
}
fmt_alpha <- function(a) formatC(a, format = "e", digits = 0)

######################################################################################
## Figure 2A: robustness (simulation_code/rev1/plot_for_point_robustness_rev1.R)
######################################################################################
q_probs <- c(0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975, 0.99, 0.995)
f2a <- rbindlist(lapply(c(0, 0.75), function(rho) {
  est <- rbindlist(lapply(c(0, 0.25, 0.5, 0.75, 0.9), function(miss) {
    aa <- as.data.table(readRDS(sprintf("%spoint_estimate/GRMr_0.5_missing_%s_rho_%s_rhogeno_siblings_her_0_negative_log10_pvals_for_power.rds",
                                        SIM, miss, rho)))
    rbind(data.table(group = "Obs-GWAS", miss = miss, est = aa$ObsG_hat_beta),
          data.table(group = "Syn-PALM", miss = miss, est = aa$SynSurrG_hat_beta))
  }))
  ## as in PlotSim(): the Obs-GWAS violin pools its estimates of all five settings
  est[group == "Obs-GWAS", miss := NA]
  est[, {
    q <- quantile(est, q_probs, names = FALSE)   # q[6] = first quartile, q[8] = third quartile
    lo <- min(est[est >= q[6] - 1.5 * (q[8] - q[6])]); hi <- max(est[est <= q[8] + 1.5 * (q[8] - q[6])])
    c(list(n_replicates = .N, mean = mean(est), sd = sd(est), min = min(est)),
      setNames(as.list(q), sprintf("q%s", sub("\\.?0+$", "", formatC(100 * q_probs, format = "f", digits = 1)))),
      list(max = max(est), box_lower_whisker = lo, box_upper_whisker = hi))
  }, by = .(group, miss)][, rho := rho][]
}))
f2a[, `:=`(panel = sprintf("rho = %g", rho),
           missing_rate = fifelse(is.na(miss), "all (0, 25, 50, 75, 90%; pooled)", sprintf("%g%%", 100 * miss)))]
f2a <- f2a[, c("panel", "rho", "group", "missing_rate", setdiff(names(f2a), c("panel", "rho", "group", "missing_rate", "miss"))), with = FALSE]
qcols <- grep("^q", names(f2a), value = TRUE)
put("Fig2a", f2a, "Figure 2A",
    "Robustness simulation (sib pairs, GRM r = 0.5, n_obs = 30,000, no genetic effect, so the estimate equals its bias): summary statistics of the distribution of the estimated genetic effect over the simulation replicates for each violin/box plot. The Obs-GWAS violin of each panel pools its estimates of the five missing-rate settings, as in the figure. Box plots: median, quartiles and whiskers (most extreme value within 1.5 x IQR of the quartiles); outliers were not drawn.",
    c(panel = "panel of the figure (left rho = 0, right rho = 0.75)", rho = "correlation of target and synthetic phenotype",
      group = "method (violin)", missing_rate = "missing rate of the target phenotype (Syn-PALM violins)",
      n_replicates = "number of simulation replicates", mean = "mean estimate (= mean bias)", sd = "standard deviation of the estimates",
      min = "minimum", setNames(sprintf("%s percentile of the estimates", sub("^q", "", qcols)), qcols),
      max = "maximum", box_lower_whisker = "lower whisker of the box plot", box_upper_whisker = "upper whisker of the box plot"))

######################################################################################
## Figure 2B: QQ plots, sib pairs (simulation_code/rev1/plot_fig2b_rev1.R)
######################################################################################
acc <- readRDS(paste0(SIM, "figs/t1e_accumulated_rev1.rds"))
BRK <- seq(0, 4, by = 0.002)
qq_points <- function(a) {
  N <- a$n; tl <- sort(a$tail, decreasing = TRUE); nt <- length(tl)
  cs <- cumsum(rev(a$bins)); mid <- rev((BRK[-1] + BRK[-length(BRK)]) / 2)
  rk <- unique(round(exp(seq(log(nt + 1), log(N), length.out = 3000))))
  o_b <- mid[findInterval(rk - nt, c(0, cs), left.open = TRUE)]
  rr <- c(seq_len(nt), rk)
  data.table(rank = rr, n_null_snps = N, expected = -log10((rr - 0.5) / N), observed = c(tl, o_b),
             clower = -log10(qbeta(0.025, rr, N - rr + 1)), cupper = -log10(qbeta(0.975, rr, N - rr + 1)))
}
ks <- names(acc)[startsWith(names(acc), "0.5|") & endsWith(names(acc), "|SynSurrG")]
f2b <- rbindlist(lapply(ks, function(kk) { z <- strsplit(kk, "|", fixed = TRUE)[[1]]
  cbind(missing_rate = as.numeric(z[2]), rho = as.numeric(z[3]), qq_points(acc[[kk]])) }))
f2b <- f2b[missing_rate %in% c(0, 0.5, 0.9)][order(missing_rate, rho, rank)]
setnames(f2b, c("expected", "observed", "clower", "cupper"),
         c("expected_neglog10P", "observed_neglog10P", "ci95_upper_neglog10P", "ci95_lower_neglog10P"))   # -log10 of the 2.5% beta quantile is the upper limit
setcolorder(f2b, c("missing_rate", "rho", "rank", "n_null_snps", "expected_neglog10P", "observed_neglog10P", "ci95_lower_neglog10P", "ci95_upper_neglog10P"))
put("Fig2b", f2b, "Figure 2B",
    "QQ plots of Syn-PALM under the null (sib pairs, GRM r = 0.5) at missing rates 0, 50% and 90%: the plotted points. All SNPs with -log10 P > 4 are given individually; the bulk is represented, as in the figure, by 3,000 log-spaced ranks whose observed value is the midpoint of a 0.002-wide bin of -log10 P. The grey 95% band of the figure is the one of the rho = 0 rows.",
    c(missing_rate = "missing rate (panel)", rho = "correlation of target and synthetic phenotype (colour)",
      rank = "rank of the P value (1 = smallest)", n_null_snps = "number of null SNP tests in the setting",
      expected_neglog10P = "expected -log10 P, -log10((rank - 0.5)/n)", observed_neglog10P = "observed -log10 P",
      ci95_lower_neglog10P = "lower limit of the pointwise 95% interval under the null", ci95_upper_neglog10P = "upper limit of the pointwise 95% interval"))

######################################################################################
## Figures 2C and 3: power (simulation_code/rev1/power_analysis_kinsetting_plots_rev1.R)
######################################################################################
alpha_level <- 5e-08
settings <- expand.grid(miss = c(0, 0.25, 0.5, 0.75, 0.9), rho = c(0, 0.25, 0.5, 0.75), pve_g = seq(0.0001, 0.001, 0.0001))
kin_grm <- c(cousin = 0.125, halfsiblings = 0.25, siblings = 0.5)
kin_lab <- c(cousin = "Cousins", halfsiblings = "Half-siblings", siblings = "Siblings")
pw <- rbindlist(lapply(names(kin_grm), function(kin) rbindlist(lapply(seq_len(nrow(settings)), function(rr) {
  s <- settings[rr, ]
  aa <- readRDS(paste0(BASE, "/", SIM, "power_kinsetting_results/GRMr_", kin_grm[[kin]], "_missing_", s$miss, "_rho_", s$rho,
                       "_rhogeno_", kin, "_her_", s$pve_g, "_negative_log10_pvals_for_power.rds"))
  data.table(relatedness = kin_lab[[kin]], GRMr = kin_grm[[kin]], miss = s$miss, rho = s$rho, pve_g = s$pve_g,
             n_replicates = nrow(aa),
             power_SynPALM = mean(10^(-aa[, "SynSurrG"]) < alpha_level), power_ObsGWAS = mean(10^(-aa[, "ObsG"]) < alpha_level),
             power_SynSurr = mean(10^(-aa[, "SynSurr"]) < alpha_level), power_IndGWAS = mean(10^(-aa[, "Obs"]) < alpha_level))
}))))
f2c <- pw[relatedness == "Siblings" & miss %in% c(0, 0.5, 0.9),
          .(missing_rate = miss, rho, per_SNP_heritability_pct = 100 * pve_g, n_replicates,
            power_SynPALM, power_ObsGWAS)][order(missing_rate, rho, per_SNP_heritability_pct)]
put("Fig2c", f2c, "Figure 2C",
    sprintf("Power at alpha = 5e-8 of Syn-PALM and Obs-GWAS (sib pairs, GRM r = 0.5, n_obs = 30,000) against the per-SNP heritability, at missing rates 0, 50%% and 90%% (panels) and rho = 0-0.75 (line types). Power = share of the replicates with P < 5e-8."),
    c(missing_rate = "missing rate (panel)", rho = "correlation of target and synthetic phenotype (line type)",
      per_SNP_heritability_pct = "per-SNP heritability (%; x axis)", n_replicates = "number of simulation replicates",
      power_SynPALM = "power of Syn-PALM", power_ObsGWAS = "power of Obs-GWAS"))

h7 <- pw[abs(pve_g - 0.0007) < 1e-12]
panels3 <- function(rel, rho) trimws(paste(ifelse(rel == "Cousins", "left (cousins);", ""),
                                           ifelse(rel == "Siblings", "middle (siblings);", ""),
                                           ifelse(rho == 0.75, "right (rho = 0.75);", "")))
f3a <- h7[, .(relatedness, GRMr, rho, missing_rate = miss, n_replicates, power_SynPALM, power_ObsGWAS,
              power_ratio = power_SynPALM / power_ObsGWAS, shown_in_panel = sub(";$", "", panels3(relatedness, rho)))][order(relatedness, rho, missing_rate)]
f3b <- h7[, .(relatedness, GRMr, rho, missing_rate = miss, n_replicates, power_SynPALM, power_SynSurr,
              power_ratio = power_SynPALM / power_SynSurr, shown_in_panel = sub(";$", "", panels3(relatedness, rho)))][order(relatedness, rho, missing_rate)]
c3 <- c(relatedness = "relatedness of the pairs", GRMr = "kinship coefficient of the pairs in the GRM",
        rho = "correlation of target and synthetic phenotype (bar colour in the left and middle panels)",
        missing_rate = "missing rate (x axis)", n_replicates = "number of simulation replicates",
        power_SynPALM = "power of Syn-PALM at alpha = 5e-8")
put("Fig3a", f3a, "Figure 3A",
    "Power ratio Syn-PALM / Obs-GWAS at a per-SNP heritability of 0.07% and alpha = 5e-8 (n_obs = 30,000) for all settings; the bars of the three panels are the rows indicated in shown_in_panel (left: cousins by rho; middle: siblings by rho; right: rho = 0.75 by relatedness). A ratio is undefined (NaN/Inf) when the power of the comparator is 0.",
    c(c3, power_ObsGWAS = "power of Obs-GWAS", power_ratio = "power ratio Syn-PALM / Obs-GWAS (bar height)",
      shown_in_panel = "panel(s) of Figure 3A in which the bar is shown (blank: not shown)"))
put("Fig3b", f3b, "Figure 3B",
    "Power ratio Syn-PALM / SynSurr at a per-SNP heritability of 0.07% and alpha = 5e-8 (n_obs = 30,000), layout as Figure 3A.",
    c(c3, power_SynSurr = "power of SynSurr (independent samples)", power_ratio = "power ratio Syn-PALM / SynSurr (bar height)",
      shown_in_panel = "panel(s) of Figure 3B in which the bar is shown (blank: not shown)"))

######################################################################################
## Table 1: type I error ratios (simulation_code/results_rev1/figs/table1_merged_rev1.csv)
######################################################################################
t1 <- fread(paste0(SIM, "figs/table1_merged_rev1.csv"))
t1[, `:=`(design = fifelse(des == "whole", "Pairs labelled as a whole", "Relatives split"),
          method = fifelse(method == "SynSurrG", "Syn-PALM", "POP-GWAS"))]
t1[, shown_in_Table1 := miss %in% c(0, 0.5, 0.75) & (des == "whole" | miss > 0)]
t1o <- t1[, .(design, GRMr, missing_rate = miss, rho, method, n_null_snps = n,
              `count_1e-05`, `ratio_1e-05`, `poisson_p_1e-05` = `p_1e-05`,
              `count_1e-07` = fifelse(des == "split", NA_real_, as.numeric(`count_1e-07`)),
              `ratio_1e-07` = fifelse(des == "split", NA_real_, as.numeric(`ratio_1e-07`)),
              `poisson_p_1e-07` = fifelse(des == "split", NA_real_, as.numeric(`p_1e-07`)), shown_in_Table1)]
put("Table1", t1o, "Table 1",
    "Type I error of Syn-PALM and POP-GWAS for sib-pair data: number of null SNP tests, number with P below alpha, ratio of the empirical type I error to alpha, and the one-sided exact Poisson P value for a ratio above 1 (bold in Table 1 if < 0.01). Rows at missing rates 0.25 and 0.9 are not in Table 1 (Supplementary Table); alpha = 1e-7 is not estimable in the split design (blank).",
    c(design = "labelling design (left / right block of Table 1)", GRMr = "kinship coefficient of the pairs", missing_rate = "missing rate pi",
      rho = "correlation rho", method = "method", n_null_snps = "number of null SNP tests",
      `count_1e-05` = "number of tests with P < 1e-5", `ratio_1e-05` = "empirical type I error / 1e-5 (Table 1 entry)",
      `poisson_p_1e-05` = "one-sided exact Poisson P value of the count at 1e-5",
      `count_1e-07` = "number of tests with P < 1e-7", `ratio_1e-07` = "empirical type I error / 1e-7 (Table 1 entry)",
      `poisson_p_1e-07` = "one-sided exact Poisson P value of the count at 1e-7", shown_in_Table1 = "row shown in Table 1"))

######################################################################################
## Table 2: height ablation (ablation_results/v2geno_main_results_crossfit_clusterfold)
######################################################################################
t2 <- fread("ablation_results/v2geno_main_results_crossfit_clusterfold/table_counts_5e-8.csv")
t2o <- t2[, .(missing_rate_pct = 100 * missing_rate, OracleGWAS_findings = OracleG,
              SynPALM_findings = SynSurrG_findings_num, SynPALM_overlap = SynSurrG_overlapped, SynPALM_overlap_pct = SynSurrG_overlapped_prob,
              SynPALM_additional = SynSurrG_add_num, ObsGWAS_findings = ObsG_findings_num, ObsGWAS_overlap = ObsG_overlapped,
              ObsGWAS_overlap_pct = ObsG_overlapped_prob, ObsGWAS_additional = ObsG_add_num, rho_oof_synthetic_height = rho_oof)]
put("Table2", t2o, "Table 2",
    "Height ablation in the UK Biobank: numbers of genome-wide significant SNPs (P < 5e-8, autosomes) of Oracle-GWAS, Syn-PALM and Obs-GWAS at each masking rate, the overlap with Oracle-GWAS and the additional SNPs.",
    c(missing_rate_pct = "share of the height values masked (%)", OracleGWAS_findings = "significant SNPs of Oracle-GWAS (all height values)",
      SynPALM_findings = "significant SNPs of Syn-PALM", SynPALM_overlap = "of which also significant in Oracle-GWAS",
      SynPALM_overlap_pct = "overlap as % of the Oracle-GWAS SNPs", SynPALM_additional = "significant in Syn-PALM but not in Oracle-GWAS",
      ObsGWAS_findings = "significant SNPs of Obs-GWAS", ObsGWAS_overlap = "of which also significant in Oracle-GWAS",
      ObsGWAS_overlap_pct = "overlap as % of the Oracle-GWAS SNPs", ObsGWAS_additional = "significant in Obs-GWAS but not in Oracle-GWAS",
      rho_oof_synthetic_height = "out-of-fold correlation of the synthetic and observed height (not shown in Table 2)"))

######################################################################################
## Figures 4 and 5A: UKB-PPP primary associations (code/post_cf3/figs/fig4_panels_cf.R,
## fig4c_top20_cf.R)
######################################################################################
gmap <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
gmap[, chr := as.character(chr)]
per <- fread(paste0(TAB, "summary_per_protein.csv"))
cnt <- per[, .(protein_short_name = tolower(gsub("-", "_", protein)), protein, rho_oof, n_labelled = n_obs,
               syn_count = `SynG_loci_5e-08`, obs_count = `ObsG_loci_5e-08`)]
cat(sprintf("primary associations at 5e-8: Syn-PALM %d, Obs-GWAS %d (%d proteins)\n", sum(cnt$syn_count), sum(cnt$obs_count), nrow(cnt)))

## ---- 4C: top-20 proteins by out-of-fold correlation
top20 <- per[order(-rho_oof)][1:20]
f4c <- top20[, .(rank = .I, protein = tolower(protein), rho_oof, primary_assoc_ObsGWAS = `ObsG_loci_5e-08`, primary_assoc_SynPALM = `SynG_loci_5e-08`)]
put("Fig4c", f4c, "Figure 4C",
    "Numbers of primary associations (P < 5e-8, clumped within +-1 Mb, MHC excluded) of Obs-GWAS and Syn-PALM for the 20 proteins with the highest out-of-fold correlation between synthetic and observed protein levels (in parentheses in the figure).",
    c(rank = "rank by out-of-fold correlation (x-axis order)", protein = "protein (UKB-PPP assay)", rho_oof = "out-of-fold correlation of synthetic and observed protein",
      primary_assoc_ObsGWAS = "primary associations of Obs-GWAS (bar)", primary_assoc_SynPALM = "primary associations of Syn-PALM (bar)"))

## ---- 4D: functional categories
cwp <- merge(cnt, gmap[, .(id, panel)], by.x = "protein_short_name", by.y = "id", all.x = TRUE)
cwp[, panel4 := sub("_II", "", panel, fixed = TRUE)]
f4d <- cwp[, .(proteins = .N, primary_assoc_SynPALM = sum(syn_count), primary_assoc_ObsGWAS = sum(obs_count)), by = .(category = panel4)]
f4d[, `:=`(difference = primary_assoc_SynPALM - primary_assoc_ObsGWAS)]
f4d[, difference_per_protein := round(difference / proteins, 2)]
f4d <- f4d[order(-primary_assoc_SynPALM)]
put("Fig4d", f4d, "Figure 4D",
    "Primary associations (P < 5e-8) of Syn-PALM and Obs-GWAS summed over the proteins of each Olink functional category; the label of each bar gives the number of proteins and the additional Syn-PALM associations per protein.",
    c(category = "Olink Explore functional category (panels I and II combined)", proteins = "number of proteins in the category",
      primary_assoc_SynPALM = "sum of the primary associations of Syn-PALM (bar)", primary_assoc_ObsGWAS = "sum of the primary associations of Obs-GWAS (bar)",
      difference = "Syn-PALM minus Obs-GWAS", difference_per_protein = "difference / proteins (in parentheses)"))

## ---- 4E: diseases and traits (protein annotations of Sun et al. 2023, data/protein_disease.csv)
suppressPackageStartupMessages({ library(readr); library(stringr) })
pdis <- suppressMessages(read_csv("data/protein_disease.csv", show_col_types = FALSE))
main_headers <- colnames(pdis); sub_headers <- as.character(pdis[1, ]); current_main <- ""
for (i in seq_along(main_headers)) {
  if (!is.na(main_headers[i]) && !str_detect(main_headers[i], "^\\.\\.\\.")) current_main <- main_headers[i]
  main_headers[i] <- current_main
}
new_colnames <- paste(main_headers, sub_headers, sep = "_"); new_colnames[1:2] <- c("UKBPPP_ProteinID", "Protein_name")
colnames(pdis) <- new_colnames
pdis <- as.data.table(pdis[-1, ])[!is.na(UKBPPP_ProteinID)]
pdis[, protein_short_name := gsub("-", "_", tolower(sub(":.*", "", UKBPPP_ProteinID)))]
disease_used <- unique(main_headers)[-c(1:2, 28)]
lp_cols <- grep("log10\\(p\\)", names(pdis), value = TRUE)
pdis[, (lp_cols) := lapply(.SD, as.numeric), .SDcols = lp_cols]
cnt_l <- copy(cnt)[, protein_short_name := tolower(protein_short_name)]
f4e_all <- rbindlist(lapply(disease_used, function(d) {
  tc <- str_subset(str_subset(names(pdis), fixed(d)), "log10\\(p\\)")
  if (!length(tc)) return(NULL)
  sig <- pdis[get(tc) > 15][order(-get(tc))]
  m <- cnt_l[protein_short_name %in% sig$protein_short_name]
  data.table(disease_trait = d, annotated_proteins = nrow(sig), proteins = nrow(m),
             primary_assoc_SynPALM = sum(m$syn_count), primary_assoc_ObsGWAS = sum(m$obs_count),
             lead_protein = if (nrow(sig)) sig$protein_short_name[1] else NA_character_,
             lead_protein_log10P_Sun2023 = if (nrow(sig)) sig[[tc]][1] else NA_real_)
}))
f4e_all[, `:=`(difference = primary_assoc_SynPALM - primary_assoc_ObsGWAS)]
f4e_all[, difference_per_protein := round(difference / proteins, 2)]
shown <- f4e_all[order(-difference)][1:15]$disease_trait
f4e_all[, shown_in_Fig4E := disease_trait %in% shown]
f4e_all <- f4e_all[order(-shown_in_Fig4E, -primary_assoc_SynPALM)]
put("Fig4e", f4e_all, "Figure 4E",
    "Primary associations (P < 5e-8) of Syn-PALM and Obs-GWAS summed over the proteins annotated to each disease or trait in Sun et al. (Nature 2023; association with the disease/trait at -log10 P > 15); the figure shows the 15 groups with the largest difference (shown_in_Fig4E), with the lead protein (most significant annotation) at the bar end. All 25 groups are listed.",
    c(disease_trait = "disease or trait", annotated_proteins = "proteins annotated in Sun et al. (-log10 P > 15)",
      proteins = "annotated proteins analysed here (in parentheses)", primary_assoc_SynPALM = "sum of the primary associations of Syn-PALM (bar)",
      primary_assoc_ObsGWAS = "sum of the primary associations of Obs-GWAS (bar)", lead_protein = "protein with the most significant annotation (label at the bar end)",
      lead_protein_log10P_Sun2023 = "-log10 P of the lead protein's association with the disease/trait in Sun et al.",
      difference = "Syn-PALM minus Obs-GWAS", difference_per_protein = "difference / proteins (in parentheses)",
      shown_in_Fig4E = "group shown in Figure 4E (15 largest differences)"))

## ---- primary associations of Syn-PALM with plotting coordinates (5A; basis of 4A, 4B)
chr_order <- as.character(1:22)
chr_off <- gmap[chr %in% chr_order, .(chr_len = max(pos)), by = chr][order(as.integer(chr))]
chr_off[, offset := cumsum(as.numeric(shift(chr_len, fill = 0)))]
loci <- fread(paste0(TAB, "loci_SynG_5e-08.csv"))
pdf_ <- loci[, .(protein_short_name = gsub("-", "_", protein), rho_oof, rsid, chr = as.character(chr), POS38,
                 A1, A2, neglog10P_SynPALM = SynSurrG_negative_log10_pvalues, beta_SynPALM = SynSurrG_hat_beta,
                 se_SynPALM = sqrt(SynSurrG_var_hat_beta), neglog10P_ObsGWAS = ObsG_negative_log10_pvalues,
                 beta_ObsGWAS = ObsG_hat_beta, se_ObsGWAS = sqrt(ObsG_var_hat_beta))]
pdf_ <- merge(pdf_, chr_off[, .(chr, offset)], by = "chr")
pdf_[, qtl_pos_cum := as.numeric(POS38) + offset][, offset := NULL]
pdf_ <- merge(pdf_, gmap[, .(id, gene_chr = chr, gene_pos = pos)], by.x = "protein_short_name", by.y = "id", all.x = TRUE)
pdf_ <- merge(pdf_, chr_off[, .(gene_chr = chr, offset_gene = offset)], by = "gene_chr", all.x = TRUE)
pdf_[, gene_pos_cum := as.numeric(gene_pos) + offset_gene][, offset_gene := NULL]
pdf_[, cis_trans := fifelse(!is.na(gene_chr) & chr == gene_chr & abs(as.numeric(POS38) - as.numeric(gene_pos)) < 1e6, "Cis", "Trans")]
stopifnot(nrow(pdf_) == nrow(loci))
cat(sprintf("Fig 5A: %d primary associations, cis %d, trans %d\n", nrow(pdf_), pdf_[cis_trans == "Cis", .N], pdf_[cis_trans == "Trans", .N]))

source("code/post_cf3/nearest_coding_gene.R")
ng <- unique(pdf_[, .(rsid, chr, POS38)])
ng <- as.data.table(nearest_coding_gene(ng$rsid, ng$chr, ng$POS38))
rsid2gene <- setNames(ng$nearest_gene, ng$rsid)
top_gene <- function(r) { g <- rsid2gene[names(sort(table(r), decreasing = TRUE))[1]]
  if (is.na(g)) names(sort(table(r), decreasing = TRUE))[1] else unname(g) }
top_rsid <- function(r) names(sort(table(r), decreasing = TRUE))[1]

## ---- 4A: 5-Mb windows (circular hotspot plot) and 5A top panel
window_size <- 5e6
genome_max <- max(pdf_$qtl_pos_cum)
pdf_[, window := floor(qtl_pos_cum / window_size)]
setorder(pdf_, protein_short_name, chr, POS38)
win <- pdf_[, .(primary_assoc = .N, distinct_proteins = uniqueN(protein_short_name), pos_cum = mean(qtl_pos_cum),
                chromosomes = paste(unique(chr[order(as.integer(chr))]), collapse = ","),
                first_variant = paste0(chr[which.min(qtl_pos_cum)], ":", POS38[which.min(qtl_pos_cum)]),
                last_variant = paste0(chr[which.max(qtl_pos_cum)], ":", POS38[which.max(qtl_pos_cum)]),
                most_frequent_lead_variant = top_rsid(rsid), label_gene = top_gene(rsid)), by = window][order(window)]
win[, `:=`(bar_height_plotted = pmin(primary_assoc, 250), window_start_cum = window * window_size, window_end_cum = (window + 1) * window_size)]
## Fig 4A labels: windows with > 250 associations, top 11; windows closer than 2.5% of the genome share a label
hl <- win[primary_assoc > 250][order(-primary_assoc)][1:min(11, sum(win$primary_assoc > 250))][order(pos_cum)]
hl[, grp := cumsum(c(1, diff(pos_cum) > 0.025 * genome_max))]
lab4a <- hl[order(grp, -primary_assoc), .(label = paste0(label_gene, " (", primary_assoc, ")", collapse = "; "), windows = list(window)), by = grp]
win[, label_Fig4A := NA_character_]
for (i in seq_len(nrow(lab4a))) win[window %in% lab4a$windows[[i]], label_Fig4A := lab4a$label[i]]
## Fig 5A top panel (Summary_of_pQTLs_ThreePanel): the 11 windows with the most associations
top11 <- win[order(-primary_assoc)][1:11]$window
win[, label_Fig5A_top := fifelse(window %in% top11, fifelse(primary_assoc > 250, paste0(label_gene, " (", primary_assoc, ")"), label_gene), NA_character_)]
f4a <- win[, .(window, window_start_cum, window_end_cum, pos_cum, chromosomes, first_variant, last_variant, primary_assoc,
               bar_height_plotted, distinct_proteins, most_frequent_lead_variant, label_gene, label_Fig4A, label_Fig5A_top)]
put("Fig4a", f4a, "Figure 4A (and top panel of Figure 5A)",
    "Syn-PALM primary associations (P < 5e-8; 24,053) counted in 5-Mb windows of the concatenated autosomes (bar height, capped at 250 in the figure). Hotspots with more than 250 associations (the 11 largest) are labelled by the nearest protein-coding gene (UCSC knownGene hg38) of the most frequent lead variant of the window, with the count. Each protein contributes at most a few lead variants to a window, so the count is close to the number of associated proteins (distinct_proteins). Cumulative coordinate = GRCh38 position + chromosome offset (Fig5a_chr_offsets).",
    c(window = "window index, floor(cumulative position / 5 Mb)", window_start_cum = "window start (cumulative bp)", window_end_cum = "window end (cumulative bp)",
      pos_cum = "mean cumulative position of the variants in the window (bar position)", chromosomes = "chromosome(s) of the variants in the window",
      first_variant = "chr:position (GRCh38) of the first variant in the window", last_variant = "chr:position (GRCh38) of the last variant in the window",
      primary_assoc = "number of Syn-PALM primary associations in the window", bar_height_plotted = "bar height in the figure, min(primary_assoc, 250)",
      distinct_proteins = "number of distinct proteins among them", most_frequent_lead_variant = "most frequent lead variant (rsid) in the window",
      label_gene = "nearest protein-coding gene of the most frequent lead variant", label_Fig4A = "label in Figure 4A (windows within 2.5% of the genome share a label)",
      label_Fig5A_top = "label in the top panel of Figure 5A (11 largest windows)"))

## ---- 4B: primary associations per protein at the gene position (and 5A right panel)
pc <- pdf_[, .(primary_assoc = .N, gene_chr = gene_chr[1], gene_pos = gene_pos[1], gene_pos_cum = mean(gene_pos_cum, na.rm = TRUE),
               top_rsid = rsid[which.max(neglog10P_SynPALM)], rho_oof = rho_oof[1]), by = protein_short_name]
pc[is.nan(gene_pos_cum), gene_pos_cum := NA]
pc[, bar_height_plotted := pmin(primary_assoc, 50)]
top20p <- pc[!is.na(gene_pos_cum)][order(-primary_assoc)][1:20]$protein_short_name
pc[, `:=`(labelled_Fig4B = protein_short_name %in% top20p, plotted = !is.na(gene_pos_cum))]
pc[, protein := toupper(protein_short_name)]
f4b <- pc[order(-primary_assoc, protein), .(protein, rho_oof, gene_chr, gene_pos, gene_pos_cum, primary_assoc, bar_height_plotted,
                                             top_rsid, labelled_Fig4B, plotted)]
put("Fig4b", f4b, "Figure 4B (and right panel of Figure 5A)",
    "Number of Syn-PALM primary associations (P < 5e-8) of each protein with at least one, plotted at the position of its encoding gene (bar height capped at 50); the 20 proteins with the most primary associations are labelled with their numbers. Proteins whose gene is not on an autosome are not plotted.",
    c(protein = "protein", rho_oof = "out-of-fold correlation of synthetic and observed protein", gene_chr = "chromosome of the encoding gene",
      gene_pos = "gene position (GRCh38)", gene_pos_cum = "cumulative gene position (bar position; offsets in Fig5a_chr_offsets)",
      primary_assoc = "number of Syn-PALM primary associations", bar_height_plotted = "bar height in the figure, min(primary_assoc, 50)",
      top_rsid = "lead variant with the smallest P value", labelled_Fig4B = "labelled in Figure 4B (20 proteins with the most primary associations)",
      plotted = "plotted (gene on an autosome)"))

## ---- 5A main panel: cis/trans map
f5a <- pdf_[order(as.integer(chr), POS38, protein_short_name),
            .(protein = toupper(protein_short_name), rsid, chr, POS38, A1, A2, gene_chr, gene_pos, cis_trans,
              qtl_pos_cum, gene_pos_cum, neglog10P_SynPALM, beta_SynPALM, se_SynPALM, neglog10P_ObsGWAS, beta_ObsGWAS, se_ObsGWAS)]
put("Fig5a", f5a, "Figure 5A (main panel)",
    "All 24,053 Syn-PALM primary associations (P < 5e-8, clumped within +-1 Mb, MHC excluded): lead variant position (x) against the position of the gene encoding the protein (y), cis (within 1 Mb of the gene, red) or trans (green), point size -log10 P of Syn-PALM. Summary statistics of Syn-PALM and Obs-GWAS at the lead variant are included. The top and right panels are those of Figure 4A (sheet Fig4a) and Figure 4B (sheet Fig4b).",
    c(protein = "protein", rsid = "lead variant", chr = "chromosome of the variant", POS38 = "position of the variant (GRCh38)",
      A1 = "effect allele", A2 = "other allele", gene_chr = "chromosome of the encoding gene", gene_pos = "position of the encoding gene (GRCh38)",
      cis_trans = "Cis (same chromosome, < 1 Mb from the gene) or Trans (colour)", qtl_pos_cum = "cumulative variant position (x)",
      gene_pos_cum = "cumulative gene position (y)", neglog10P_SynPALM = "-log10 P of Syn-PALM (point size)",
      beta_SynPALM = "effect estimate of Syn-PALM (per effect allele, SD units)", se_SynPALM = "standard error of Syn-PALM",
      neglog10P_ObsGWAS = "-log10 P of Obs-GWAS", beta_ObsGWAS = "effect estimate of Obs-GWAS", se_ObsGWAS = "standard error of Obs-GWAS"))
put("Fig5a_chr_offsets", chr_off[, .(chr, chr_length_used = chr_len, offset_bp = offset)], "Figures 4A, 4B, 5A (coordinates)",
    "Offsets of the concatenated autosomes used for the cumulative coordinates of Figures 4A, 4B and 5A (cumulative = GRCh38 position + offset; chromosome lengths as the largest gene position in the protein map).",
    c(chr = "chromosome", chr_length_used = "length used for the chromosome (bp)", offset_bp = "offset added to the GRCh38 position"))

######################################################################################
## Figure 4F: SNPsea, Syn-PALM (code/post_cf3/figs/fig4f_snpsea.R)
######################################################################################
gct <- fread(cmd = paste("zcat", path.expand("~/tools/snpsea_data/GeneAtlas2004.gct.gz")), skip = 2)
E <- as.matrix(gct[, -(1:2)])
sea <- fread(paste0(TAB, "snpsea/SynG/condition_pvalues.txt"))
bonf <- -log10(0.05 / ncol(E))
x <- sea[order(pvalue)][1:30]
Cm <- cor(E[, x$condition]); o <- hclust(as.dist(1 - Cm), method = "average")$order; lev <- x$condition[o]
sea[, `:=`(neglog10P = -log10(pvalue), bonferroni_significant = -log10(pvalue) > bonf, rank_by_P = frank(pvalue, ties.method = "first"))]
sea[, shown_in_Fig4F := condition %in% x$condition]
sea[, order_in_Fig4F := match(condition, rev(lev))]   # 1 = top row of the panel (ggplot draws the first factor level at the bottom)
f4f <- sea[order(order_in_Fig4F, rank_by_P), .(tissue = condition, pvalue, neglog10P, nulls_observed, nulls_tested, score,
                                               bonferroni_significant, rank_by_P, shown_in_Fig4F, order_in_Fig4F)]
put("Fig4f", f4f, sprintf("Figure 4F (bars)"),
    sprintf("SNPsea tissue enrichment of the genes implicated by the Syn-PALM lead variants (Gene Atlas 2004, %d tissues/cell types): empirical P value per tissue; the figure shows the 30 tissues with the smallest P (rows ordered by UPGMA clustering of their expression profiles), with the Bonferroni threshold 0.05/%d (-log10 = %.2f). P = 1e-6 is the resolution limit of the permutation test.", ncol(E), ncol(E), bonf),
    c(tissue = "Gene Atlas tissue / cell type", pvalue = "SNPsea empirical P value", neglog10P = "-log10 P (bar length)",
      nulls_observed = "null SNP sets with a score at least as large", nulls_tested = "null SNP sets tested", score = "SNPsea specificity score",
      bonferroni_significant = "P < 0.05 / 79 (blue bar, asterisk)", rank_by_P = "rank of the P value", shown_in_Fig4F = "among the 30 tissues shown",
      order_in_Fig4F = "row of the tissue in Figure 4F from the top (blank: not shown)"))
hm <- as.data.table(as.table(Cm)); setnames(hm, c("tissue_1", "tissue_2", "pearson_r"))
hm[, `:=`(i1 = match(tissue_1, lev), i2 = match(tissue_2, lev))]
hm <- hm[i1 >= i2][order(i1, i2)][, .(tissue_1, tissue_2, row_index = i1, column_index = i2, pearson_r)]
put("Fig4f_heatmap", hm, "Figure 4F (heatmap)",
    "Pearson correlations of the Gene Atlas 2004 expression profiles of the 30 tissues shown in Figure 4F (lower triangle drawn as the heatmap; tissues in the UPGMA order of the figure, index 1 = bottom row of the panel).",
    c(tissue_1 = "tissue (row)", tissue_2 = "tissue (column)", row_index = "position of tissue_1 in the clustering order",
      column_index = "position of tissue_2 in the clustering order", pearson_r = "Pearson correlation of the expression profiles (colour)"))

######################################################################################
## Figure 5B: circos of the SuSiE fine-mapping (code/post_cf3/figs/fig5b_circos_susie.R)
######################################################################################
rho_all <- as.data.table(readRDS("results/rho_oof_cf_noPRS_noPC.rds"))[order(-rho_oof)]
sel <- rho_all$protein[1:25]
gmap[, key := gsub("-", "_", id)]
S <- lapply(sel, function(p) readRDS(sprintf("%ssusie/%s.rds", TAB, p))$SynG); names(S) <- sel
var <- rbindlist(lapply(S, `[[`, "var"), fill = TRUE)
cs  <- rbindlist(lapply(S, `[[`, "cs"), fill = TRUE)
var[, `:=`(nlp = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10), in_cs = !is.na(cs))]
f5b_t <- var[order(protein, as.integer(chr), POS38), .(protein = toupper(gsub("_", "/", protein)), lead_rsid, chr, rsid, POS38, z,
                                                       neglog10P_SynPALM = nlp, pip, credible_set = cs, in_95pct_credible_set = in_cs)]
put("Fig5b_track", f5b_t, "Figure 5B (inner track)",
    "Variants with posterior inclusion probability (PIP) >= 0.01 in the SuSiE fine-mapping of the Syn-PALM primary-association loci of the 25 proteins with the highest out-of-fold correlation: the inner track draws a line of height -log10 P (truncated at the top of the track, about 50) at each variant, coloured by protein if the variant is in a 95% credible set (otherwise pink/blue by the sign of z).",
    c(protein = "protein", lead_rsid = "lead variant of the fine-mapped locus", chr = "chromosome", rsid = "variant", POS38 = "position (GRCh38)",
      z = "Syn-PALM z statistic", neglog10P_SynPALM = "-log10 P of Syn-PALM (line height)", pip = "SuSiE posterior inclusion probability",
      credible_set = "95% credible set of the locus the variant belongs to (blank: none)", in_95pct_credible_set = "variant in a 95% credible set"))
gp <- gmap[key %in% sel][!duplicated(key)][as.character(chr) %in% as.character(1:22)]
top <- cs[, .(protein, chr, pos = top_pos, rsid = top_rsid, pip = top_pip, lead_rsid, cs, cs_size, purity_min_abs_corr, variants)]
## as in the figure: best = highest-PIP top variant per protein and chromosome (first in PIP order)
best <- top[order(-pip)][!duplicated(paste(protein, chr))]
lab <- best[pip >= 0.95][order(paste0("chr", chr), pos)][!duplicated(rsid)]
links <- merge(best[pip >= 0.95], gp[, .(protein = key, gene_chr = as.character(chr), gene_pos = pos)], by = "protein")
links[, is_cis := gene_chr == as.character(chr) & abs(gene_pos - pos) < 1e6]
top[, `:=`(labelled = paste(protein, chr, rsid) %in% paste(lab$protein, lab$chr, lab$rsid),
           trans_link = paste(protein, chr, rsid) %in% links[is_cis == FALSE, paste(protein, chr, rsid)],
           cis_next_to_gene = paste(protein, chr, rsid) %in% links[is_cis == TRUE, paste(protein, chr, rsid)])]
top <- merge(top, gp[, .(protein = key, gene_chr = as.character(chr), gene_pos = pos)], by = "protein", all.x = TRUE)
f5b_c <- top[order(protein, as.integer(chr), pos), .(protein = toupper(gsub("_", "/", protein)), gene_chr, gene_pos, lead_rsid, credible_set = cs,
                                                     cs_size, purity_min_abs_corr, chr, top_rsid = rsid, top_pos = pos, top_pip = pip,
                                                     credible_set_variants = variants, labelled_rsid = labelled, trans_link = trans_link,
                                                     cis_next_to_gene = cis_next_to_gene)]
cat(sprintf("Fig 5B: %d credible sets, %d labelled rsids, %d trans links (figure script csv: %d)\n", nrow(f5b_c), sum(f5b_c$labelled_rsid),
            sum(f5b_c$trans_link), nrow(fread(paste0(TAB, "fig5b_circos_susie_links.csv")))))
put("Fig5b_credible_sets", f5b_c, "Figure 5B (labels and links)",
    "The 95% credible sets of the 25 proteins of Figure 5B with their top variant. For each protein and chromosome the top variant with the highest PIP is labelled with its rsid if PIP >= 0.95 (labelled_rsid) and, if it lies on another chromosome or more than 1 Mb from the gene, linked by an arc from the gene of the protein (trans_link); cis signals sit next to the protein label.",
    c(protein = "protein", gene_chr = "chromosome of the encoding gene (protein label)", gene_pos = "gene position (GRCh38)",
      lead_rsid = "lead variant of the fine-mapped locus", credible_set = "credible set (SuSiE effect)", cs_size = "number of variants in the credible set",
      purity_min_abs_corr = "purity (minimum absolute correlation) of the credible set", chr = "chromosome of the credible set",
      top_rsid = "variant with the highest PIP in the credible set", top_pos = "its position (GRCh38)", top_pip = "its PIP",
      credible_set_variants = "variants of the credible set", labelled_rsid = "rsid labelled in the outer ring", trans_link = "arc drawn from the gene to the variant",
      cis_next_to_gene = "PIP >= 0.95 cis signal (no arc; next to the gene label)"))

######################################################################################
## Figure 5C: GDF15 Manhattan plots (code/post_cf3/figs/manhattan_protein_cf.R)
######################################################################################
suppressPackageStartupMessages(source("code/gwas_plot_function.R"))
pv <- as.data.table(readRDS("results/web_SynGWAS_summary_results_cf_noPRS_noPC/gdf15_negative_log10_pvals_v2geno_crossfit.rds"))
L  <- readRDS(paste0(TAB, "gdf15_loci.rds"))
lead_lab <- rbindlist(lapply(c(ObsG = "Obs-GWAS", SynG = "Syn-PALM"), function(m) {
  ld <- L$by_threshold[["5e-08"]][[if (m == "Obs-GWAS") "ObsG" else "SynG"]]$loci
  g <- get_nearest_gene_from_pos(data.frame(rsid = ld$rsid, chr = ld$chr, pos = ld$POS38), build = "hg38")
  data.table(rsid = g$rsid, gene = g$nearest_gene, method = m)
}))
print(lead_lab)
f5c <- pv[ObsG_negative_log10_pvalues >= 3 | SynSurrG_negative_log10_pvalues >= 3 | rsid %in% lead_lab$rsid,
          .(rsid, chr, POS38, REF, ALT, neglog10P_SynPALM = SynSurrG_negative_log10_pvalues, beta_SynPALM = SynSurrG_hat_beta,
            se_SynPALM = sqrt(SynSurrG_var_hat_beta), neglog10P_ObsGWAS = ObsG_negative_log10_pvalues, beta_ObsGWAS = ObsG_hat_beta,
            se_ObsGWAS = sqrt(ObsG_var_hat_beta))][order(chr, POS38)]
f5c[, `:=`(lead_SynPALM_label = lead_lab[method == "Syn-PALM"]$gene[match(rsid, lead_lab[method == "Syn-PALM"]$rsid)],
           lead_ObsGWAS_label = lead_lab[method == "Obs-GWAS"]$gene[match(rsid, lead_lab[method == "Obs-GWAS"]$rsid)])]
cat(sprintf("Fig 5C: %d of %d variants kept\n", nrow(f5c), nrow(pv)))
put("Fig5c", f5c, "Figure 5C",
    sprintf("Manhattan plots of GDF15 (Syn-PALM, left; Obs-GWAS, right): summary statistics of the %s variants with P < 1e-3 in either analysis (the lead variants of the genome-wide significant loci, P < 5e-8, are labelled by their nearest gene). The plots show all %s directly genotyped variants analysed; the full summary statistics of GDF15 (and of all 2,916 proteins) are available on the Syn-PALM web portal %s.",
            format(nrow(f5c), big.mark = ","), format(nrow(pv), big.mark = ","), PORTAL),
    c(rsid = "variant", chr = "chromosome", POS38 = "position (GRCh38)", REF = "reference allele", ALT = "alternative (effect) allele",
      neglog10P_SynPALM = "-log10 P of Syn-PALM (y, left plot)", beta_SynPALM = "effect estimate of Syn-PALM", se_SynPALM = "standard error of Syn-PALM",
      neglog10P_ObsGWAS = "-log10 P of Obs-GWAS (y, right plot)", beta_ObsGWAS = "effect estimate of Obs-GWAS", se_ObsGWAS = "standard error of Obs-GWAS",
      lead_SynPALM_label = "gene label of a Syn-PALM lead variant (red, left plot)", lead_ObsGWAS_label = "gene label of an Obs-GWAS lead variant (red, right plot)"))

######################################################################################
## Figure 6: validation (code/post_cf3/figs/fig_validation_cf.R)
######################################################################################
xv <- fread(paste0(TAB, "external_validation_loci.csv"))
xv[, grp := fifelse(method == "ObsG", "Obs-GWAS", fifelse(only, "Syn-PALM, additional", "Syn-PALM, shared"))]
pa <- rbind(xv[catalog_records > 0, .(replicated = sum(in_catalog), total = .N), by = grp][, criterion := "Reported in the GWAS Catalog (proteins with Catalog records)"],
            xv[!is.na(sun_sig), .(replicated = sum(sun_sig), total = .N), by = grp][, criterion := "Lead variant P < 5e-8 in Sun et al. (N = 51,507)"])
pa[, share_pct := 100 * replicated / total]
f6a <- pa[order(criterion != "Reported in the GWAS Catalog (proteins with Catalog records)", match(grp, c("Obs-GWAS", "Syn-PALM, shared", "Syn-PALM, additional"))),
          .(criterion, group = grp, replicated, total, share_pct)]
put("Fig6a", f6a, "Figure 6A",
    "Share of the primary associations (P < 5e-8) reported in the GWAS Catalog for the same protein within +-1 Mb (top; proteins with Catalog records) and whose lead variant has P < 5e-8 in the combined UKB-PPP analysis of Sun et al. (bottom), pooled over proteins, for Obs-GWAS and for the Syn-PALM associations shared with Obs-GWAS or additional (no Obs-GWAS primary association within +-1 Mb).",
    c(criterion = "panel (top / bottom)", group = "group of primary associations (bar)", replicated = "associations meeting the criterion",
      total = "associations in the group (with data)", share_pct = "share (%; bar length)"))
so <- xv[method == "SynG" & only == TRUE][order(protein, chr, POS38)]
so[, sun := fifelse(is.na(sun_sig), "no data", fifelse(sun_sig, "yes", "no"))]
cat(sprintf("Fig 6B-D: %d additional associations; Obs P<1e-3 %.1f%%, <1e-5 %.1f%%; median SE ratio %.3f\n", nrow(so),
            100 * mean(so$ObsG_nlp >= 3), 100 * mean(so$ObsG_nlp >= 5), median(so$SynSurrG_se / so$ObsG_se)))
idc <- c(protein = "protein", rsid = "Syn-PALM lead variant", chr = "chromosome", POS38 = "position (GRCh38)")
put("Fig6b", so[, .(protein, rsid, chr, POS38, neglog10P_SynPALM = SynSurrG_nlp, neglog10P_ObsGWAS = ObsG_nlp, P_lt_5e8_in_Sun2023 = sun)],
    "Figure 6B",
    "The 4,538 additional primary associations of Syn-PALM (no Obs-GWAS primary association within +-1 Mb), at the Syn-PALM lead variant: -log10 P of Obs-GWAS against Syn-PALM, coloured by P < 5e-8 in Sun et al. (Syn-PALM axis on a log scale, limited to 30 in the figure).",
    c(idc, neglog10P_SynPALM = "-log10 P of Syn-PALM (x)", neglog10P_ObsGWAS = "-log10 P of Obs-GWAS (y)",
      P_lt_5e8_in_Sun2023 = "lead variant with P < 5e-8 in Sun et al. (colour; 'no data': variant not available, grey)"))
put("Fig6c", so[, .(protein, rsid, chr, POS38, beta_SynPALM = SynSurrG_beta, beta_ObsGWAS = ObsG_beta)], "Figure 6C",
    "The 4,538 additional primary associations of Syn-PALM: effect estimate of Obs-GWAS against Syn-PALM at the Syn-PALM lead variant (the figure is limited to the 99.8th percentile of the absolute estimates).",
    c(idc, beta_SynPALM = "effect estimate of Syn-PALM (x)", beta_ObsGWAS = "effect estimate of Obs-GWAS (y)"))
put("Fig6d", so[, .(protein, rsid, chr, POS38, se_SynPALM = SynSurrG_se, se_ObsGWAS = ObsG_se, se_ratio = SynSurrG_se / ObsG_se)], "Figure 6D",
    "The 4,538 additional primary associations of Syn-PALM: ratio of the standard errors of Syn-PALM and Obs-GWAS at the Syn-PALM lead variant (histogram, bin width 0.0125; red line, median).",
    c(idc, se_SynPALM = "standard error of Syn-PALM", se_ObsGWAS = "standard error of Obs-GWAS", se_ratio = "SE ratio Syn-PALM / Obs-GWAS (histogram)"))
sp <- fread("results/split_validation_cf_noPRS_noPC_summary.csv")
sp[, Method := fifelse(method == "SynSurrG", "Syn-PALM", "Obs-GWAS")]
pe <- rbind(sp[, .(level = "GWS SNPs", criterion = "P < 0.05, same sign", replicated = sum(snp_nom), discoveries = sum(snp_disc)), by = Method],
            sp[, .(level = "GWS SNPs", criterion = "P < 5e-8", replicated = sum(snp_overlap), discoveries = sum(snp_disc)), by = Method],
            sp[, .(level = "Primary associations", criterion = "P < 0.05, same sign", replicated = sum(locus_nom), discoveries = sum(locus_disc)), by = Method],
            sp[, .(level = "Primary associations", criterion = "P < 5e-8", replicated = sum(locus_overlap), discoveries = sum(locus_disc)), by = Method])
pe[, share_pct := 100 * replicated / discoveries]
put("Fig6e", pe[order(level, criterion, Method != "Syn-PALM"), .(level, criterion, method = Method, replicated, discoveries, share_pct)], "Figure 6E",
    "Split-sample analysis of the 100 proteins with the highest prediction accuracy (two halves of whole relatedness clusters): share of the genome-wide significant SNPs and of the primary associations of the discovery half replicated in the validation half (P < 0.05 with the same sign, or P < 5e-8), summed over proteins.",
    c(level = "genome-wide significant SNPs or primary associations (panel)", criterion = "replication criterion in the validation half (x)",
      method = "method (bar colour)", replicated = "discoveries replicated", discoveries = "discoveries in the discovery half", share_pct = "share replicated (%; bar height)"))
pf <- dcast(sp, protein + rho_disc + n_obs_disc ~ Method, value.var = "snp_disc")
setnames(pf, c("Obs-GWAS", "Syn-PALM", "n_obs_disc", "rho_disc"), c("GWS_SNPs_ObsGWAS", "GWS_SNPs_SynPALM", "n_labelled_discovery_half", "rho_oof_discovery_half"))
pf <- pf[order(-GWS_SNPs_SynPALM), .(protein, rho_oof_discovery_half, n_labelled_discovery_half, GWS_SNPs_ObsGWAS, GWS_SNPs_SynPALM)]
cat(sprintf("Fig 6F: Syn %d, Obs %d; more Syn %d of %d\n", sum(pf$GWS_SNPs_SynPALM), sum(pf$GWS_SNPs_ObsGWAS),
            sum(pf$GWS_SNPs_SynPALM > pf$GWS_SNPs_ObsGWAS), nrow(pf)))
put("Fig6f", pf, "Figure 6F",
    "Split-sample analysis: number of genome-wide significant SNPs (P < 5e-8) per protein in the discovery half, Syn-PALM against Obs-GWAS (plotted as count + 1 on log scales).",
    c(protein = "protein", rho_oof_discovery_half = "out-of-fold correlation of synthetic and observed protein in the discovery half",
      n_labelled_discovery_half = "labelled individuals (with protein measurement) in the discovery half",
      GWS_SNPs_ObsGWAS = "genome-wide significant SNPs of Obs-GWAS (x, plus 1)", GWS_SNPs_SynPALM = "genome-wide significant SNPs of Syn-PALM (y, plus 1)"))

######################################################################################
## README and Excel workbook
######################################################################################
rd <- rbindlist(readme)
sheet_order <- c("Fig2a", "Fig2b", "Fig2c", "Fig3a", "Fig3b", "Table1", "Table2", "Fig4a", "Fig4b", "Fig4c", "Fig4d", "Fig4e",
                 "Fig4f", "Fig4f_heatmap", "Fig5a", "Fig5a_chr_offsets", "Fig5b_track", "Fig5b_credible_sets", "Fig5c",
                 "Fig6a", "Fig6b", "Fig6c", "Fig6d", "Fig6e", "Fig6f")
stopifnot(setequal(sheet_order, rd$Sheet))
rd <- rd[match(sheet_order, Sheet)]
notes <- data.table(Sheet = c("", "", "", "", ""), `Display item` = "", Content = c(
  "Source Data for: Syn-PALM (Nature Communications, revision R1). One sheet per figure panel or table; this README lists each sheet, the display item it underlies, its content and the meaning of its columns.",
  "Figure 1 is a schematic and contains no data.",
  "Only aggregate / summary-level data are provided: simulation summaries, GWAS summary statistics of individual variants, and per-protein, per-locus or per-group counts. No individual-level UK Biobank data (participant identifiers, phenotypes, protein measurements or genotypes) are included.",
  sprintf("The full genome-wide summary statistics of Syn-PALM and Obs-GWAS for all 2,916 proteins (the basis of Figures 4-6) are available on the Syn-PALM web portal %s; sheet Fig5c gives a reduced set (P < 1e-3) of the GDF15 summary statistics plotted in Figure 5C.", PORTAL),
  "Syn-PALM = method of this paper (SynSurrG in the code); Obs-GWAS = linear mixed model GWAS of the observed phenotype; primary association = lead variant of a +-1 Mb clump of genome-wide significant variants (P < 5e-8; MHC chr6:25.5-34.0 Mb excluded); rho_oof = out-of-fold correlation of the synthetic and the observed protein level."),
  Rows = NA_integer_, Columns = "")
fwrite(rbind(notes, rd), file.path(CSV, "README.csv"))
order_file <- file.path(CSV, "sheet_order.txt")
writeLines(c("README", rd$Sheet), order_file)
py <- "python3"
st <- system2(py, c("-P", file.path(BASE, "code/post_cf3/source_data/csv_to_xlsx.py"), CSV, order_file, file.path(OUT, "Source_Data.xlsx")))
if (st != 0) stop("csv_to_xlsx.py failed; the CSV files are in ", CSV)
old <- setwd(OUT); unlink("Source_Data_csv.zip")
utils::zip("Source_Data_csv.zip", file.path("Source_Data", c("README.csv", paste0(rd$Sheet, ".csv"))), flags = "-q9X")
setwd(old)
cat("wrote", file.path(OUT, "Source_Data.xlsx"), "and", file.path(OUT, "Source_Data_csv.zip"), "\n")
