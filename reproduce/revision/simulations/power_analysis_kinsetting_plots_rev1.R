## Revision copy of simulation_code/power_analysis_kinsetting_plots.R: reads simulation_code/results_rev1/ (rho = 0.75 rerun,
## other rho unchanged; see merge_power_bias_rev1.R) and writes to simulation_code/results_rev1/figs/.
####### combine results ##########

library(purrr)
library(dplyr)
library(tidyr)

miss_vec <- c(0,0.25,0.5,0.75,0.9)
rho_vec <- c(0,0.25,0.5,0.75)
snp_heritability_vec <- seq(0.0001, 0.001, 0.0001)

settings <- expand.grid(miss = miss_vec, rho = rho_vec, pve_g = snp_heritability_vec)

# 定义一个根据 GRMr rho_geno 读取并拼接所有 setting 的函数
read_power_results <- function(GRMr, kinsetting) {
  map_dfr(
    .x = seq_len(nrow(settings)),
    .f = function(rr) {
      miss  <- settings$miss[rr]
      rho   <- settings$rho[rr]
      pve_g <- settings$pve_g[rr]
      
      path <- paste0("<PROJECT_DIR>/simulation_code/results_rev1/power_kinsetting_results/GRMr_",
                     GRMr,"_missing_",miss,"_rho_",rho,"_rhogeno_",kinsetting,"_her_",pve_g, "_negative_log10_pvals_for_power.rds")
      
      aa <- readRDS(path)
      
      # 保持原来的列结构：setting_num, miss, rho, pve_g, 然后是 aa 中的列
      data.frame(
        setting_num = rr,
        miss        = miss,
        rho         = rho,
        pve_g       = pve_g,
        aa
      )
    }
  )
}

### simulation settings we have :
## 1. correlation between Y and hatY 0 0.25 0.5 0.75
## 2. missing rate 0 0.25 0.5 0.75 0.9
## 3. GRMr 0.125 0.25 0.5
## 4. kinsettings "cousin" 0.125 "halfsiblings" 0.25 "siblings" 0.5

alpha_level <- 5e-08

library(ggplot2)

## ---- common look (theme_rev1.R): method colours / shapes, rho line types and blues,
## relatedness purples; restyled 2026-10-16 (previous version: code/before/
## power_analysis_kinsetting_plots_rev1_2026-10-16_pre_restyle.R)
source("<PROJECT_DIR>/simulation_code/rev1/theme_rev1.R")
FIG <- "<PROJECT_DIR>/simulation_code/results_rev1/figs/"
kin_grm <- c(cousin = 0.125, halfsiblings = 0.25, siblings = 0.5)

## power curves against the per-SNP heritability, one panel per missing rate.
## d: miss, rho, pve_g, method (label), power. by = "method": colour/shape by method and
## line type by rho; by = "rho": a single method, colour by rho.
power_curve_plot <- function(d, by = c("method", "rho"), nrow = 2, base_size = 12) {
  by <- match.arg(by)
  d$rho_f <- factor(d$rho, names(rho_lty))
  d$panel <- factor(miss_lab(d$miss), miss_lab(miss_vec))
  if (by == "method") {
    meth <- levels(d$method)
    p <- ggplot(d, aes(pve_g * 100, power, colour = method, shape = method, linetype = rho_f,
                       group = interaction(rho_f, method))) +
      geom_line(linewidth = 0.55) + geom_point(size = 1.5) +
      scale_colour_manual(values = method_cols[meth], breaks = meth, name = "Method") +
      scale_shape_manual(values = method_shapes[meth], breaks = meth, name = "Method") +
      scale_linetype_manual(values = rho_lty, name = expression(rho), drop = TRUE) +
      guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1),
             linetype = guide_legend(order = 2, override.aes = list(colour = "grey20", shape = NA)))
  } else {
    p <- ggplot(d, aes(pve_g * 100, power, colour = rho_f, group = rho_f)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
      scale_colour_manual(values = rho_cols, name = expression(rho))
  }
  p + facet_wrap(~ panel, nrow = nrow) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
    scale_x_continuous(breaks = if (nrow == 1) scales::breaks_width(0.04, offset = 0.02) else scales::breaks_width(0.02)) +
    labs(x = "Per-SNP heritability (%)", y = power_lab(alpha_level)) +
    theme_rev1(base_size, legend = "right") +
    theme(legend.key.width = unit(2.2, "lines"))
}
## legend in the empty sixth cell of a 2 x 3 grid of five panels
## (stacked; side by side for the four-method figure, whose method legend is taller)
legend_in_empty_cell <- function(box = "vertical")
  theme(legend.position = "inside", legend.position.inside = c(if (box == "vertical") 0.855 else 0.835, if (box == "vertical") 0.23 else 0.2),
        legend.justification = c(0.5, 0.5), legend.box = box, legend.box.just = c("left", "top"),
        legend.spacing.y = unit(6, "pt"), legend.spacing.x = unit(10, "pt"),
        legend.key.width = unit(if (box == "vertical") 2.2 else 1.6, "lines"))

## bar chart of power ratios at one per-SNP heritability (Fig 3 panels)
## d: miss, ratio and the fill column; the y scale (common to the panels of a row) is
## added by save_ratio()
ratio_bar_plot <- function(d, fill, values, fill_name, title, ylab_txt, sub_txt) {
  d$miss_pct <- factor(d$miss * 100, levels = sort(unique(d$miss * 100)))
  ggplot(d, aes(miss_pct, ratio, fill = .data[[fill]])) +
    geom_col(position = position_dodge(width = 0.82), width = 0.75, colour = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.45) +
    scale_fill_manual(values = values, name = fill_name) +
    labs(x = "Missing rate (%)", y = ylab_txt, title = title, subtitle = sub_txt) +
    theme_rev1(12.5, legend = "bottom") +
    theme(panel.grid.major.x = element_blank(), axis.text = element_text(size = rel(0.92)),
          legend.key.width = unit(0.85, "lines"), legend.key.height = unit(0.85, "lines"),
          legend.text = element_text(size = rel(0.8), margin = margin(l = 2.5)),
          legend.key.spacing.x = unit(3, "pt"), legend.margin = margin(0, 0, 0, 0), legend.location = "plot",
          plot.margin = margin(4, 6, 2, 3))
}
save_ratio <- function(p, name, ymax) ggsave(
  plot = p + scale_y_continuous(limits = c(0, ymax), breaks = ratio_breaks_rev1, expand = expansion(mult = c(0, 0.03))),
  filename = paste0(FIG, name), device = cairo_pdf, height = 3, width = 3, units = "in")

##### 1. power plot for SynSurrG  and SynSurrG compared with others ################

for (kinsetting in c("cousin","halfsiblings","siblings")) {
  
  GRMr <- kin_grm[[kinsetting]]
  
  com_results <- read_power_results(GRMr,kinsetting)
  
  com_results$power_SynSurrG <- 10^{-com_results$SynSurrG} < alpha_level
  com_results$power_ObsG <- 10^{-com_results$ObsG} < alpha_level
  com_results$power_OracleG <- 10^{-com_results$OracleG} < alpha_level
  com_results$power_Obs <- 10^{-com_results$Obs} < alpha_level
  com_results$power_Oracle <- 10^{-com_results$Oracle} < alpha_level
  com_results$power_SynSurr <- 10^{-com_results$SynSurr} < alpha_level
  
  ## mean power per setting and method
  meth_lab <- c(power_SynSurrG = "Syn-PALM", power_ObsG = "Obs-GWAS", power_SynSurr = "SynSurr", power_Obs = "Ind-GWAS")
  pw <- com_results %>%
    pivot_longer(cols = c(power_ObsG, power_SynSurrG, power_SynSurr, power_Obs),
                 names_to = "method", values_to = "power") %>%
    group_by(miss, rho, pve_g, method) %>%
    summarise(power = mean(power), .groups = "drop") %>%
    mutate(method = factor(meth_lab[method], meth_lab))
  
  ## (a) Syn-PALM only, colour by rho
  power_plot <- power_curve_plot(droplevels(filter(pw, method == "Syn-PALM")), by = "rho") + legend_in_empty_cell()
  ggsave(plot = power_plot, filename = paste0(FIG, "power_Syn-PALM_kinsetting_", kinsetting, ".pdf"),
         device = cairo_pdf, height = 5, width = 7.2, units = "in")
  
  ## (b) Syn-PALM, Obs-GWAS, SynSurr and Ind-GWAS
  power_plot <- power_curve_plot(pw, by = "method") + legend_in_empty_cell("horizontal")
  ggsave(plot = power_plot, filename = paste0(FIG, "power_Syn-PALM_others_kinsetting_", kinsetting, ".pdf"),
         device = cairo_pdf, height = 5, width = 7.2, units = "in")
  
  ## (c) Syn-PALM and Obs-GWAS (Supp Figs S5-S7)
  pw2 <- droplevels(filter(pw, method %in% c("Syn-PALM", "Obs-GWAS")))
  power_plot <- power_curve_plot(pw2, by = "method") + legend_in_empty_cell()
  ggsave(plot = power_plot, filename = paste0(FIG, "power_Syn-PALM_ObsG_kinsetting_", kinsetting, ".pdf"),
         device = cairo_pdf, height = 5, width = 7.2, units = "in")
  
  ## Fig 2C of the main text: sib pairs, missing rates 0, 50% and 90% in one row
  ## (layout of the submitted figure)
  if (kinsetting == "siblings") {
    p2c <- power_curve_plot(filter(pw2, miss %in% c(0, 0.5, 0.9)), by = "method", nrow = 1)
    ggsave(plot = p2c, filename = paste0(FIG, "fig2c_power_siblings_rev1.pdf"),
           device = cairo_pdf, height = 2.45, width = 7, units = "in")
  }
  
}


####### 2. power ratio compare across different settings #################
## Fig 3: bar charts at one per-SNP heritability; the panels of one comparison share
## the y range (saved after all of them are made)

ind_pve_g <- 7
rho_name <- expression(rho)
ylab_obs <- "Power ratio\n(Syn-PALM / Obs-GWAS)"
ylab_ss  <- "Power ratio\n(Syn-PALM / SynSurr)"

##### compare SynSurrG with ObsG

## given a heritability and every rho
## for each GRMr and kinsetting

plots_obs <- list()
for (kinsetting in c("cousin","halfsiblings","siblings")) {
  
  GRMr <- kin_grm[[kinsetting]]
  kinsetting_label <- kin_lab[[kinsetting]]
  
  com_results <- read_power_results(GRMr,kinsetting)
  
  com_results$power_ObsG <- 10^{-com_results$ObsG} < alpha_level
  com_results$power_SynSurrG <- 10^{-com_results$SynSurrG} < alpha_level
  
  data_re0 <- com_results[, c("setting_num","miss","rho","pve_g","power_SynSurrG","power_ObsG")] %>%
    group_by(miss, rho, pve_g) %>% summarise_at(vars(power_SynSurrG, power_ObsG), mean) %>%
    mutate(ratio = power_SynSurrG / power_ObsG)
  
  data_re <- data_re0[which(data_re0$pve_g %in% unique(data_re0$pve_g)[c(ind_pve_g)]),]
  sub_txt <- sprintf("Per-SNP heritability %g%%", 100 * unique(data_re$pve_g))
  data_re$rho_f <- factor(data_re$rho, names(rho_cols))
  
  plots_obs[[paste0("kinsetting_", kinsetting)]] <-
    ratio_bar_plot(data_re, "rho_f", rho_cols, rho_name, kinsetting_label, ylab_obs, sub_txt)
}

## given a heritability and every GRMr
## for each rho and kinsetting

com_results_grm01 <- read_power_results(GRMr = 0.125,kinsetting = "cousin")
com_results_grm02 <- read_power_results(GRMr = 0.25,kinsetting = "halfsiblings")
com_results_grm03 <- read_power_results(GRMr = 0.5,kinsetting = "siblings")

com_results_grm01$kinsettings <- "cousin"
com_results_grm02$kinsettings <- "halfsiblings"
com_results_grm03$kinsettings <- "siblings"

com_results_all <- rbind(com_results_grm01,com_results_grm02,com_results_grm03)

com_results_all$power_SynSurrG <- 10^{-com_results_all$SynSurrG} < alpha_level
com_results_all$power_ObsG <- 10^{-com_results_all$ObsG} < alpha_level


data_re0 <- com_results_all[, c("setting_num","miss","rho","pve_g","power_SynSurrG","power_ObsG","kinsettings")] %>%
  group_by(miss, rho, pve_g,kinsettings) %>%
  summarise_at(vars(power_SynSurrG, power_ObsG), mean) %>%
  mutate(ratio = power_SynSurrG / power_ObsG)

for (rho in c(0,0.25,0.5,0.75)) {
  data_re <- data_re0[which(data_re0$pve_g %in% unique(data_re0$pve_g)[c(ind_pve_g)] & data_re0$rho %in% c(rho)),]
  sub_txt <- sprintf("Per-SNP heritability %g%%", 100 * unique(data_re$pve_g))
  data_re$kin_f <- factor(kin_lab[data_re$kinsettings], kin_lab)
  
  plots_obs[[paste0("rho_", rho)]] <-
    ratio_bar_plot(data_re, "kin_f", kin_cols, NULL, sprintf("\u03c1 = %g", rho), ylab_obs, sub_txt)
}

## common y range: all ratios at this heritability
sel <- data_re0$pve_g %in% unique(data_re0$pve_g)[ind_pve_g]
ymax_obs <- max(data_re0$ratio[sel][is.finite(data_re0$ratio[sel])], 1)
for (nm in names(plots_obs))
  save_ratio(plots_obs[[nm]], ymax = ymax_obs, paste0("power_ratio_Syn-PALM_vs_obsG_", nm, ".pdf"))

##### compare SynSurrG with SynSurr

## given a heritability and every rho
## for each GRMr and kinsetting

plots_ss <- list()
for (kinsetting in c("cousin","halfsiblings","siblings")) {
  
  GRMr <- kin_grm[[kinsetting]]
  kinsetting_label <- kin_lab[[kinsetting]]
  
  com_results <- read_power_results(GRMr,kinsetting)
  
  com_results$power_SynSurr <- 10^{-com_results$SynSurr} < alpha_level
  com_results$power_SynSurrG <- 10^{-com_results$SynSurrG} < alpha_level
  
  data_re0 <- com_results[, c("setting_num","miss","rho","pve_g","power_SynSurrG","power_SynSurr")] %>%
    group_by(miss, rho, pve_g) %>%
    summarise_at(vars(power_SynSurrG, power_SynSurr), mean) %>%
    mutate(ratio = power_SynSurrG / power_SynSurr)
  
  data_re <- data_re0[which(data_re0$pve_g %in% unique(data_re0$pve_g)[c(ind_pve_g)]),]
  sub_txt <- sprintf("Per-SNP heritability %g%%", 100 * unique(data_re$pve_g))
  data_re$rho_f <- factor(data_re$rho, names(rho_cols))
  
  plots_ss[[paste0("kinsetting_", kinsetting)]] <-
    ratio_bar_plot(data_re, "rho_f", rho_cols, rho_name, kinsetting_label, ylab_ss, sub_txt)
}

## given a heritability and every GRMr
## for each rho and kinsetting

com_results_grm01 <- read_power_results(GRMr = 0.125,kinsetting = "cousin")
com_results_grm02 <- read_power_results(GRMr = 0.25,kinsetting = "halfsiblings")
com_results_grm03 <- read_power_results(GRMr = 0.5,kinsetting = "siblings")

com_results_grm01$kinsettings <- "cousin"
com_results_grm02$kinsettings <- "halfsiblings"
com_results_grm03$kinsettings <- "siblings"

com_results_all <- rbind(com_results_grm01,com_results_grm02,com_results_grm03)

com_results_all$power_SynSurrG <- 10^{-com_results_all$SynSurrG} < alpha_level
com_results_all$power_SynSurr <- 10^{-com_results_all$SynSurr} < alpha_level

data_re0 <- com_results_all[, c("setting_num","miss","rho","pve_g","power_SynSurrG","power_SynSurr","kinsettings")] %>%
  group_by(miss, rho, pve_g,kinsettings) %>%
  summarise_at(vars(power_SynSurrG, power_SynSurr), mean) %>%
  mutate(ratio = power_SynSurrG / power_SynSurr)

for (rho in c(0,0.25,0.5,0.75)) {
  data_re <- data_re0[which(data_re0$pve_g %in% unique(data_re0$pve_g)[c(ind_pve_g)] & data_re0$rho %in% c(rho)),]
  sub_txt <- sprintf("Per-SNP heritability %g%%", 100 * unique(data_re$pve_g))
  data_re$kin_f <- factor(kin_lab[data_re$kinsettings], kin_lab)
  
  plots_ss[[paste0("rho_", rho)]] <-
    ratio_bar_plot(data_re, "kin_f", kin_cols, NULL, sprintf("\u03c1 = %g", rho), ylab_ss, sub_txt)
}

sel <- data_re0$pve_g %in% unique(data_re0$pve_g)[ind_pve_g]
ymax_ss <- max(data_re0$ratio[sel][is.finite(data_re0$ratio[sel])], 1)
for (nm in names(plots_ss))
  save_ratio(plots_ss[[nm]], ymax = ymax_ss, paste0("power_ratio_Syn-PALM_vs_SynSurr_", nm, ".pdf"))
