## Summaries and figures of S2 (alphaG_rev1.R) and S7 (variance_components_rev1.R).
##   Rscript analyze_alphaG_vc_rev1.R
## Writes to simulation_code/results_rev1/figs/:
##   S2_typeIerror.csv, S2_qq.{pdf,png}      type I error of SNPs that affect only S
##   S7_vc_summary.csv, S7_vc_boxplot_*.pdf, S7_time.{csv,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
R <- "simulation_code/results_rev1/"; FIG <- paste0(R, "figs/"); dir.create(FIG, showWarnings = FALSE)
arr_ok <- function(f) as.integer(sub(".*_array_(\\d+)\\.rds", "\\1", f)) <= 500

## ---- S2 -------------------------------------------------------------------
fs <- list.files(paste0(R, "S2"), "\\.rds$", full.names = TRUE); fs <- fs[arr_ok(fs)]
if (length(fs)) {
  d <- rbindlist(lapply(fs, function(f) { x <- as.data.table(readRDS(f)); s <- attr(readRDS(f), "settings")
    x[, `:=`(setting = as.integer(sub(".*setting_(\\d+)_.*", "\\1", f)), miss = s$miss, h_a = s$h_a, K = s$K,
             inv = s$invariance_max_abs_diff)] }))
  m <- melt(d, id.vars = c("setting", "miss", "h_a", "K", "type", "rep", "snp", "inv"),
            measure.vars = c("SynSurrG", "SynSurrG_refit", "ObsG", "SynSurr"), variable.name = "method", value.name = "nlp")
  m <- m[is.finite(nlp)]
  al <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
  t1 <- m[, c(list(n = .N, lambda_gc = median(qchisq(10^-nlp, 1, lower.tail = FALSE)) / qchisq(0.5, 1)),
              setNames(lapply(al, function(a) mean(nlp >= -log10(a))), formatC(al, format = "e", digits = 0))),
          by = .(setting, miss, h_a, K, type, method)]
  setorder(t1, setting, type, method)
  t1[, max_invariance_diff := d[, max(inv, na.rm = TRUE), by = setting][match(t1$setting, setting), V1]]
  fwrite(t1, paste0(FIG, "S2_typeIerror.csv")); print(t1, digits = 3)
  vcs <- d[snp == 1 & type == "alpha", .(tau_T2 = mean(tau_T2), tau_TS = mean(tau_TS), tau_S2 = mean(tau_S2)),
           by = .(setting, miss, h_a, K)]
  fwrite(vcs, paste0(FIG, "S2_variance_components.csv"))
  q <- m[type == "alpha" & method %in% c("SynSurrG", "ObsG")]
  q <- q[order(-nlp), .(observed = nlp, expected = -log10((seq_len(.N) - 0.5) / .N), N = .N, r = seq_len(.N)),
         by = .(setting, miss, h_a, K, method)]
  q <- q[r <= 2000 | r %% max(1, floor(N / 3000)) == 0]
  ## common look (theme_rev1.R): two strip rows (missing rate; h_alpha and K), red identity
  ## line, grey 95% band with dashed borders (pointwise, from the ranks of the Syn-PALM points)
  source("simulation_code/rev1/theme_rev1.R")
  num_pm <- function(x) { e <- floor(log10(x) + 1e-9); m <- signif(x / 10^e, 3)   # plotmath of a number
    ifelse(m == 1, sprintf("10^%d", e), sprintf("%g %%*%% 10^%d", m, e)) }
  q[, method := factor(c(SynSurrG = "Syn-PALM", ObsG = "Obs-GWAS")[as.character(method)], c("Syn-PALM", "Obs-GWAS"))]
  q[, miss_p := factor(sprintf('"Missing rate %g%%"', 100 * miss), sprintf('"Missing rate %g%%"', 100 * sort(unique(miss))))]
  q[, hk_p := sprintf('italic(h)[alpha] == %s * "," ~ italic(K) == "%s"', num_pm(h_a), sapply(K, format, big.mark = ",", scientific = FALSE))]
  q[, hk_p := factor(hk_p, unique(hk_p[order(h_a)]))]
  band <- q[method == "Syn-PALM", .(expected, clower = -log10(qbeta(0.025, r, N - r + 1)),
                                    cupper = -log10(qbeta(0.975, r, N - r + 1))), by = .(miss_p, hk_p)]
  p <- ggplot(q, aes(expected, observed)) + qq_ref_layers(band) +
    geom_point(aes(colour = method, shape = method), size = 1, stroke = 0.25) +
    facet_wrap(~ miss_p + hk_p, scales = "free", labeller = label_parsed, nrow = 2) +
    scale_colour_manual(values = method_cols, name = NULL) +
    scale_shape_manual(values = c("Syn-PALM" = 16, "Obs-GWAS" = 2), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.4))) +
    scale_x_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.03)) +
    labs(x = qq_xlab, y = qq_ylab) +
    theme_rev1(12, legend = "bottom") + theme(strip.text = element_text(size = rel(0.85)))
  ggsave(paste0(FIG, "S2_qq.pdf"), p, width = 7.4, height = 5.6, device = cairo_pdf)
  ggsave(paste0(FIG, "S2_qq.png"), p, width = 7.4, height = 5.6, dpi = 300)
}

## ---- S7 -------------------------------------------------------------------
## completed arrays plus the (replicate, rho) parts of arrays that are still running
## (variance_components_rev1.R saves every part as soon as it is fitted)
fs <- list.files(paste0(R, "variance_components"), "\\.rds$", full.names = TRUE); fs <- fs[arr_ok(fs)]
fp <- list.files(paste0(R, "variance_components/parts"), "\\.rds$", full.names = TRUE)
if (length(fs) + length(fp)) {
  v <- unique(rbindlist(lapply(c(fs, fp), function(f) as.data.table(readRDS(f)))),
              by = c("miss", "rho", "array", "rep", "method", "param"))
  print(dcast(unique(v[, .(miss, rho, array, rep)])[, .N, by = .(miss, rho)], miss ~ rho, value.var = "N"))
  v[, method := factor(method, c("MoM", "REML", "BOLT-REML"))]
  sm <- v[, .(n_rep = .N, truth = truth[1], mean = mean(estimate), bias = mean(estimate - truth),
              sd = sd(estimate), mse = mean((estimate - truth)^2)), by = .(miss, rho, param, method)]
  fwrite(sm, paste0(FIG, "S7_vc_summary.csv"))
  source("simulation_code/rev1/theme_rev1.R")
  par_lab <- c(tau_T2 = "tau[T]^2", tau_S2 = "tau[S]^2", tau_TS = "tau[TS]",
               sigma_T2 = "sigma[T]^2", sigma_S2 = "sigma[S]^2", sigma_TS = "sigma[TS]")
  v[, par := factor(par_lab[param], par_lab)]
  for (r in sort(unique(v$rho))) for (mi in sort(unique(v$miss))) {
    d <- v[rho == r & miss == mi]
    if (!nrow(d)) next
    nrep <- uniqueN(d[, .(array, rep)])
    p <- ggplot(d, aes(method, estimate, fill = method)) +
      geom_boxplot(width = 0.62, outlier.size = 0.5, outlier.alpha = 0.5, linewidth = 0.35, colour = "grey20") +
      geom_hline(aes(yintercept = truth), linetype = "dashed", colour = "#C0392B", linewidth = 0.5) +
      facet_wrap(~ par, scales = "free_y", nrow = 2, labeller = label_parsed) +
      scale_fill_manual(values = method_cols[c("MoM", "REML", "BOLT-REML")]) +
      labs(x = NULL, y = "Estimate",
           title = bquote("Missing rate" ~ .(100 * mi) * "%," ~ rho == .(r)),
           subtitle = sprintf("%d replicates; dashed line: true value", nrep)) +
      theme_rev1(10, legend = "none") + theme(panel.grid.major.x = element_blank())
    ggsave(sprintf("%sS7_vc_boxplot_missing_%s_rho_%s.pdf", FIG, mi, r), p, width = 6.4, height = 4.6, device = cairo_pdf)
  }
  tm <- unique(v[, .(miss, rho, array, rep, method, seconds, n)])[, .(n_rep = .N, mean_seconds = mean(seconds), sd_seconds = sd(seconds),
                                                                      median_seconds = median(seconds), n = n[1]),
                                                                   by = .(miss, rho, method)]
  fwrite(tm, paste0(FIG, "S7_time.csv")); print(tm, digits = 3)
  ## points (mean) and lines per method; a log-scale bar chart would draw sub-second bars downwards
  p <- ggplot(tm, aes(factor(100 * miss), mean_seconds, colour = method, shape = method, group = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_wrap(~ rho, nrow = 1, labeller = label_bquote(rho == .(rho))) +
    scale_y_log10(breaks = 10^(-2:5), labels = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000", "100,000")) +
    scale_colour_manual(values = method_cols[c("MoM", "REML", "BOLT-REML")], name = NULL) +
    scale_shape_manual(values = c(MoM = 16, REML = 15, `BOLT-REML` = 17), name = NULL) +
    labs(x = "Missing rate (%)", y = "Seconds per fit (log scale)") +
    theme_rev1(10, legend = "bottom")
  ggsave(paste0(FIG, "S7_time.pdf"), p, width = 8.4, height = 3.5, device = cairo_pdf)
}
cat("done\n")
