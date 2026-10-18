## SuSiE fine-mapping of the primary associations of the no-PRS / no-PC UKB-PPP
## analysis (Fig 5B; replaces code/finemap_results/, which used single-variant
## Wakefield ABFs). For every primary association (lead variant of a +-1 Mb clump
## at P < 5e-8, from post_protein.R) of Syn-PALM and of Obs-GWAS, all autosomal
## genotyped variants within +-1 Mb of the lead variant (MHC chr6:25.5-34.0 Mb
## excluded) are fine-mapped with
## susieR::susie_rss (L = 10, 95% credible sets, minimum purity 0.5) from the
## z-scores beta/SE of that method and an in-sample LD matrix computed from 10,000
## unrelated individuals of the GWAS analysis dataset (seed 20261013; missing
## genotypes mean-imputed). n is the effective sample size implied by the SEs,
## median over the window of 1 / {2p(1-p) SE^2} (the phenotypes are
## inverse-normal transformed, variance 1).
##   Rscript code/post_cf3/finemap_susie.R <task> <n_task> [ncores=1]   per protein
##   Rscript code/post_cf3/finemap_susie.R summarize                    pooled tables
## env FINEMAP_PROTEINS=a,b,c restricts to these proteins.
## Writes results/table_results_cf_noPRS_noPC/susie/<protein>.rds (list with one
## data.table per method: loci, credible sets and variants with PIP >= 0.01) and,
## with summarize, susie_credible_sets.csv / susie_variants_pip.csv / susie_summary.txt.
suppressPackageStartupMessages({ library(data.table); library(susieR) })
setwd("<PROJECT_DIR>")
D   <- "results/table_results_cf_noPRS_noPC/"
OUT <- paste0(D, "susie/"); dir.create(OUT, showWarnings = FALSE)
args <- commandArgs(TRUE)
WIN <- 1e6; N_LD <- 10000L

if (args[1] != "summarize") {
  task <- as.integer(args[1]); ntask <- as.integer(args[2])
  ncores <- if (length(args) >= 3) as.integer(args[3]) else 1L
  ## LD reference panel
  FD <- readRDS("data/crossfit_folds_K5.rds")
  set.seed(20261013)
  ld_eid <- as.character(FD$analysis_eid[sort(sample(FD$independent_indices, N_LD))])
  G0 <- BEDMatrix::BEDMatrix("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
  ld_rows <- match(ld_eid, rownames(G0)); stopifnot(!anyNA(ld_rows))
  files <- sort(list.files(D, "_loci\\.rds$", full.names = TRUE))
  files <- files[seq(task, length(files), by = ntask)]
  if (nzchar(Sys.getenv("FINEMAP_PROTEINS")))
    files <- files[sub("_loci\\.rds$", "", basename(files)) %in% strsplit(Sys.getenv("FINEMAP_PROTEINS"), ",")[[1]]]
  one_locus <- function(lead, S, bcol, vcol, m, p) {
    w <- S[chr == lead$chr & abs(POS38 - lead$POS38) <= WIN & !(chr == "6" & POS38 >= 25.5e6 & POS38 <= 34e6)]
    w <- w[is.finite(get(bcol)) & is.finite(get(vcol)) & get(vcol) > 0]
    G <- G0[ld_rows, w$rsid, drop = FALSE]
    af <- colMeans(G, na.rm = TRUE) / 2
    for (j in which(colSums(is.na(G)) > 0)) G[is.na(G[, j]), j] <- 2 * af[j]
    ok <- af > 0 & af < 1 & apply(G, 2, sd) > 0
    w <- w[ok]; G <- G[, ok, drop = FALSE]; af <- af[ok]
    z <- w[[bcol]] / sqrt(w[[vcol]])
    R <- cor(G)
    n_eff <- median(1 / (2 * af * (1 - af) * w[[vcol]]))
    fit <- tryCatch(susie_rss(z = z, R = R, n = n_eff, L = 10, coverage = 0.95, min_abs_corr = 0.5, max_iter = 500),
                    error = function(e) e)
    base <- data.table(protein = p, method = m, lead_rsid = lead$rsid, chr = lead$chr, lead_pos = lead$POS38,
                       n_snps = nrow(w), n_eff = round(n_eff))
    if (inherits(fit, "error"))
      return(list(locus = base[, `:=`(converged = FALSE, n_cs = NA_integer_, error = conditionMessage(fit))], cs = NULL, var = NULL))
    pip <- fit$pip
    cs <- fit$sets$cs
    cs_dt <- if (length(cs)) rbindlist(lapply(seq_along(cs), function(k) {
      idx <- cs[[k]]; top <- idx[which.max(pip[idx])]
      data.table(protein = p, method = m, lead_rsid = lead$rsid, chr = lead$chr, cs = names(cs)[k], cs_size = length(idx),
                 purity_min_abs_corr = fit$sets$purity[k, "min.abs.corr"], top_rsid = w$rsid[top], top_pos = w$POS38[top],
                 top_pip = pip[top], top_z = z[top], lead_in_cs = lead$rsid %in% w$rsid[idx],
                 variants = paste(w$rsid[idx], collapse = ";"))
    })) else NULL
    keep <- which(pip >= 0.01)
    var_dt <- data.table(protein = p, method = m, lead_rsid = lead$rsid, chr = lead$chr, rsid = w$rsid[keep],
                         POS38 = w$POS38[keep], z = z[keep], pip = pip[keep],
                         cs = vapply(keep, function(i) { k <- which(vapply(cs, function(s) i %in% s, TRUE)); if (length(k)) names(cs)[k[1]] else NA_character_ }, ""))
    list(locus = base[, `:=`(converged = isTRUE(fit$converged), n_cs = length(cs), error = NA_character_,
                             lead_pip = if (lead$rsid %in% w$rsid) pip[match(lead$rsid, w$rsid)] else NA_real_)],
         cs = cs_dt, var = var_dt)
  }
  run_protein <- function(f) {
    L <- readRDS(f); p <- L$protein
    o <- paste0(OUT, p, ".rds"); if (file.exists(o)) return(invisible(NULL))
    S <- as.data.table(readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p)))
    S[, chr := as.character(chr)]
    res <- list()
    for (m in c("SynG", "ObsG")) {
      lo <- L$by_threshold[["5e-08"]][[m]]$loci
      if (!nrow(lo)) next
      lo <- as.data.table(lo)[, chr := as.character(chr)]
      pre <- if (m == "SynG") "SynSurrG" else "ObsG"
      r <- lapply(seq_len(nrow(lo)), function(i) one_locus(lo[i], S, paste0(pre, "_hat_beta"), paste0(pre, "_var_hat_beta"), m, p))
      res[[m]] <- list(loci = rbindlist(lapply(r, `[[`, "locus"), fill = TRUE),
                       cs = rbindlist(lapply(r, `[[`, "cs"), fill = TRUE),
                       var = rbindlist(lapply(r, `[[`, "var"), fill = TRUE))
    }
    attr(res, "rho_oof") <- L$rho_oof
    saveRDS(res, o)
  }
  invisible(parallel::mclapply(files, function(f) tryCatch(run_protein(f), error = function(e) message("FAILED ", f, ": ", conditionMessage(e))),
                               mc.cores = ncores, mc.preschedule = FALSE))
  quit(save = "no")
}

## ---- summary ---------------------------------------------------------------
fs <- list.files(OUT, "\\.rds$", full.names = TRUE)
R <- lapply(fs, readRDS)
grab <- function(m, what) rbindlist(lapply(R, function(x) x[[m]][[what]]), fill = TRUE)
sink(paste0(D, "susie_summary.txt"), split = TRUE)
cat("proteins fine-mapped:", length(R), "\n")
for (m in c("SynG", "ObsG")) {
  lo <- grab(m, "loci"); cs <- grab(m, "cs"); v <- grab(m, "var")
  cat(sprintf("\n[%s] loci %s (failed %d, not converged %d); loci with >= 1 credible set %s (%.1f%%); credible sets %s; per locus: mean %.2f\n",
              m, format(nrow(lo), big.mark = ","), lo[!is.na(error), .N], lo[converged == FALSE & is.na(error), .N],
              format(lo[n_cs > 0, .N], big.mark = ","), 100 * lo[, mean(n_cs > 0, na.rm = TRUE)],
              format(nrow(cs), big.mark = ","), lo[, mean(n_cs, na.rm = TRUE)]))
  if (nrow(cs)) cat(sprintf("  credible set size: median %d, single-variant sets %s (%.1f%%), sets with <= 5 variants %.1f%%; sets containing the lead variant %.1f%%\n",
              as.integer(median(cs$cs_size)), format(cs[cs_size == 1, .N], big.mark = ","), 100 * cs[, mean(cs_size == 1)],
              100 * cs[, mean(cs_size <= 5)], 100 * cs[, mean(lead_in_cs)]))
  if (nrow(v)) cat(sprintf("  variants with PIP >= 0.95: %s; >= 0.5: %s\n", format(v[pip >= 0.95, .N], big.mark = ","), format(v[pip >= 0.5, .N], big.mark = ",")))
  fwrite(cs, sprintf("%ssusie_credible_sets_%s.csv", D, m)); fwrite(v, sprintf("%ssusie_variants_pip_%s.csv.gz", D, m))
  fwrite(lo, sprintf("%ssusie_loci_%s.csv", D, m))
}
sink()
