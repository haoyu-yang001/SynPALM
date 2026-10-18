## Proteins to rerun after the 2026-10-13 VC fix in funcs_v2.R (residual
## components re-solved after the genetic block is projected).
## The fix changes results only when the genetic block was projected
## (vc$projected_G) and the scan was run before the fix (no vc$rematch_sigma).
##   Rscript code/post_cf3/list_vc_rerun.R [ncores=4]
## Writes code/sh_run/run_noPRS_noPC/order_vcfix_main.txt  (lines of order_cor_vec7_new.txt)
##        code/sh_run/run_noPRS_noPC/order_vcfix_split.txt (lines of order_split_top100.txt)
## and results/vc_check_cf_noPRS_noPC.csv.
suppressPackageStartupMessages({ library(data.table); library(parallel) })
setwd("<PROJECT_DIR>")
NC  <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 4L
SH  <- "code/sh_run/run_noPRS_noPC/"
RES <- "results/allukb_protein_pQTL_european_discovery/"
names_all <- strsplit(readLines("data/df_allUKB_eur1.tsv", n = 1), "\t")[[1]]
pidx  <- as.integer(readLines("code/sh_run/order_cor_vec7_new.txt"))
pname <- names_all[48 + pidx]

vc_row <- function(vc) data.table(
  projected_G = isTRUE(vc$projected_G), projected_E = isTRUE(vc$projected_E),
  rematched = isTRUE(vc$rematch_sigma),
  raw_tau_T2 = vc$raw[["tau_T2"]], raw_tau_TS = vc$raw[["tau_TS"]], raw_tau_S2 = vc$raw[["tau_S2"]],
  raw_sigma_T2 = vc$raw[["sigma_T2"]])

## main scans
main <- rbindlist(mclapply(seq_along(pname), function(i) {
  f <- sprintf("%s%s/%s_cf_noPRS_noPC_scan_allSNPs.rds", RES, pname[i], pname[i])
  if (!file.exists(f)) return(data.table(set = "main", line = i, protein = pname[i], done = FALSE))
  cbind(data.table(set = "main", line = i, protein = pname[i], done = TRUE), vc_row(attr(readRDS(f), "vc")))
}, mc.cores = NC), fill = TRUE)

## split scans (vc holds one entry per half)
sp_idx <- as.integer(readLines(paste0(SH, "order_split_top100.txt")))
split <- rbindlist(mclapply(seq_along(sp_idx), function(i) {
  p <- names_all[48 + sp_idx[i]]
  f <- sprintf("results/split_validation_cf_noPRS_noPC/%s_split_scan.rds", p)
  if (!file.exists(f)) return(data.table(set = "split", line = i, protein = p, done = FALSE))
  vc <- attr(readRDS(f), "vc")
  r  <- rbindlist(lapply(vc, function(v) vc_row(if (!is.null(v$SynSurrG)) v$SynSurrG else v)))
  data.table(set = "split", line = i, protein = p, done = TRUE,
             projected_G = any(r$projected_G), projected_E = any(r$projected_E),
             rematched = all(r$rematched), raw_tau_T2 = min(r$raw_tau_T2),
             raw_tau_TS = NA_real_, raw_tau_S2 = min(r$raw_tau_S2), raw_sigma_T2 = min(r$raw_sigma_T2))
}, mc.cores = NC), fill = TRUE)

all <- rbind(main, split, fill = TRUE)
fwrite(all, "results/vc_check_cf_noPRS_noPC.csv")
rr <- function(s) all[set == s & done & projected_G & !rematched, line]
writeLines(as.character(pidx[rr("main")]),    paste0(SH, "order_vcfix_main.txt"))
writeLines(as.character(sp_idx[rr("split")]), paste0(SH, "order_vcfix_split.txt"))
for (s in c("main", "split"))
  cat(sprintf("%s: %d done, %d projected_G, %d to rerun (not yet run with the fix)\n", s,
              all[set == s & done == TRUE, .N], all[set == s & done & projected_G, .N], length(rr(s))))
