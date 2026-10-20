## Reply to Reviewer 1, minor comment (f) (revision 2026-10-20): are the labelled
## individuals (UKB-PPP randomised baseline subcohort, after the exclusion of the
## consortium-selected and COVID-19 imaging participants) a random subsample of the
## 402,086 analysed individuals? Standardised mean differences (SMD) of baseline
## characteristics between labelled (GDF15 measured) and unlabelled individuals, the share
## of individuals with a relative in the cohort, and the SMD of the synthetic phenotype of
## GDF15 and LEP.
##   Rscript code/post_cf3/labelled_balance.R
## Writes results/table_results_cf_noPRS_noPC/labelled_balance.csv
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
FD <- readRDS("data/crossfit_folds_K5.rds")
v  <- c("age", "sex", "bmi", "height", "sbp", "chol", "HbA1c", "TDI", "smoking", "alcohol", "geneticPC1", "geneticPC2", "UKBcentre")
x  <- fread("data/df_allUKB_eur1.tsv", select = c("eid", v, "gdf15", "lep"))
x  <- x[match(FD$analysis_eid, eid)]
x[, labelled := !is.na(gdf15)]
rel <- integer(length(FD$fold)); for (b in FD$blocks) if (length(b) > 1) rel[b] <- 1L
x[, has_relative := rel]
for (p in c("gdf15", "lep")) {
  s1 <- readRDS(sprintf("results/allukb_protein_pQTL_european_discovery/%s/%s_synthetic_phenotype_cf_noPRS_noPC.rds", p, p))
  stopifnot(identical(as.character(s1$eid), as.character(FD$analysis_eid)))
  x[, (paste0("Yhat_", p)) := s1$protein_hat]
}
num <- c(setdiff(v, c("sex", "smoking", "UKBcentre")), "has_relative", "Yhat_gdf15", "Yhat_lep")
smd <- function(a, b) (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
res <- rbindlist(lapply(num, function(k) {
  a <- as.numeric(x[labelled == TRUE][[k]]); b <- as.numeric(x[labelled == FALSE][[k]])
  data.table(variable = k, labelled = mean(a, na.rm = TRUE), unlabelled = mean(b, na.rm = TRUE), smd = smd(a, b)) }))
## categorical: sex (male share), smoking categories, centre (max |difference| in share)
cat_smd <- function(k) {
  t <- prop.table(table(x[[k]], x$labelled), 2); d <- t[, "TRUE"] - t[, "FALSE"]
  p <- (t[, "TRUE"] + t[, "FALSE"]) / 2
  data.table(variable = paste0(k, " (largest category difference)"), labelled = NA_real_, unlabelled = NA_real_,
             smd = d[which.max(abs(d))] / sqrt(p[which.max(abs(d))] * (1 - p[which.max(abs(d))]))) }
res <- rbind(res, cat_smd("sex"), cat_smd("smoking"), cat_smd("UKBcentre"))
res <- rbind(data.table(variable = "N", labelled = sum(x$labelled), unlabelled = sum(!x$labelled), smd = NA), res)
fwrite(res, "results/table_results_cf_noPRS_noPC/labelled_balance.csv")
print(res, digits = 3)
cat(sprintf("max |SMD| over the baseline characteristics: %.3f\n", max(abs(res$smd), na.rm = TRUE)))
