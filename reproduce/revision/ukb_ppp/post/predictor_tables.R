## Supplementary Data for Reviewer 2 comment 10: candidate predictors of the
## random forest (with UK Biobank field IDs) and the predictors selected for
## every protein in every cross-fitting fold (step1_protein_crossfit.R, no PRS /
## no genetic PCs).
##   Rscript code/post_cf3/predictor_tables.R
## Writes results/table_results_cf_noPRS_noPC/:
##   predictors_candidates.csv          all candidates: variable, UKB field ID, type, description
##   predictors_selected_by_fold.csv.gz protein, fold, variable, kind (continuous / fixed / ICD-10), rank
##   predictors_selection_frequency.csv share of protein-folds in which each candidate was selected
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
OUT <- "results/table_results_cf_noPRS_noPC/"
TAG <- "_cf_noPRS_noPC"
hdr  <- names(fread("data/extracted_pheno_data528_eur.tsv", nrows = 0))
cand <- setdiff(hdr[8:length(hdr)], "prs"); cand <- cand[!grepl("^p22009_", cand)]
fixed <- c("p20116_i0", "p21000_i0", "p31")
icd  <- setdiff(names(fread("data/disease_data_level1_eur.tsv", nrows = 0)), "eid")
dict <- fread("data/Predictive_Lab_Values528.csv")
field <- function(v) ifelse(grepl("^p\\d+", v), sub("^p(\\d+).*", "\\1", v), NA)
cands <- rbind(
  data.table(variable = cand,  kind = "continuous (screened by |r|, top 100 per fold)"),
  data.table(variable = fixed, kind = "always included"),
  data.table(variable = icd,   kind = "ICD-10 level-1 category (screened by Wilcoxon, top 10 per fold)"))
cands[, ukb_field_id := field(variable)]
cands[dict, on = .(variable = variable_name), `:=`(type = i.Type, description = i.Description)]
cands[kind != "always included" & startsWith(kind, "ICD"), `:=`(type = "binary (any ICD-10 diagnosis in the block)",
                                                                  description = variable)]
fwrite(cands, paste0(OUT, "predictors_candidates.csv"))

prot <- names(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2964]
sel <- rbindlist(lapply(prot, function(p) {
  f <- sprintf("results/allukb_protein_pQTL_european_discovery/%s/%s_selected_features_by_fold%s.rds", p, p, TAG)
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  rbindlist(lapply(seq_along(x), function(k) {
    cp <- x[[k]]$cov_pred
    rbind(data.table(protein = p, fold = k, variable = setdiff(cp, fixed),
                     kind = "continuous", rank = seq_along(setdiff(cp, fixed))),
          data.table(protein = p, fold = k, variable = fixed, kind = "fixed", rank = NA_integer_),
          data.table(protein = p, fold = k, variable = x[[k]]$cov_dis, kind = "ICD-10", rank = seq_along(x[[k]]$cov_dis)))
  }))
}))
stopifnot(all(sel[kind == "continuous", variable] %in% cand), all(sel[kind == "ICD-10", variable] %in% icd))
fwrite(sel, paste0(OUT, "predictors_selected_by_fold.csv.gz"))
nf <- uniqueN(sel[, .(protein, fold)])
fq <- sel[, .(share_of_protein_folds = .N / nf), by = .(variable, kind)][order(-share_of_protein_folds)]
fq[cands, on = "variable", `:=`(ukb_field_id = i.ukb_field_id, description = i.description)]
fwrite(fq, paste0(OUT, "predictors_selection_frequency.csv"))
cat("candidates:", nrow(cands), "(continuous", length(cand), ", ICD-10", length(icd), ")\n",
    "proteins:", uniqueN(sel$protein), " protein-folds:", nf, "\n")
print(head(fq, 15))
