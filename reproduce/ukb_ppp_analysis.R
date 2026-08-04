## ---------------------------------------------------------------------------
## Provided for inspection, not as a runnable script.
##
## This is the analysis code used to produce the proteome-wide results reported
## in the manuscript. It cannot be executed outside an approved compute
## environment: it requires individual-level UK Biobank phenotype, proteomic and
## genotype data, a pre-computed sparse GRM for the full cohort, and an HPC
## cluster. It is included so that the exact call sequence, input formats,
## covariate set, chunking scheme and output structure can be inspected.
##
## Site-specific directory paths have been replaced with placeholders:
##   <PROJECT_DIR>  project root holding data/, genotype_data/ and results/
##   <GRM_DIR>      directory holding the pre-computed sparse GRM
##   <R_LIBS>       R library path used on the cluster
##   <LAB_DIR>      shared lab directory
##
## Note that this script sourced the analysis functions directly. The same
## functions are now installed as the SynPALM package, so the equivalent line is
##   library(SynPALM)
## together with library(withr), which this script relies on for with_seed().
##
## For a self-contained, runnable demonstration that needs no external data, see
## inst/examples/quickstart.R.
## ---------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(Matrix)
library(ranger)

#####################
setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")

#protein_index #1:2920
args <- commandArgs(trailingOnly = TRUE)

arraynum <- as.integer(args[3])
# The first argument is protein_index (renamed from selected_protein)
protein_index <- as.integer(args[1])
# The second argument is the chunk number
chunk_array <- as.integer(args[2])

cat("Processing protein index:", protein_index, "\n")
cat("Processing chunk array:", chunk_array, "\n")
cat("Processing arraynum:", arraynum, "\n")

### read the grm matrix ################
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm*2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
#obtain the eid
eids_sGRM <- colnames(sgrm)
length(eids_sGRM)

## if sparse
nonzero_count <- nnzero(sgrm)
total_count <- prod(dim(sgrm))
nonzero_count / total_count

### read phenotype and protein data ################
## df_allUKB_eur do not need reorder
df_allUKB_eur <- fread(
  "<PROJECT_DIR>/data/df_allUKB_eur1.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)
## the following need to be reorder
df_all_variables_eur_orignal <- fread(
  "<PROJECT_DIR>/data/extracted_pheno_data528_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)
disease_data_level1_eur_orignal <- fread(
  "<PROJECT_DIR>/data/disease_data_level1_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

all(disease_data_level1_eur_orignal$eid == df_all_variables_eur_orignal$eid)

common_ids <- eids_sGRM[eids_sGRM %in% df_all_variables_eur_orignal$eid]
idx <- match(common_ids, df_all_variables_eur_orignal$eid)
df_all_variables_eur <- df_all_variables_eur_orignal[idx, ]
disease_data_level1_eur <- disease_data_level1_eur_orignal[idx,]

df_allUKB_eur_protein <- df_allUKB_eur[,c(1,49:2968)]

all(disease_data_level1_eur$eid == df_allUKB_eur$eid)
all(df_all_variables_eur$eid == df_allUKB_eur$eid)
all(disease_data_level1_eur$eid == df_allUKB_eur$eid)
all(df_allUKB_eur_protein$eid == df_allUKB_eur$eid)

protein_names <- colnames(df_allUKB_eur_protein)[2:ncol(df_allUKB_eur_protein)]

obs_index <- which(!is.na(df_allUKB_eur_protein[,protein_names[protein_index]]))

G0 <- BEDMatrix::BEDMatrix(path = paste0("<PROJECT_DIR>/genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed",sep = ""), simple_names = T) 
geno_samp.id <- rownames(G0)
matched_indices_prs <- match(df_allUKB_eur_protein$eid, geno_samp.id)

### read the protein PRS model
path_PRS <- "<PROJECT_DIR>/data/proteinPRSmodel"
dirs <- tolower(list.dirs(path_PRS, full.names = FALSE, recursive = FALSE))
#length(dirs)

G_PRS <- data.frame(prs = rep(0,length(df_allUKB_eur_protein$eid)))
rownames(G_PRS) <- df_allUKB_eur_protein$eid

if(protein_names[protein_index] %in% dirs){
  target_path <- file.path(path_PRS, toupper(protein_names[protein_index]))
  txt_files <- list.files(target_path, pattern = "\\.txt$", full.names = TRUE)
  txt_list <- lapply(txt_files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  combined_data <- do.call(rbind, txt_list)
  intercept <- as.numeric(combined_data[1, 3])
  snp_weights_all <- combined_data[-1, ] # 去掉第一行的权重表
  common_snps <- intersect(snp_weights_all[, 1], colnames(G0))
  matched_weights <- snp_weights_all[match(common_snps, snp_weights_all[, 1]), ]
  G0_selected <- G0[matched_indices_prs, common_snps, drop = FALSE]
  G0_selected[is.na(G0_selected)] <- 0
  weights_vector <- as.numeric(matched_weights[, 3])
  G_PRS <- G0_selected %*% weights_vector/10 + intercept/5
  colnames(G_PRS) <- "prs"
}

rm(df_all_variables_eur_orignal,disease_data_level1_eur_orignal)
gc()

all(df_all_variables_eur$eid == rownames(G_PRS))
df_all_variables_eur <- data.frame(df_all_variables_eur,PRS =G_PRS)

### rf prediction for the protein #############
X <- as.matrix(df_all_variables_eur[obs_index,8:(ncol(df_all_variables_eur)-1)])
Y <- as.matrix(df_allUKB_eur_protein[obs_index, -1])
Z <- as.matrix(disease_data_level1_eur[obs_index, -1])

cor_results <- cor(X, Y[,protein_index])
Y_vec <- Y[, protein_index]

m <- ncol(Z)
sig_icd10_log10 <- numeric(m)

sig_icd10_log10 <- vapply(seq_len(m),
                          function(i) {
                            zi <- Z[, i]
                            vals <- unique(zi) 
                            if (length(vals) == 2) {
                              idx1 <- (zi == vals[1])
                              idx0 <- (zi == vals[2])
                              pval <- wilcox.test(
                                x= Y_vec[idx0],y= Y_vec[idx1],alternative = "two.sided",
                                paired= FALSE,exact= FALSE)$p.value
                              return(-log10(pval))
                            } else {
                              return(0)
                            }
                          },
                          numeric(1)
)

train_size <- floor(0.1 * length(obs_index))
sample_for_train <- with_seed(25324, sample(df_allUKB_eur_protein[obs_index,"eid"], 
                                            size = train_size, replace = FALSE))
train_id <- data.frame(eid = sample_for_train)
test_id <- data.frame(eid = setdiff(df_allUKB_eur_protein$eid, train_id$eid))

cov_pred <- c(colnames(X)[order(abs(cor_results),decreasing = T)[1:100]],"p20116_i0","p21000_i0","p31","prs")
if(length(which(sig_icd10_log10>5) <= 5)){
  cov_dis <- colnames(Z)[order(abs(sig_icd10_log10),decreasing = T)[1:5]]
}else{
  cov_dis <- colnames(Z)[order(abs(sig_icd10_log10),decreasing = T)[1:10]]
}

train_rf1 <- (df_all_variables_eur %>% filter(eid %in% train_id$eid))[,cov_pred]
train_rf2 <- (disease_data_level1_eur %>% filter(eid %in% train_id$eid))[,cov_dis]
train_rf2[] <- lapply(train_rf2, function(x) 
  factor(x, levels = c(0, 1), labels = c("0", "1"))
)
train_rf <- cbind(train_rf1, train_rf2)
cols_idx <- grep("[ -]", names(train_rf))
old_names <- names(train_rf)[cols_idx]
new_names <- gsub("[- ]", "_", old_names)
names(train_rf)[cols_idx] <- new_names

train_rf$protein <- (df_allUKB_eur_protein %>% filter(eid %in% train_id$eid))[,protein_names[protein_index]]

test_rf1 <- (df_all_variables_eur %>% filter(eid %in% test_id$eid))[,cov_pred]
#test_rf2 <- (disease_data_level0_eur %>% filter(eid %in% test_id$eid))[,cov_dis]
test_rf2 <- (disease_data_level1_eur %>% filter(eid %in% test_id$eid))[,cov_dis]
test_rf2[] <- lapply(test_rf2, function(x) 
  factor(x, levels = c(0, 1), labels = c("0", "1"))
)
test_rf <- cbind(test_rf1,test_rf2)
cols_idx <- grep("[ -]", names(test_rf))
old_names <- names(test_rf)[cols_idx]
new_names <- gsub("[- ]", "_", old_names)
names(test_rf)[cols_idx] <- new_names

#rf model
rf_model <- ranger(protein ~ ., data = train_rf,num.trees = 300, importance = "impurity", 
                   keep.inbag = TRUE, save.memory = FALSE,write.forest = TRUE,seed = 260611)

rf_pred <- predict(rf_model, data = test_rf)$predictions

#cov_id <- c("age",paste0("geneticPC",1:20),"sex","batch","UKBcentre")
cov_id <- c("p21022",paste0("p22009_a",1:20),"p31","p22000")
test_df1 <- (df_all_variables_eur %>% filter(eid %in% test_id$eid))[,c("eid",cov_id)]
colnames(test_df1) <- c("eid","age",paste0("geneticPC",1:20),"sex","batch")

test_df2 <- test_id %>% left_join(
  df_allUKB_eur %>% dplyr::select(eid, UKBcentre),by = "eid")

test_df <- test_df1 %>% left_join(test_df2,by = "eid")

test_df$age2 <- test_df$age^2
test_df$age_sex <- test_df$age*(as.numeric(test_df$sex)-1)
test_df$age2_sex <- test_df$age2*(as.numeric(test_df$sex)-1)

test_df$protein <- (df_allUKB_eur_protein %>% filter(eid %in% test_id$eid))[,protein_names[protein_index]]
test_df$protein_hat <- rf_pred

cor_results <- cor(test_df$protein[which(!is.na(test_df$protein))],test_df$protein_hat[which(!is.na(test_df$protein))])

saveRDS(cor_results, file = paste0("<PROJECT_DIR>/results/allukb_protein_pQTL_european_discovery/",
                                   protein_names[protein_index],"/",protein_names[protein_index],"_correlation_more_vars_withICD_level1_noPRS_select_for_each_protein.rds"))

saveRDS(length(obs_index), file = paste0("<PROJECT_DIR>/results/allukb_protein_pQTL_european_discovery/",
                                         protein_names[protein_index],"/",protein_names[protein_index],"_number_of_observed_protein.rds"))


rm(df_all_variables_eur,df_allUKB_eur,df_allUKB_eur_protein,disease_data_level1_eur)
gc()

#### using reml to estimate the variances components ###############

## inverse normal transform of the outcome Y and Yhat
test_df <- INT(test_df, "protein") 
n <- nrow(test_df)
r <- rank(test_df %>% pull(protein_hat))
test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat

test_df_dmm <- model.matrix(~ sex + batch + UKBcentre, data = test_df)
continuous_vars <- test_df[, setdiff(colnames(test_df),c("sex", "batch", "UKBcentre"))]
test_df_transformed <- cbind(continuous_vars, test_df_dmm)
test_df_transformed$`(Intercept)` <- NULL # fitting model without Intercept
#head(test_df_transformed)
colnames(test_df_transformed) <- gsub("-", "_", colnames(test_df_transformed))

## obtain the submatrix of GRM and obtain the independent samples
index_all <- match(as.character(test_df_transformed$eid), eids_sGRM)
GRM <- sgrm[index_all, index_all] #398318 not excetly the same for each protein

hh <- find_blocks_vectorized(GRM)
independent_indices <- with_seed(25324, {
  sapply(hh, function(x) x[sample.int(length(x), 1)])
}) #256148

obs_protein_index <- which(!is.na(test_df_transformed$protein))

##### trans test_df_transformed to mydf form 
X_all_names <- setdiff(colnames(test_df_transformed), 
                       c("eid", "protein", "protein_hat","int", "yhat_int"))

mydf <- list(X_all = test_df_transformed[,X_all_names],
             S = test_df_transformed$yhat_int,
             Y_obs = test_df_transformed$int,
             GRM = GRM)

SynSurrG_step1_pars <- SynSurrG_ablation_estimate(mydf)
ObsG_step1_pars <- ObsG_ablation_estimate(mydf)

SynSurr_step1_pars <- SynSurr_ablation_estimate(mydf,independent_indices)
Obs_step1_pars <- Obs_ablation_estimate(mydf,independent_indices)

############# v2 genotypes ##################


cat("Processing v2geno protein name", protein_names[protein_index], "chunk array:", chunk_array, "\n")

rs_index <- colnames(G0)

# Split the vector into the specified number of sets
rs_index_list <- split(rs_index, cut(seq_along(rs_index), breaks = arraynum, labels = FALSE))

chunk_size <- 200  # Load 200 SNPs at a time
num_chunks <- ceiling(length(rs_index_list[[chunk_array]]) / chunk_size)  # Total chunks

SynSurrG_results <- list()
ObsG_results <- list()

SynSurr_results <- list()
Obs_results <- list()

for (chunk_index in seq_len(num_chunks)) {
  
  matched_indices <- match(test_df_transformed$eid, geno_samp.id)
  if (any(is.na(matched_indices))) {
    stop("Some eids in test_df are not found in G.")
  }
  
  chunk_start <- (chunk_index - 1) * chunk_size + 1
  chunk_end <- min(chunk_index * chunk_size, length(rs_index_list[[chunk_array]]))  # Ensure last chunk does not exceed length
  chunk_idx <- rs_index_list[[chunk_array]][chunk_start:chunk_end]  # Get SNP indices for this chunk
  
  # Read only the SNPs in this chunk and store in list
  Genotype_matrix <- as.matrix(G0[matched_indices, chunk_idx])
  Genotype_matrix <- replace(Genotype_matrix, is.na(Genotype_matrix), 0)
  Genotype_matrix <- fix_constant_columns(Genotype_matrix, intersect(obs_protein_index,independent_indices))
  
  SynSurrG_results[[chunk_index]] <- score_test_SynSurrG_multiply(g_matrix = Genotype_matrix, step1_pars = SynSurrG_step1_pars)
  ObsG_results[[chunk_index]] <- score_test_ObsG_multiply(g_matrix = Genotype_matrix,step1_pars = ObsG_step1_pars)
  
  SynSurr_results[[chunk_index]] <- score_test_SynSurr_multiply(g_matrix = Genotype_matrix,step1_pars = SynSurr_step1_pars,independent_indices)
  Obs_results[[chunk_index]] <- score_test_Obs_multiply(g_matrix = Genotype_matrix,step1_pars = Obs_step1_pars,independent_indices)
  
  SynSurrG_results[[chunk_index]] <- lapply(SynSurrG_results[[chunk_index]], function(vec) {
    names(vec) <- colnames(Genotype_matrix)
    vec
  })
  ObsG_results[[chunk_index]] <- lapply(ObsG_results[[chunk_index]], function(vec) {
    names(vec) <- colnames(Genotype_matrix)
    vec
  })
  SynSurr_results[[chunk_index]] <- lapply(SynSurr_results[[chunk_index]], function(vec) {
    names(vec) <- colnames(Genotype_matrix)
    vec
  })
  Obs_results[[chunk_index]] <- lapply(Obs_results[[chunk_index]], function(vec) {
    names(vec) <- colnames(Genotype_matrix)
    vec
  })
  
  if (chunk_index %% 5 == 0) {
    print(chunk_index)
  }
}

merged_SynSurrG <- merge_results(SynSurrG_results)
merged_ObsG     <- merge_results(ObsG_results)
merged_SynSurr  <- merge_results(SynSurr_results)
merged_Obs      <- merge_results(Obs_results)

SynSurrG_negative_log10_pvalues <- merged_SynSurrG$negative_log10_pval_SynSurrG
SynSurrG_hat_beta <- merged_SynSurrG$hat_beta_SynSurrG
SynSurrG_var_hat_beta <- merged_SynSurrG$var_hat_beta_SynSurrG

ObsG_negative_log10_pvalues <- merged_ObsG$negative_log10_pval_ObsG
ObsG_hat_beta <- merged_ObsG$hat_beta_ObsG
ObsG_var_hat_beta <- merged_ObsG$var_hat_beta_ObsG

SynSurr_negative_log10_pvalues <- merged_SynSurr$negative_log10_pval_SynSurr
SynSurr_hat_beta <- merged_SynSurr$hat_beta_SynSurr
SynSurr_var_hat_beta <- merged_SynSurr$var_hat_beta_SynSurr

Obs_negative_log10_pvalues <- merged_Obs$negative_log10_pval_Obs
Obs_hat_beta <- merged_Obs$hat_beta_Obs
Obs_var_hat_beta <- merged_Obs$var_hat_beta_Obs

results <- cbind(SynSurrG_negative_log10_pvalues,
                 SynSurrG_hat_beta,SynSurrG_var_hat_beta,
                 ObsG_negative_log10_pvalues,
                 ObsG_hat_beta,ObsG_var_hat_beta,
                 SynSurr_negative_log10_pvalues,
                 SynSurr_hat_beta,SynSurr_var_hat_beta,
                 Obs_negative_log10_pvalues,
                 Obs_hat_beta,Obs_var_hat_beta)

saveRDS(results, file = paste0("<PROJECT_DIR>/results/allukb_protein_pQTL_european_discovery/",
                               protein_names[protein_index],"/",protein_names[protein_index],
                               "_chunkarray",chunk_array, "_negative_log10_pvals_v2geno_noPRS.rds"))

print(paste("Protein",protein_names[protein_index],"chunk_array",chunk_array, "v2geno processing completed."))
