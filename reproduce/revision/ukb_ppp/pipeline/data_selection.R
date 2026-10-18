library(data.table)
library(dplyr)
library(Matrix)

df_allUKB_eur <- fread(
  "<PROJECT_DIR>/data/df_allUKB_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

df_all_variables <- fread(
  "<PROJECT_DIR>/data/extracted_pheno_data528.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)
colnames(df_all_variables) <- gsub("^participant\\.", "", colnames(df_all_variables))

df_all_variables_names <- fread(
  "<PROJECT_DIR>/data/Predictive_Lab_Values528.csv",
  sep       = ",", header    = TRUE,data.table = FALSE
)

# in df_all_variables_names, add a column for the correlation between each variable

ids <- df_allUKB_eur[,1]
ix <- match(ids, df_all_variables[,1])
df_extracted_pheno_data528_sub <- df_all_variables[ix, ]

# check
all(df_extracted_pheno_data528_sub[,1] == df_allUKB_eur[,1])

## impute the missing values, 
aa <- apply(df_extracted_pheno_data528_sub,2,function(x){length(which(is.na(x)))})

## remove the vauralve that are all missing
which(aa == nrow(df_extracted_pheno_data528_sub))

df_all_variables_names <- df_all_variables_names[-which(aa == nrow(df_extracted_pheno_data528_sub)),]
df_extracted_pheno_data528_sub <- df_extracted_pheno_data528_sub[,-which(aa == nrow(df_extracted_pheno_data528_sub))]

all(colnames(df_extracted_pheno_data528_sub) == df_all_variables_names[,1])

## using the median to impute
for (i in 1:ncol(df_extracted_pheno_data528_sub)) {
  df_extracted_pheno_data528_sub[which(is.na(df_extracted_pheno_data528_sub[,i])),i] <- median(as.numeric(df_extracted_pheno_data528_sub[,i]),na.rm = T)
}

df_allUKB_eur_protein <- df_allUKB_eur[,c(1,49:2968)]

obs_index <- which(!is.na(df_allUKB_eur_protein[,2]))

####caulate the correlation between each protein and each variable
X <- as.matrix(df_extracted_pheno_data528_sub[obs_index,8:ncol(df_extracted_pheno_data528_sub)])
Y <- as.matrix(df_allUKB_eur_protein[obs_index, -1])

# this computes cor(X[,i], Y[,j]) for all i,j
cor_results <- cor(X, Y)


h1 <- apply(cor_results,1,mean)

df_all_variables_names <- df_all_variables_names[-(which(is.na(h1))+7),]
df_extracted_pheno_data528_sub <- df_extracted_pheno_data528_sub[,-(which(is.na(h1))+7)]

all(colnames(df_extracted_pheno_data528_sub) == df_all_variables_names[,1])

X <- as.matrix(df_extracted_pheno_data528_sub[obs_index,8:ncol(df_extracted_pheno_data528_sub)])
Y <- as.matrix(df_allUKB_eur_protein[obs_index, -1])

# this computes cor(X[,i], Y[,j]) for all i,j
cor_results <- cor(X, Y)

h1 <- apply(cor_results,1,mean)

df_all_variables_names <- cbind(df_all_variables_names,c(rep("NA",7),h1))
colnames(df_all_variables_names)[4] <- "cor"


write.table(
  df_all_variables_names,
  file = "<PROJECT_DIR>/data/Predictive_Lab_Values528_cor.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

#### read disease data ######################

df_ICD10_level0coding_info <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/df_ICD10_long_level0_info.tsv",data.table = F)
df_ICD10_level0coding_info[,1] <- gsub(" ", "_", df_ICD10_level0coding_info[,1])
disease_data_level0 <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/extracted_data_ICD10_20231220_matrix_level0.txt",data.table = F)

df_ICD10_level1coding_info <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/df_ICD10_long_level1_info.tsv",data.table = F)
df_ICD10_level1coding_info[,1] <- gsub(" ", "_", df_ICD10_level1coding_info[,1])
disease_data_level1 <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/extracted_data_ICD10_20231220_matrix_level1.txt",data.table = F)

df_ICD10_level2coding_info <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/df_ICD10_long_level2_info.tsv",data.table = F)
df_ICD10_level2coding_info[,1] <- gsub(" ", "_", df_ICD10_level2coding_info[,1])
disease_data_level2 <- fread("<UKB_DIR>/phenotypes/ICD_Phenos/extracted_data_ICD10_20231220_matrix_level2.txt",data.table = F)


ids <- df_allUKB_eur[,1]
ix <- match(ids, disease_data_level0[,1])
disease_data_level0_eur <- disease_data_level0[ix, ]
# check
all(disease_data_level0_eur[,1] == df_allUKB_eur[,1])

ix <- match(ids, disease_data_level1[,1])
disease_data_level1_eur <- disease_data_level1[ix, ]
# check
all(disease_data_level1_eur[,1] == df_allUKB_eur[,1])

ix <- match(ids, disease_data_level2[,1])
disease_data_level2_eur <- disease_data_level2[ix, ]
# check
all(disease_data_level2_eur[,1] == df_allUKB_eur[,1])

write.table(
  disease_data_level0_eur,
  file = "<PROJECT_DIR>/data/disease_data_level0_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  disease_data_level1_eur,
  file = "<PROJECT_DIR>/data/disease_data_level1_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  disease_data_level2_eur,
  file = "<PROJECT_DIR>/data/disease_data_level2_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

############ re-order samples ####################

library(data.table)
library(dplyr)
library(Matrix)

df_allUKB_eur_protein <- fread(
  "<PROJECT_DIR>/data/df_allUKB_eur_protein.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

obs_index <- which(!is.na(df_allUKB_eur_protein[,2]))
miss_index <- which(is.na(df_allUKB_eur_protein[,2]))

df_allUKB_eur_protein <- df_allUKB_eur_protein[c(obs_index,miss_index),]

df_extracted_pheno_data528_sub <- fread(
  "<PROJECT_DIR>/data/extracted_pheno_data528_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

df_extracted_pheno_data528_sub <- df_extracted_pheno_data528_sub[c(obs_index,miss_index),]

disease_data_level2_eur <- fread(
  "<PROJECT_DIR>/data/disease_data_level2_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

disease_data_level2_eur <- disease_data_level2_eur[c(obs_index,miss_index),]

disease_data_level1_eur <- fread(
  "<PROJECT_DIR>/data/disease_data_level1_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

disease_data_level1_eur <- disease_data_level1_eur[c(obs_index,miss_index),]

disease_data_level0_eur <- fread(
  "<PROJECT_DIR>/data/disease_data_level0_eur.tsv",
  sep       = "\t", header    = TRUE, data.table = FALSE)

disease_data_level0_eur <- disease_data_level0_eur[c(obs_index,miss_index),]


write.table(
  disease_data_level0_eur,
  file = "<PROJECT_DIR>/data/disease_data_level0_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  disease_data_level1_eur,
  file = "<PROJECT_DIR>/data/disease_data_level1_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  disease_data_level2_eur,
  file = "<PROJECT_DIR>/data/disease_data_level2_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  df_allUKB_eur_protein,
  file = "<PROJECT_DIR>/data/df_allUKB_eur_protein.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)

write.table(
  df_extracted_pheno_data528_sub,
  file = "<PROJECT_DIR>/data/extracted_pheno_data528_eur.tsv",
  sep = "\t",        # Use tab as the column separator
  quote = FALSE,     # Do not wrap character or factor data in quotes
  row.names = FALSE, # Do not write out row names
  col.names = TRUE   # Write column headers
)
