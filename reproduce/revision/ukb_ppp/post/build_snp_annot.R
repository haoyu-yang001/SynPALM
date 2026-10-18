## One table of the 591,558 QC'd array SNPs: chr, rsid, GRCh37 position and
## alleles from the .bim, GRCh38 position and REF/ALT from the Olink rsid map.
## BEDMatrix counts copies of the .bim A1 allele, so A1 is the effect allele
## of every beta in the scan output (REF/ALT are carried only for lookup).
library(data.table)
setwd("<PROJECT_DIR>")
bim <- fread("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bim", header = FALSE,
             col.names = c("chr", "rsid", "cm", "POS37", "A1", "A2"))[, cm := NULL]
maps <- list.files("<UKB_DIR>/gwas_data/olink_rsid_map/", "\\.gz$", full.names = TRUE)
mp <- rbindlist(lapply(maps, function(f) fread(f, select = c("rsid", "POS38", "REF", "ALT"))))
mp <- mp[rsid %in% bim$rsid]
mp <- unique(mp, by = "rsid")
ann <- merge(bim, mp, by = "rsid", all.x = TRUE, sort = FALSE)
ann <- ann[match(bim$rsid, ann$rsid)]
stopifnot(identical(ann$rsid, bim$rsid))
ann[, EA_is := fifelse(is.na(REF), NA_character_, fifelse(A1 == ALT, "ALT", fifelse(A1 == REF, "REF", "other")))]
cat("SNPs:", nrow(ann), " with POS38:", sum(!is.na(ann$POS38)), "\n"); print(table(ann$EA_is, useNA = "ifany"))
saveRDS(ann, "data/snp_annot_v2geno.rds")
