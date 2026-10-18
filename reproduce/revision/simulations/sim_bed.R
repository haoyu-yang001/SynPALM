## Helpers to write simulated genotypes as PLINK 1 binary files (for BOLT-REML).
## Genotypes are A1 allele counts (0/1/2), SNP-major, no missing values.

## open a .bed connection and write the magic number + SNP-major flag
bed_open <- function(prefix) {
  con <- file(paste0(prefix, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  con
}

## append one SNP (integer vector of A1 counts, length n)
## PLINK codes: 2 copies of A1 -> 00, het -> 10, 0 copies -> 11 (01 = missing),
## four genotypes per byte, first individual in the lowest two bits
bed_write_snp <- function(con, g) {
  code <- c(3L, 2L, 0L)[g + 1L]
  pad  <- (-length(code)) %% 4L
  if (pad) code <- c(code, integer(pad))
  m <- matrix(code, nrow = 4L)
  writeBin(as.raw(m[1, ] + 4L * m[2, ] + 16L * m[3, ] + 64L * m[4, ]), con)
}

bed_write_fam_bim <- function(prefix, n, m, chr = 1L) {
  id <- paste0("i", seq_len(n))
  data.table::fwrite(data.table::data.table(id, id, 0, 0, 0, -9),
                     paste0(prefix, ".fam"), sep = " ", col.names = FALSE)
  data.table::fwrite(data.table::data.table(chr, paste0("s", seq_len(m)), 0, seq_len(m) * 1000L, "A", "G"),
                     paste0(prefix, ".bim"), sep = "\t", col.names = FALSE)
  invisible(id)
}
