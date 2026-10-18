## Variance-component estimation: method of moments (Syn-PALM) vs REML (GMMAT)
## vs BOLT-REML (Supp Figs S12, S13; replaces variance_estimate_compare*.R, whose
## BOLT-REML results had no generating code).
## BOLT-REML needs genotypes, so the genetic components are generated from
## simulated SNPs: per array, M = 10,000 SNPs (MAF ~ U(0.05, 0.5)) are dropped
## through sib pairs, Z = standardised genotypes, and for replicate k
##     (u_T, u_S) = Z W_k / sqrt(M) %*% chol(Gmat),  W_k ~ N(0, I_{M x 2}),
## so Cov(u) = Gmat (x) ZZ'/M. Residuals e ~ N(0, Emat) iid; Y = u_T + e_T,
## S = u_S + e_S, covariate x ~ N(0, 1) with no effect (as before). The first
## n_obs = 30,000 individuals are labelled. Gmat / Emat: tauT = 0.7, tauS = 0.4,
## sigmaT = 0.7, sigmaS = 0.5 (SDs) and the covariances of sim_cov_TS(rho).
## Method of moments and GMMAT use the sparse GRM of the sib pairs (empirical
## kinship ZZ'/M within pairs, 0 elsewhere), as Syn-PALM does in UKB; BOLT-REML
## uses all M SNPs. All three fit, as in the submitted figure, the variances of
## Y on the labelled and of S on all individuals, and the covariances on the
## labelled individuals (bivariate fit). MoM estimates are the raw (unconstrained)
## moment solutions. Times are wall-clock seconds of the estimation only, one thread.
##   Rscript variance_components_rev1.R <miss> <array> [n_rep=10] [M=10000] [rep_first=1] [rep_last=n_rep]
## Each (replicate, rho) fit is saved on its own in <OUT>/parts/ as soon as it is done,
## and parts that exist are skipped, so a preempted or timed-out task resumes; a task
## can also be restricted to replicates rep_first..rep_last (the RNG state of every
## (replicate, rho) is recorded in the original order, so the results are identical to a
## single run, also when VC_RHO_FIRST changes the order of the fits). When all
## n_rep x 4 parts of an array exist they are combined into the array file.
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("simulation_code/rev1/sim_bed.R"); library(GMMAT) })
options(synsurrg.verbose = FALSE)
args     <- commandArgs(trailingOnly = TRUE)
miss     <- as.numeric(args[1])
array_id <- as.integer(args[2])
n_rep    <- if (length(args) >= 3) as.integer(args[3]) else 10L
M        <- if (length(args) >= 4) as.integer(args[4]) else 10000L
rep_first <- if (length(args) >= 5) as.integer(args[5]) else 1L
rep_last  <- if (length(args) >= 6) as.integer(args[6]) else n_rep
BOLT     <- path.expand("~/bin/bolt")
Sys.setenv(LD_LIBRARY_PATH = paste("<INTEL_ONEAPI_DIR>/compiler/2024.2/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))

rho_vec <- c(0, 0.25, 0.5, 0.75)
n_obs <- 3e4; n <- n_obs + round(n_obs * miss / (1 - miss))
tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5
OUT <- "simulation_code/results_rev1/variance_components/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("%smissing_%s_array_%d.rds", OUT, miss, array_id)
if (file.exists(out_file)) quit(save = "no")
PART <- paste0(OUT, "parts/"); dir.create(PART, showWarnings = FALSE)
part_file <- function(k, rho) sprintf("%smissing_%s_array_%d_rep_%d_rho_%s.rds", PART, miss, array_id, k, rho)
## env VC_RHO_ONLY (e.g. "0.75"): fit only these rho in this run (the others are left for a later run)
rho_only <- suppressWarnings(as.numeric(strsplit(Sys.getenv("VC_RHO_ONLY", ""), ",")[[1]]))
rho_only <- rho_only[!is.na(rho_only)]
rho_run  <- if (length(rho_only)) intersect(rho_vec, rho_only) else rho_vec
todo_any <- !all(file.exists(mapply(part_file, rep(rep_first:rep_last, each = length(rho_run)), rho_run)))
if (!todo_any) {          # nothing to fit in this range: combine if the whole array is done, then stop
  pf <- mapply(part_file, rep(seq_len(n_rep), each = 4), rho_vec)
  if (all(file.exists(pf))) { saveRDS(do.call(rbind, lapply(pf, readRDS)), out_file); invisible(file.remove(pf)) }
  quit(save = "no")
}
seed <- 20261014L + 1000L * match(miss, c(0, 0.25, 0.5, 0.75, 0.9)) + array_id
set.seed(seed)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

## ---- genotypes (streamed to a .bed), genetic effects and sparse GRM --------
tmp <- tempfile("vc_", tmpdir = Sys.getenv("TMPDIR", tempdir())); dir.create(tmp)
prefix <- file.path(tmp, "g")
con <- bed_open(prefix)
U  <- matrix(0, n, 2 * n_rep)               # Z W / sqrt(M) for all replicates
kd <- numeric(n); ko <- numeric(n / 2)      # diagonal and within-pair kinship
odd <- seq(1, n, 2); evn <- odd + 1
drop1 <- function(g) rbinom(length(g), 1L, g / 2)
B <- 100L                                   # SNPs per block for U += Z W
for (b0 in seq(1L, M, by = B)) {
  nb <- min(B, M - b0 + 1L); Zb <- matrix(0, n, nb)
  for (j in seq_len(nb)) {
    p  <- runif(1, 0.05, 0.5)
    gf <- rbinom(n / 2, 2L, p); gm <- rbinom(n / 2, 2L, p)
    g  <- integer(n); g[odd] <- drop1(gf) + drop1(gm); g[evn] <- drop1(gf) + drop1(gm)
    bed_write_snp(con, g)
    z  <- (g - 2 * p) / sqrt(2 * p * (1 - p))
    Zb[, j] <- z; kd <- kd + z^2; ko <- ko + z[odd] * z[evn]
  }
  U <- U + Zb %*% matrix(rnorm(nb * 2 * n_rep), nb)
}
rm(Zb)
close(con); bed_write_fam_bim(prefix, n, M)
U <- U / sqrt(M); kd <- kd / M; ko <- ko / M
GRM <- sparseMatrix(i = c(seq_len(n), odd, evn), j = c(seq_len(n), evn, odd),
                    x = c(kd, ko, ko), dims = c(n, n))
colnames(GRM) <- rownames(GRM) <- paste0("i", seq_len(n))
obs <- seq_len(n_obs)
GRM_oo <- GRM[obs, obs]
tick("genotypes: n = ", n, ", M = ", M, "; mean within-pair kinship ", round(mean(ko), 4))

## ---- estimators -------------------------------------------------------------
mom <- function(y, s, x) {
  rT <- residuals(lm(y[obs] ~ x[obs])); rS <- residuals(lm(s ~ x))
  Go <- GRM_oo; diag(Go) <- 0; Ga <- GRM; diag(Ga) <- 0
  tT <- as.numeric(t(rT) %*% Go %*% rT) / sum(Go^2)
  tS <- as.numeric(t(rS) %*% Ga %*% rS) / sum(Ga^2)
  tX <- as.numeric(t(rS[obs]) %*% Go %*% rT) / sum(Go^2)
  c(tau_T2 = tT, tau_S2 = tS, sigma_T2 = (sum(rT^2) - tT * sum(diag(GRM_oo))) / n_obs,
    sigma_S2 = (sum(rS^2) - tS * sum(diag(GRM))) / n, tau_TS = tX,
    sigma_TS = (sum(rS[obs] * rT) - tX * sum(diag(GRM_oo))) / n_obs)
}
reml <- function(y, s, x) {
  d  <- data.frame(id = rownames(GRM), y = y, s = s, x = x)
  f2 <- glmmkin(cbind(y, s) ~ x, data = d[obs, ], kins = GRM_oo, id = "id", family = gaussian())
  fT <- glmmkin(y ~ x, data = d[obs, ], kins = GRM_oo, id = "id", family = gaussian())
  fS <- glmmkin(s ~ x, data = d, kins = GRM, id = "id", family = gaussian())
  c(tau_T2 = fT$theta[[2]], tau_S2 = fS$theta[[2]], sigma_T2 = fT$theta[[1]], sigma_S2 = fS$theta[[1]],
    tau_TS = f2$theta$kins1[1, 2], sigma_TS = f2$theta$residuals[1, 2])
}
## BOLT drops individuals with a missing phenotype itself (and builds its GRM
## from the individuals analysed), so fits on Y use the labelled individuals.
bolt_run <- function(cols) {
  log <- file.path(tmp, "bolt.log")
  args <- c(paste0("--bfile=", prefix), paste0("--phenoFile=", file.path(tmp, "ph.txt")),
            paste0("--phenoCol=", cols), paste0("--covarFile=", file.path(tmp, "ph.txt")),
            "--qCovarCol=x", "--reml", "--numThreads=1")
  st <- system2(BOLT, args, stdout = log, stderr = log)
  L <- readLines(log)
  if (st != 0 || !any(grepl("^Total elapsed time", L))) stop("BOLT failed:\n", paste(tail(L, 20), collapse = "\n"))
  v  <- as.numeric(sub(".*sigma2: ([-0-9.e+]+) .*", "\\1", grep("variance sigma2:", L, value = TRUE)))
  vc <- grep("^Variance component [01]: \\[", L, value = TRUE)
  vc <- tail(vc, 2)                                          # final estimates
  num <- lapply(vc, function(s) {
    a <- as.numeric(strsplit(gsub("[^-0-9.e+,]", "", sub(".*\\]", "", s)), ",")[[1]])
    if (length(a) == 4) a[c(1, 2, 4)] else a })                 # 2x2 -> (T, TS, S)
  ## BOLT reports components of the standardised phenotypes
  sc <- if (length(v) == 2) c(v[1], sqrt(v[1] * v[2]), v[2]) else v[1]
  list(env = num[[1]] * sc, gen = num[[2]] * sc)                 # (T, TS, S) or (T)
}
bolt <- function(y, s, x) {
  data.table::fwrite(data.table::data.table(FID = rownames(GRM), IID = rownames(GRM), y = y, s = s, x = x),
                     file.path(tmp, "ph.txt"), sep = " ", na = "NA", quote = FALSE)
  b2 <- bolt_run(c("y", "s"))             # labelled: T, TS, S
  bS <- bolt_run("s")                     # all: S
  bT <- bolt_run("y")                     # labelled: T
  c(tau_T2 = bT$gen[1], tau_S2 = bS$gen[1], sigma_T2 = bT$env[1], sigma_S2 = bS$env[1],
    tau_TS = b2$gen[2], sigma_TS = b2$env[2])
}

## ---- replicates ---------------------------------------------------------------
## The noise of replicate k and rho is drawn in the fixed order k = 1..n_rep, rho = rho_vec;
## the RNG state at the start of every (k, rho) is recorded first, so that the fits can be
## run in another order (env VC_RHO_FIRST, e.g. "0.75": that rho first for all replicates)
## with exactly the same draws (the fits themselves do not use the R RNG).
rho_first <- suppressWarnings(as.numeric(strsplit(Sys.getenv("VC_RHO_FIRST", ""), ",")[[1]]))
keys   <- expand.grid(rho = rho_vec, k = seq_len(n_rep))      # rho varies fastest, as before
states <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  states[[i]] <- get(".Random.seed", envir = globalenv())
  invisible(rnorm(2 * n)); invisible(rnorm(n))                # e and x below
}
for (i in order(!(keys$rho %in% rho_first), seq_len(nrow(keys)))) {
  k <- keys$k[i]; rho <- keys$rho[i]
  if (k < rep_first || k > rep_last || !(rho %in% rho_run) || file.exists(part_file(k, rho))) next
  assign(".Random.seed", states[[i]], envir = globalenv())
  cv <- sim_cov_TS(rho, tauT, tauS, sigmaT, sigmaS)
  Gm <- matrix(c(tauT^2, cv[["tauTS"]], cv[["tauTS"]], tauS^2), 2)
  Em <- matrix(c(sigmaT^2, cv[["sigmaTS"]], cv[["sigmaTS"]], sigmaS^2), 2)
  u  <- U[, 2 * k - 1:0] %*% chol(Gm)
  e  <- matrix(rnorm(2 * n), n) %*% chol(Em)
  x  <- rnorm(n)
  y  <- u[, 1] + e[, 1]; y[-obs] <- NA
  s  <- u[, 2] + e[, 2]
  truth <- c(tau_T2 = tauT^2, tau_S2 = tauS^2, sigma_T2 = sigmaT^2, sigma_S2 = sigmaS^2,
             tau_TS = cv[["tauTS"]], sigma_TS = cv[["sigmaTS"]])
  rows <- list()
  for (m in c("MoM", "REML", "BOLT-REML")) {
    f  <- switch(m, MoM = mom, REML = reml, `BOLT-REML` = bolt)
    tm <- system.time(est <- tryCatch(f(y, s, x), error = function(e) { message(m, ": ", conditionMessage(e)); NULL }))[["elapsed"]]
    if (is.null(est)) next
    rows[[length(rows) + 1]] <- data.frame(miss = miss, rho = rho, array = array_id, rep = k, method = m,
                                           param = names(truth), truth = truth, estimate = est[names(truth)],
                                           seconds = tm, n = n, n_obs = n_obs, M = M, seed = seed)
  }
  ## a part is saved only if all three methods succeeded (e.g. BOLT killed on a bad node -> redone on rerun)
  if (length(rows) == 3L) saveRDS(do.call(rbind, rows), part_file(k, rho)) else message("rep ", k, " rho ", rho, ": a fit failed, not saved")
  tick("rep ", k, " rho ", rho, " done")
}
unlink(tmp, recursive = TRUE)
pf <- mapply(part_file, rep(seq_len(n_rep), each = 4), rho_vec)
if (all(file.exists(pf))) {
  saveRDS(do.call(rbind, lapply(pf, readRDS)), out_file)
  invisible(file.remove(pf))
  tick("saved ", out_file)
} else tick(sum(file.exists(pf)), " of ", length(pf), " parts of array ", array_id, " present; ", out_file, " not yet written")
