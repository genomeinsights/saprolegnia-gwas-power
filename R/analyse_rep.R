## analyse_rep.R -- one replicate of the two-stage clone-aware GWAS design.
##
##   Rscript R/analyse_rep.R <prefix> <seed> [out.rds]
##
## <prefix>.geno.u8 / .sites.tsv / .inds.tsv come from py/ts_to_geno.py. They are
## intermediate files: run_one.sh writes them to out/tmp/, deletes them once this
## script has finished, and they are not in the repository. Regenerate them for any
## replicate with ./run_one.sh C REP (seeded), keeping them by commenting out the
## final rm -f lines there. Only the output .rds (out/res/) is kept.
##
## The Simes unit test (ldscnr_simes) is computed below but not reported in the
## figure or summaries, which use emmax_snp and ldscnr_consensus only.
##
##  1. Shallow sequencing of all isolates (depth ~ Pois(DEPTH_LOW), error ERR):
##     single-read IBS distances (as ANGSD -doIBS 1), corrected for each
##     isolate's own read-pair mismatch rate, so clone-mates sit near 0.
##  2. Clones = average-linkage clusters below CLONE_CUT of the median
##     unrelated distance. One random isolate per inferred clone.
##  3. Deep sequencing of the representatives (depth ~ Pois(DEPTH_HIGH)):
##     genotype likelihoods -> called genotypes used for all scans. With
##     DEPTH_HIGH < 8 (shallow test): posterior mean dosages under HWE (allele
##     frequencies by EM) for EMMAX and consensus units; rounded dosages for the
##     GDS (LD decay, stage-1 clustering, GRM).
##  4. Phenotypes (one value per clone): for each PVE in PVE_MAJOR and each
##     chromosome, one major QTL + polygenic background + deme effect + noise.
##  5. Scans, reusing the LDscnR simulation pipeline settings
##     (module_sim_3sp53: rho 0.5, size floor 2, GCTA GRM on stage-1 pruned
##     markers, BH alpha 0.05, stage-2 region assembly):
##       emmax_snp        single-SNP EMMAX, BH, hits merged within GAP bp
##       ldscnr_consensus EMMAX on stage-1 consensus dosages (statistic unit)
##       ldscnr_simes     Simes-combined EMMAX SNP p-values per stage-1 unit
##  6. Scoring: see section 5 (LD-aware truth around the major QTL).
suppressMessages({library(data.table); library(SNPRelate); library(LDscnR)})

args   <- commandArgs(TRUE)
prefix <- args[1]
seed   <- as.integer(args[2])
outf   <- if (length(args) >= 3) args[3] else paste0(prefix, ".res.rds")
set.seed(seed)

DEPTH_LOW  <- as.numeric(Sys.getenv("DEPTH_LOW", 2))
DEPTH_HIGH <- as.numeric(Sys.getenv("DEPTH_HIGH", 15))
ERR        <- 0.01
CLONE_CUT  <- as.numeric(Sys.getenv("CLONE_CUT", 0.25))
N_IBS      <- 50000L
MAF_MIN    <- 0.05
PVE_MAJOR  <- c(0.05, 0.10, 0.20, 0.30)
PVE_BG     <- 0.10; N_BG <- 100L
PVE_DEME   <- 0.05
ALPHA      <- 0.05
GAP        <- 3e5     # ld_outlier_test()'s default physical gap
WIN        <- 1e5     # truth window (= DISTANCE_THRESHOLD in the sim pipeline)

## ---- read ---------------------------------------------------------------------
inds  <- fread(paste0(prefix, ".inds.tsv"))
sites <- fread(paste0(prefix, ".sites.tsv"))
n <- nrow(inds); S <- nrow(sites)
G <- matrix(readBin(paste0(prefix, ".geno.u8"), "raw", n * S), nrow = n, byrow = TRUE)
G <- matrix(as.integer(G), nrow = n)
cat(sprintf("[read] %d isolates x %d sites ; true distinct clones %d\n", n, S, uniqueN(inds$clone_id)))

## ---- reads ----------------------------------------------------------------------
sim_reads <- function(g, depth) {
  d <- matrix(rpois(length(g), depth), nrow(g))
  p <- (g / 2) * (1 - ERR) + (1 - g / 2) * ERR
  k <- matrix(rbinom(length(g), d, p), nrow(g))
  list(d = d, k = k)
}

## ---- 1. shallow stage: clone assignment --------------------------------------
js <- sort(sample.int(S, min(N_IBS, S)))
rl <- sim_reads(G[, js], DEPTH_LOW)
cov <- rl$d > 0
X <- matrix(0, n, length(js)); X[cov] <- rbinom(sum(cov), 1, rl$k[cov] / rl$d[cov])
C <- cov * 1
mism <- X %*% t(C - X) + (C - X) %*% t(X)
nn   <- C %*% t(C)
D <- mism / nn
## self mismatch: two reads (without replacement) at sites with depth >= 2
two <- rl$d >= 2
pd  <- ifelse(two, 2 * rl$k * (rl$d - rl$k) / (rl$d * (rl$d - 1)), NA)
self <- rowMeans(pd, na.rm = TRUE)
Dc <- D - outer(self, self, "+") / 2
Dr <- Dc / median(Dc[upper.tri(Dc)])
diag(Dr) <- 0
hc <- hclust(as.dist(pmax(Dr, 0)), method = "average")
est_clone <- cutree(hc, h = CLONE_CUT)
same_true <- outer(inds$clone_id, inds$clone_id, "==")[upper.tri(Dr)]
same_est  <- outer(est_clone, est_clone, "==")[upper.tri(Dr)]
ari <- {
  tab <- table(inds$clone_id, est_clone); ch2 <- function(x) sum(x * (x - 1) / 2)
  a <- ch2(tab); b <- ch2(rowSums(tab)); c <- ch2(colSums(tab)); t <- n * (n - 1) / 2
  (a - b * c / t) / ((b + c) / 2 - b * c / t)
}
clone_diag <- list(
  n_true = uniqueN(inds$clone_id), n_est = uniqueN(est_clone), ari = ari,
  pair_sens = mean(same_est[same_true]), pair_fpr = mean(same_est[!same_true]),
  d_same = if (any(same_true)) quantile(Dr[upper.tri(Dr)][same_true], c(.5, .99)) else NA,
  d_diff = quantile(Dr[upper.tri(Dr)][!same_true], c(.001, .01, .5)))
cat(sprintf("[clones] true %d, inferred %d, ARI %.3f, pair sens %.3f, pair FPR %.2e\n",
            clone_diag$n_true, clone_diag$n_est, ari, clone_diag$pair_sens, clone_diag$pair_fpr))

## ---- 2. one isolate per inferred clone ---------------------------------------------
rep_idx <- sort(vapply(split(seq_len(n), est_clone), function(i) i[sample.int(length(i), 1)], 1L))
nr <- length(rep_idx)
Gt <- G[rep_idx, , drop = FALSE]                              # true genotypes

## ---- 3. deep stage: GL -> genotype calls ------------------------------------------------
rh <- sim_reads(Gt, DEPTH_HIGH)
ll <- sapply(0:2, function(g) {
  p <- (g / 2) * (1 - ERR) + (1 - g / 2) * ERR
  dbinom(rh$k, rh$d, p, log = TRUE)
})
dim(ll) <- c(nr, S, 3)
if (DEPTH_HIGH >= 8) {
  ## deep: maximum-likelihood genotype calls
  Gc <- matrix(max.col(matrix(ll, ncol = 3), ties.method = "random") - 1L, nr, S)
  Gc[rh$d == 0] <- NA
  af <- colMeans(Gc, na.rm = TRUE) / 2
  Gd <- Gc
} else {
  ## shallow: posterior mean dosages under HWE, allele frequencies by EM
  L0 <- exp(ll[, , 1]); L1 <- exp(ll[, , 2]); L2 <- exp(ll[, , 3])
  af <- pmin(pmax(colSums(rh$k) / pmax(colSums(rh$d), 1), 0.01), 0.99)
  for (it in 1:20) {
    w0 <- sweep(L0, 2, (1 - af)^2, "*"); w1 <- sweep(L1, 2, 2 * af * (1 - af), "*")
    w2 <- sweep(L2, 2, af^2, "*"); tot <- w0 + w1 + w2
    Gd <- (w1 + 2 * w2) / tot
    af <- pmin(pmax(colMeans(Gd) / 2, 1e-4), 1 - 1e-4)
  }
  rm(L0, L1, L2, w0, w1, w2, tot)
  Gc <- matrix(as.integer(round(Gd)), nr, S)    # integer calls, for GDS (LD, clustering, GRM)
}
call_err <- mean(Gc != Gt, na.rm = TRUE)
maf <- pmin(af, 1 - af)
keep <- which(maf >= MAF_MIN & colSums(is.na(Gc)) == 0)
Gc <- Gc[, keep]; Gd <- Gd[, keep]; Gt <- Gt[, keep]
map <- data.table(marker = sprintf("s%d", keep), Chr = paste0("Chr", sites$chr[keep]),
                  Pos = sites$pos[keep], gpos = sites$gpos[keep])
colnames(Gc) <- colnames(Gd) <- map$marker
cat(sprintf("[deep] %d representatives, %d SNPs (MAF >= %.2f), call error %.2e\n",
            nr, nrow(map), MAF_MIN, call_err))

## ---- 4. genotype-only LDscnR set-up (shared by all phenotypes) --------------------------
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(Gc, map, gds_path)
on.exit({ try(snpgdsClose(gds), silent = TRUE); unlink(gds_path) }, add = TRUE)
LD_decay <- suppressMessages(compute_LD_decay(gds, n_win_decay = 20, slide = 500, seed = seed))
stage1 <- suppressMessages(ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = 0.5, gds = gds))
grm_markers <- unique(na.omit(stage1$pruned))
GRM <- snpgdsGRM(gds, snp.id = grm_markers, method = "GCTA", verbose = FALSE,
                 autosome.only = FALSE, missing.rate = 1)$grm
Pm <- emmax_setup(Gd, GRM)
um <- ld_unit_matrix(Gd, stage1, map, size_floor = 2L, repr = "consensus_dosage")
Pu <- emmax_setup(um, GRM)
common <- list(alpha = ALPHA, assembly = "stage2_discovered", GTs = Gc, LD_decay = LD_decay,
               score_threshold = 0.80, distance_threshold = 1e5, size_floor = 2L)
merge_phys <- function(r, gap) {
  if (!nrow(r)) return(data.table(Chr = character(), from = numeric(), to = numeric()))
  setorder(r, Chr, from)
  r[, grp := cumsum(c(TRUE, diff(from) > gap | Chr[-1] != Chr[-.N]))]
  r[, .(from = min(from), to = max(to)), by = .(Chr, grp)][, grp := NULL][]
}

## ---- 5. phenotypes: one major QTL each, on true genotypes ----------------------------------
## For each PVE and each chromosome, a QTL at a random SNP with true MAF >= 0.1,
## plus polygenic background (N_BG SNPs, PVE_BG), deme effect (PVE_DEME), noise.
## A region is a true discovery if it is on the QTL chromosome and either lies
## within WIN bp of the QTL or contains a SNP with r2 >= R2_TRUE to it (true
## genotypes of the representatives). The QTL is detected if any region is true.
R2_TRUE <- 0.2
tmaf <- pmin(colMeans(Gt), 2 - colMeans(Gt)) / 2
scale_to <- function(x, v) { x <- x - mean(x); if (var(x) == 0) x else x * sqrt(v / var(x)) }
deme <- inds$deme[rep_idx]
chrs <- sort(unique(map$Chr))
is_true <- function(r, q) {
  if (!nrow(r)) return(logical())
  vapply(seq_len(nrow(r)), function(i) {
    if (r$Chr[i] != map$Chr[q]) return(FALSE)
    if (map$Pos[q] >= r$from[i] - WIN && map$Pos[q] <= r$to[i] + WIN) return(TRUE)
    idx <- which(map$Chr == r$Chr[i] & map$Pos >= r$from[i] & map$Pos <= r$to[i])
    any(suppressWarnings(cor(Gt[, idx, drop = FALSE], Gt[, q]))^2 >= R2_TRUE, na.rm = TRUE)
  }, TRUE)
}
res <- list(); regs <- list(); k <- 0L
for (pve in PVE_MAJOR) for (ch in chrs) {
  cand <- which(map$Chr == ch & tmaf >= 0.1)
  q  <- cand[sample.int(length(cand), 1)]
  bg <- sample(setdiff(seq_len(nrow(map)), q), N_BG)
  y <- scale_to(Gt[, q], pve) +
       scale_to(as.numeric(Gt[, bg] %*% rnorm(N_BG)), PVE_BG) +
       scale_to(rnorm(uniqueN(deme))[match(deme, sort(unique(deme)))], PVE_DEME)
  y <- y + scale_to(rnorm(nr), 1 - pve - PVE_BG - PVE_DEME)

  p_snp <- emmax_fast(Pm, y); names(p_snp) <- map$marker
  regions <- list(emmax_snp = merge_phys(map[p.adjust(p_snp, "BH") <= ALPHA, .(Chr, from = Pos, to = Pos)], GAP))
  pu <- emmax_fast(Pu, y)
  regions$ldscnr_consensus <- suppressMessages(do.call(ld_outlier_test,
    c(list(stage1, map, pu, statistic = "unit"), common)))$regions[, .(Chr, from, to)]
  regions$ldscnr_simes <- suppressMessages(do.call(ld_outlier_test,
    c(list(stage1, map, p_snp, statistic = "simes"), common)))$regions[, .(Chr, from, to)]
  for (m in names(regions)) {
    tr <- is_true(regions[[m]], q)
    k <- k + 1L
    res[[k]] <- data.table(method = m, pve = pve, qtl_chr = ch, qtl_pos = map$Pos[q], qtl_maf = tmaf[q],
                           detected = any(tr), n_regions = length(tr), n_false = sum(!tr))
    if (length(tr)) regs[[k]] <- regions[[m]][, `:=`(method = m, pve = pve, qtl_chr = ch,
                                                    qtl_pos = map$Pos[q], true = tr)]
  }
}
res <- rbindlist(res)
print(res[, .(power = mean(detected), regions = sum(n_regions), false = sum(n_false)), by = .(method, pve)])

saveRDS(list(res = res, regions = rbindlist(regs), clone = clone_diag, n_isolates = n, n_rep = nr, n_snp = nrow(map),
             n_grm = length(grm_markers), n_units = nrow(stage1$clusters), call_err = call_err,
             seed = seed), outf)
cat("[done] ", outf, "\n")
