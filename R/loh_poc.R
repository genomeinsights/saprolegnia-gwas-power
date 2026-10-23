## loh_poc.R -- proof of concept: recover mitotic LOH tracts from the shallow (~2x)
## first-stage data of the two-stage design.
##
##   Rscript R/loh_poc.R <prefix> <seed> [out.rds]
##
## Needs <prefix>.geno.u8 / .sites.tsv / .inds.tsv / .loh.tsv from
## py/ts_to_geno.py --loh (true LOH tracts from the tree sequence).
##
##  1. 2x reads for all 300 isolates (as in R/analyse_rep.R) and clone assignment
##     with the same single-read IBS clustering.
##  2. Within each inferred clone family (>= 2 isolates), one isolate is deep-
##     sequenced (~15x, as in the design) and its heterozygous calls define the
##     clone's heterozygous sites. A clone-mate that lost heterozygosity over a tract
##     never shows both alleles there; a two-state HMM along each clone-mate's
##     informative sites (depth >= 2) calls tracts. The representative's own tracts
##     are called from sites where clone-mates show both alleles but it is homozygous.
##  3. Shared tracts: LOH calls are binned (BIN bp) into an isolates x bins matrix,
##     treated like a genotype matrix. Adjacent bins are grouped by complete-linkage
##     clustering on r2 between their LOH-state vectors across isolates (as in
##     LDscnR's stage-1 clustering); each block is one shared LOH event, and its
##     carriers are the isolates with LOH over most of the block.
##  4. Scored against the true tracts: recall by tract length, precision, length
##     estimates, and recovery of each shared event's carriers.
## LOH present in every sampled member of a family leaves nothing to compare against
## and is not detectable this way; such tracts are scored separately (not 'testable').
suppressMessages({library(data.table)})

args   <- commandArgs(TRUE)
prefix <- args[1]
seed   <- as.integer(args[2])
outf   <- if (length(args) >= 3) args[3] else paste0(prefix, ".loh_res.rds")
set.seed(seed)

DEPTH_LOW <- 2; DEPTH_HIGH <- 15; ERR <- 0.01; N_IBS <- 50000L
CLONE_CUT <- 0.1    # tighter than the GWAS (0.25): clones founded by a mating between members
                    # of one clone (genetically selfing) must not be merged with the parent clone
RATE_IN   <- 2e-7   # HMM prior: LOH tract starts per bp
LAMBDA    <- 2e5    # HMM prior: mean tract length (bp); not the simulated value
BIN       <- 1e4
R2_BLOCK  <- 0.8
MIN_SITES <- 8      # min informative sites in a called tract
MIN_LLR   <- 15     # min log-likelihood ratio (LOH vs heterozygous) of a called tract
MIN_DENS  <- 0.8    # min density of lost het sites in a tract, relative to the genome-wide
                    # density: sparse 'het' sites are mutations private to the reference
                    # isolate in regions that were already homozygous in the clone's ancestor
MEMBER_QC <- 0.8    # min observed/expected both-allele rate for a family member
CHRLEN    <- 15750000

## ---- read ---------------------------------------------------------------------
inds  <- fread(paste0(prefix, ".inds.tsv"))
sites <- fread(paste0(prefix, ".sites.tsv"))
truth <- fread(paste0(prefix, ".loh.tsv"))
n <- nrow(inds); S <- nrow(sites)
G <- matrix(as.integer(readBin(paste0(prefix, ".geno.u8"), "raw", n * S)), nrow = n, byrow = TRUE)
truth[, iso := match(pedigree_id, inds$pedigree_id)]
truth[, len := right - left]

## ---- 1. 2x reads, clone assignment (as R/analyse_rep.R) --------------------------
d <- matrix(rpois(length(G), DEPTH_LOW), n)
k <- matrix(rbinom(length(G), d, (G / 2) * (1 - ERR) + (1 - G / 2) * ERR), n)
js <- sort(sample.int(S, min(N_IBS, S)))
dj <- d[, js]; kj <- k[, js]
cov <- dj > 0
X <- matrix(0, n, length(js)); X[cov] <- rbinom(sum(cov), 1, kj[cov] / dj[cov])
C <- cov * 1
D <- (X %*% t(C - X) + (C - X) %*% t(X)) / (C %*% t(C))
pd <- ifelse(dj >= 2, 2 * kj * (dj - kj) / (dj * (dj - 1)), NA)
self <- rowMeans(pd, na.rm = TRUE)
Dc <- D - outer(self, self, "+") / 2
Dr <- Dc / median(Dc[upper.tri(Dc)]); diag(Dr) <- 0
fam <- cutree(hclust(as.dist(pmax(Dr, 0)), method = "average"), h = CLONE_CUT)
rm(X, C, D, Dc, dj, kj, cov, pd)
cat(sprintf("[clones] true %d, inferred %d\n", uniqueN(inds$clone_id), uniqueN(fam)))

## ---- 2. per-member HMM within clone families -------------------------------------------
## The two-stage design sequences one isolate per clone deeply (~15x). Its confident
## heterozygous calls are the clone's reference heterozygosity. A clone-mate (2x)
## that has lost heterozygosity over a tract never shows both alleles at those sites.
## The representative's own LOH tracts are found the other way round: sites where
## other clone-mates show both alleles but the representative is homozygous.
viterbi <- function(pos, le) {
  ## states 1 = heterozygous, 2 = LOH; le = per-site log emission probabilities (m x 2)
  m <- length(pos); V <- matrix(0, m, 2); B <- matrix(1L, m, 2)
  V[1, ] <- le[1, ] + log(c(1 - 1e-3, 1e-3))
  for (i in seq_len(m)[-1]) {
    gap <- pos[i] - pos[i - 1]
    p01 <- 1 - exp(-gap * RATE_IN); p10 <- 1 - exp(-gap / LAMBDA)
    if (gap < 0 || gap > CHRLEN) { p01 <- 1e-3; p10 <- 0.5 }      # new chromosome
    s1 <- c(V[i - 1, 1] + log(1 - p01), V[i - 1, 2] + log(p10))
    s2 <- c(V[i - 1, 1] + log(p01),     V[i - 1, 2] + log(1 - p10))
    B[i, 1] <- which.max(s1); B[i, 2] <- which.max(s2)
    V[i, ] <- c(max(s1), max(s2)) + le[i, ]
  }
  st <- integer(m); st[m] <- which.max(V[m, ])
  for (i in (m - 1):1) st[i] <- B[i + 1, st[i + 1]]
  st
}
## run the HMM on one isolate's informative sites and return its LOH tracts
## ref_pos: positions of the reference heterozygous sites, used for each tract's
## density of lost heterozygous sites relative to the genome-wide density
call_tracts <- function(i, f, inf, p_het, p_loh, obs, ref_pos) {
  if (length(inf) < 10) return(NULL)
  le  <- cbind(ifelse(obs, log(p_het), log(1 - p_het)), ifelse(obs, log(p_loh), log(1 - p_loh)))
  st  <- viterbi(sites$gpos[inf], le)
  llr <- le[, 2] - le[, 1]
  r <- rle(st); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  lo <- which(r$values == 2)
  if (!length(lo)) return(NULL)
  g <- sites$gpos[inf]            # tract spans midway to the flanking informative sites
  L <- ifelse(starts[lo] > 1, (g[starts[lo]] + g[pmax(starts[lo] - 1, 1)]) / 2, g[starts[lo]])
  R <- ifelse(ends[lo] < length(g), (g[ends[lo]] + g[pmin(ends[lo] + 1, length(g))]) / 2, g[ends[lo]])
  rel_dens <- (findInterval(R, ref_pos) - findInterval(L, ref_pos)) / (R - L) /
              (length(ref_pos) / (4 * CHRLEN))
  data.table(iso = i, fam = f, left = L, right = R, n_sites = r$lengths[lo],
             llr = vapply(lo, function(z) sum(llr[starts[z]:ends[z]]), 0), rel_dens = rel_dens)
}
calls <- list()
for (f in unique(fam[duplicated(fam)])) {
  mem <- which(fam == f)
  rep_i <- mem[sample.int(length(mem), 1)]
  dh <- rpois(S, DEPTH_HIGH)
  kh <- rbinom(S, dh, (G[rep_i, ] / 2) * (1 - ERR) + (1 - G[rep_i, ] / 2) * ERR)
  rep_het <- kh >= 3 & dh - kh >= 3
  het <- which(rep_het)
  ## clone-mates whose rate of showing both alleles at these sites is far below
  ## expectation are not clone-mates (assignment error) and are dropped
  oth <- setdiff(mem, rep_i)
  ok <- vapply(oth, function(i) {
    inf <- het[d[i, het] >= 2]; dep <- d[i, inf]
    mean(k[i, inf] > 0 & k[i, inf] < dep) / mean(1 - 2 * 0.5^dep)
  }, 0)
  q <- pmin(ok, 1)                  # share of the reference het sites the clone-mate carries
  keep <- ok >= MEMBER_QC; oth <- oth[keep]; q <- q[keep]
  if (!length(oth)) next
  for (j in seq_along(oth)) {
    i <- oth[j]; inf <- het[d[i, het] >= 2]; dep <- d[i, inf]
    calls[[length(calls) + 1]] <- call_tracts(i, f, inf, p_het = q[j] * (1 - 2 * 0.5^dep),
      p_loh = 1 - (1 - ERR)^dep, obs = k[i, inf] > 0 & k[i, inf] < dep, ref_pos = sites$gpos[het])
  }
  ## representative: sites where >= 1 clone-mate shows both alleles (>= 2 minor reads
  ## pooled); observation = representative heterozygous at ~15x
  Ko <- colSums(k[oth, , drop = FALSE]); Do <- colSums(d[oth, , drop = FALSE])
  nbo <- colSums(k[oth, , drop = FALSE] > 0 & k[oth, , drop = FALSE] < d[oth, , drop = FALSE])
  inf <- which(nbo >= 1 & Ko >= 2 & Do - Ko >= 2)
  calls[[length(calls) + 1]] <- call_tracts(rep_i, f, inf, p_het = rep(0.8, length(inf)),
    p_loh = rep(0.02, length(inf)), obs = rep_het[inf], ref_pos = sites$gpos[inf])
}
calls <- rbindlist(calls)
calls <- calls[n_sites >= MIN_SITES & llr >= MIN_LLR & rel_dens >= MIN_DENS]
calls[, len := right - left]
cat(sprintf("[hmm] %d LOH tracts called in %d isolates\n", nrow(calls), uniqueN(calls$iso)))

## ---- 3. shared tracts: LD-style blocks of the isolates x bins LOH matrix ---------------------
nb <- ceiling(4 * CHRLEN / BIN)
Lm <- matrix(0L, n, nb)
for (j in seq_len(nrow(calls))) {
  b <- seq(floor(calls$left[j] / BIN) + 1, floor(calls$right[j] / BIN) + 1)
  Lm[calls$iso[j], b] <- 1L
}
poly <- which(colSums(Lm) > 0)
blocks <- list(); cur <- integer()
r2 <- function(a, b) { if (sd(a) == 0 || sd(b) == 0) return(0); cor(a, b)^2 }
for (b in poly) {
  same_chr <- length(cur) && (cur[1] - 1) %/% (CHRLEN / BIN) == (b - 1) %/% (CHRLEN / BIN)
  if (length(cur) && same_chr && b - cur[length(cur)] == 1 &&
      all(vapply(cur, function(c0) r2(Lm[, c0], Lm[, b]), 0) >= R2_BLOCK)) {
    cur <- c(cur, b)
  } else {
    if (length(cur)) blocks[[length(blocks) + 1]] <- cur
    cur <- b
  }
}
if (length(cur)) blocks[[length(blocks) + 1]] <- cur
ev_hat <- rbindlist(lapply(seq_along(blocks), function(i) {
  bb <- blocks[[i]]; carr <- which(rowMeans(Lm[, bb, drop = FALSE]) >= 0.5)
  data.table(block = i, left = (min(bb) - 1) * BIN, right = max(bb) * BIN,
             carriers = list(carr), n_carriers = length(carr))
}))
cat(sprintf("[blocks] %d LOH blocks, %d shared by >= 2 isolates\n", nrow(ev_hat), sum(ev_hat$n_carriers >= 2)))

## ---- 4. score ------------------------------------------------------------------------------------
## truth tracts are 'segregating' if another sampled member of the same true clone lacks LOH there
truth[, clone := inds$clone_id[iso]]
truth[, mid := (left + right) / 2]
csize <- inds[, .N, by = clone_id]
truth[, fam_n := csize$N[match(clone, csize$clone_id)]]
truth[, n_carry := vapply(seq_len(.N), function(j)
  truth[clone == clone[j] & left <= mid[j] & right >= mid[j], uniqueN(iso)], 1L)]
## ... and the region has heterozygous sites to lose: >= 5 SNPs truly het in some family member
truth[, n_famhet := vapply(seq_len(.N), function(j) {
  s <- which(sites$gpos >= left[j] & sites$gpos < right[j]); m <- which(inds$clone_id == clone[j])
  sum(colSums(G[m, s, drop = FALSE] == 1) > 0) }, 0L)]
truth[, testable := fam_n >= 2 & n_carry < fam_n & n_famhet >= 5]
ov <- function(l1, r1, l2, r2) pmax(0, pmin(r1, r2) - pmax(l1, l2))
truth[, c("called_bp", "called_len") := {
  cc <- calls[iso == .BY$iso]
  if (!nrow(cc)) list(0, NA_real_) else {
    o <- ov(left, right, cc$left, cc$right)
    list(sum(o), if (any(o > 0)) cc$len[which.max(o)] else NA_real_)
  }
}, by = .(iso, left)]
truth[, detected := called_bp >= 0.5 * len]
calls[, true_bp := vapply(seq_len(.N), function(j) {
  tt <- truth[iso == calls$iso[j]]; sum(ov(calls$left[j], calls$right[j], tt$left, tt$right))
}, 0)]
calls[, correct := true_bp >= 0.5 * len]

## true events: tracts in one clone with the same coalescence time that overlap
truth[, event := .GRP, by = .(clone, tmrca)]
truth[, event := paste(event, cumsum(c(TRUE, left[-1] > cummax(right)[-.N]))), by = event]
ev_true <- truth[, .(left = min(left), right = max(right), carriers = list(sort(unique(iso))),
                    n_carriers = uniqueN(iso), fam_n = fam_n[1], testable = all(testable)),
                 by = .(clone, event)]
ev_true[, jaccard := vapply(seq_len(.N), function(j) {
  o <- ov(left[j], right[j], ev_hat$left, ev_hat$right)
  if (!any(o > 0)) return(0)
  h <- ev_hat$carriers[[which.max(o)]]
  h <- h[inds$clone_id[h] == clone[j]]
  t <- carriers[[j]]; length(intersect(h, t)) / length(union(h, t))
}, 0)]

summ <- list(
  n_truth = nrow(truth), n_testable = sum(truth$testable),
  recall = mean(truth[testable == TRUE]$detected), precision = mean(calls$correct),
  n_calls = nrow(calls),
  shared_events = sum(ev_true$testable & ev_true$n_carriers >= 2),
  carriers_exact = mean(ev_true[testable & n_carriers >= 2]$jaccard == 1),
  carriers_jaccard = mean(ev_true[testable & n_carriers >= 2]$jaccard))
str(summ)
saveRDS(list(summary = summ, truth = truth, calls = calls, ev_true = ev_true, ev_hat = ev_hat,
             n_fam_multi = sum(table(fam) >= 2), seed = seed), outf)
cat("[done] ", outf, "\n")
