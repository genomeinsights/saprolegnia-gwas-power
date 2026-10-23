## sexrate_rep.R -- LD summaries that may carry the rate of sex, for one simulation.
##   Rscript R/sexrate_rep.R <prefix> <seed> <out.rds>
## One isolate per true clone (clone-corrected, as for I_A), true genotypes, MAF >= 0.05.
## Summaries:
##   b_ldscnr   LDscnR background LD (95% quantile of r2 between chromosomes)
##   a_pred     LDscnR LD-decay rate (per bp), mean over chromosomes
##   c_pred     LDscnR short-range LD
##   r2_between mean r2 between SNPs on different chromosomes, minus the 1/n sampling floor
##   r2_10kb    mean r2 of SNP pairs < 10 kb apart, minus the same floor
##   rbarD      standardised index of association (Agapow & Burt 2001), 1000 random SNPs
suppressMessages({library(data.table); library(SNPRelate); library(LDscnR)})
args <- commandArgs(TRUE); prefix <- args[1]; seed <- as.integer(args[2]); outf <- args[3]
set.seed(seed)
inds  <- fread(paste0(prefix, ".inds.tsv")); sites <- fread(paste0(prefix, ".sites.tsv"))
n <- nrow(inds); S <- nrow(sites)
G <- matrix(as.integer(readBin(paste0(prefix, ".geno.u8"), "raw", n * S)), nrow = n, byrow = TRUE)
rep_idx <- sort(vapply(split(seq_len(n), inds$clone_id), function(i) i[sample.int(length(i), 1)], 1L))
G <- G[rep_idx, , drop = FALSE]; nr <- nrow(G)
af <- colMeans(G) / 2; keep <- which(pmin(af, 1 - af) >= 0.05)
G <- G[, keep]; chr <- sites$chr[keep]; pos <- sites$gpos[keep]
map <- data.table(marker = sprintf("s%d", keep), Chr = paste0("Chr", chr), Pos = sites$pos[keep])
colnames(G) <- map$marker

## direct r2 summaries
pairs_between <- function(k) { i <- sample(ncol(G), k, TRUE); j <- sample(ncol(G), k, TRUE); ok <- chr[i] != chr[j]; cbind(i[ok], j[ok]) }
r2p <- function(ij) vapply(seq_len(nrow(ij)), function(z) suppressWarnings(cor(G[, ij[z, 1]], G[, ij[z, 2]]))^2, 0)
floor_r2 <- 1 / (nr - 1)
r2_between <- mean(r2p(pairs_between(5000)), na.rm = TRUE) - floor_r2
i <- sample(ncol(G) - 1, 20000, TRUE)
j <- pmin(i + sample(1:20, 20000, TRUE), ncol(G))
ok <- chr[i] == chr[j] & pos[j] - pos[i] < 1e4 & i != j
r2_10kb <- mean(r2p(cbind(i[ok], j[ok])), na.rm = TRUE) - floor_r2

## rbarD on 1000 random SNPs
loc <- sample(ncol(G), min(1000, ncol(G)))
pr <- which(upper.tri(matrix(0, nr, nr)), arr.ind = TRUE)
dj <- abs(G[pr[, 1], loc] - G[pr[, 2], loc])            # pairs x loci distances
vj <- apply(dj, 2, var); VO <- var(rowSums(dj))
rbarD <- (VO - sum(vj)) / (sum(sqrt(vj))^2 - sum(vj))

## LDscnR decay
gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(G, map, gds_path)
LD <- tryCatch(suppressWarnings(suppressMessages(compute_LD_decay(gds, n_win_decay = 20, slide = 500, seed = seed))),
               error = function(e) NULL)
snpgdsClose(gds); unlink(gds_path)
ds <- if (!is.null(LD)) as.data.table(LD$decay_sum) else NULL
out <- data.table(n_clones = nr, n_snp = ncol(G), het = mean(G == 1),
                  r2_between = r2_between, r2_10kb = r2_10kb, rbarD = rbarD,
                  b_ldscnr = if (!is.null(LD)) LD$by_chr[[1]]$b else NA_real_,
                  a_pred = if (!is.null(ds)) mean(ds$a_pred) else NA_real_,
                  c_pred = if (!is.null(ds)) mean(ds$c_pred) else NA_real_,
                  n_chr_fit = if (!is.null(ds)) nrow(ds) else 0L)
print(out)
saveRDS(out, outf)
