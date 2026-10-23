## heatmap_ldcr.R -- does LDscnR's stage-1 LD complexity reduction work on clonal data?
##   Rscript R/heatmap_ldcr.R <prefix> <chr> <from_bp> <to_bp>
## Uses one isolate per true clone (the GWAS input; true genotypes, MAF >= 0.05),
## runs LD decay + stage-1 clustering genome-wide as in R/analyse_rep.R, and draws
## genotype heatmaps of one region: (A) SNPs by position, (B) SNPs clustered on
## 1 - r2, (C) only the stage-1 representative ('central') SNP of each cluster.
suppressMessages({library(data.table); library(SNPRelate); library(LDscnR); library(ComplexHeatmap)})
args <- commandArgs(TRUE)
prefix <- args[1]; chr <- paste0("Chr", args[2]); from <- as.numeric(args[3]); to <- as.numeric(args[4])
set.seed(1)
inds  <- fread(paste0(prefix, ".inds.tsv")); sites <- fread(paste0(prefix, ".sites.tsv"))
n <- nrow(inds); S <- nrow(sites)
G <- matrix(as.integer(readBin(paste0(prefix, ".geno.u8"), "raw", n * S)), nrow = n, byrow = TRUE)
rep_idx <- sort(vapply(split(seq_len(n), inds$clone_id), function(i) i[sample.int(length(i), 1)], 1L))
G <- G[rep_idx, ]
af <- colMeans(G) / 2; keep <- which(pmin(af, 1 - af) >= 0.05)
G <- G[, keep]
map <- data.table(marker = sprintf("s%d", keep), Chr = paste0("Chr", sites$chr[keep]), Pos = sites$pos[keep])
colnames(G) <- map$marker; rownames(G) <- sprintf("iso%d", rep_idx)
cat(sprintf("%d clone representatives x %d SNPs\n", nrow(G), ncol(G)))

gds_path <- tempfile(fileext = ".gds")
gds <- create_gds_from_geno(G, map, gds_path)
LD_decay <- suppressMessages(compute_LD_decay(gds, n_win_decay = 20, slide = 500, seed = 1))
stage1 <- suppressMessages(ld_complexity_reduction(map = map, LD_decay = LD_decay, rho = 0.5, gds = gds))
snpgdsClose(gds); unlink(gds_path)
cl <- as.data.table(stage1$clusters)
print(LD_decay$decay_sum)
nl <- if ("n_loci" %in% names(cl)) cl$n_loci else lengths(cl$members)
cat(sprintf("genome-wide: %d SNPs -> %d stage-1 clusters (%d representatives); cluster size median %g, mean %.1f, max %d; singletons %.0f%%\n",
            nrow(map), nrow(cl), length(unique(na.omit(stage1$pruned))), median(nl), mean(nl), max(nl), 100 * mean(nl == 1)))

ms <- as.data.table(stage1$map_snp)
reg <- ms[Chr == chr & Pos >= from & Pos <= to]
setorder(reg, Pos)
cat(sprintf("region %s:%g-%g: %d SNPs in %d clusters\n", chr, from, to, nrow(reg), uniqueN(reg$CL_id)))
Gr <- G[, reg$marker]
## within-cluster vs between-cluster r2 in the region
R2 <- cor(Gr)^2
same <- outer(reg$CL_id, reg$CL_id, "=="); ut <- upper.tri(R2)
cat(sprintf("region r2: within clusters mean %.2f, between adjacent-ish clusters mean %.2f\n",
            mean(R2[ut & same]), mean(R2[ut & !same])))

row_ord <- rownames(Gr)[hclust(dist(Gr))$order]
deme <- setNames(paste0("site ", inds$deme[rep_idx]), rownames(G))
cl_col <- setNames(as.character(reg$CL_id), reg$marker)
reps <- intersect(unique(na.omit(stage1$pruned)), reg$marker)
reps <- reg$marker[reg$marker %in% reps]
common <- list(row_order = row_ord, row_annotation = deme[row_ord], row_annotation_name = "Site",
               col_annotation_name = "Stage-1 cluster", col_annotation_legend = FALSE,
               width = 9, height = 6)
lab <- sprintf("%s %.1f-%.1f Mb", chr, from / 1e6, to / 1e6)
do.call(plot_genotype_heatmap, c(list(Gr, cluster_columns = FALSE, col_annotation = cl_col,
  title = sprintf("A. %s: %d SNPs by position", lab, ncol(Gr)), out_file = "figures/heatmap_A_position.png"), common))
do.call(plot_genotype_heatmap, c(list(Gr, cluster_columns = TRUE, col_annotation = cl_col,
  title = sprintf("B. Same SNPs clustered on 1 - r2"), out_file = "figures/heatmap_B_clustered.png"), common))
do.call(plot_genotype_heatmap, c(list(Gr[, reps, drop = FALSE], cluster_columns = FALSE, col_annotation = cl_col[reps],
  title = sprintf("C. Stage-1 representative SNPs only: %d of %d", length(reps), ncol(Gr)),
  out_file = "figures/heatmap_C_representatives.png"), common))
