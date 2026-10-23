## depth_plot.R -- GWAS power and false loci when the clone representatives are
## sequenced at 2x, 4x or 15x (out/res_d2, out/res_d4, out/res). Run from the repo root.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
dirs <- c("2" = "out/res_d2", "4" = "out/res_d4", "15" = "out/res")
x <- rbindlist(lapply(names(dirs), function(dp) {
  f <- list.files(dirs[[dp]], "\\.res\\.rds$", full.names = TRUE)
  rbindlist(lapply(f, function(ff) {
    r <- readRDS(ff); C <- as.numeric(sub("^C([0-9.]+)_r.*", "\\1", basename(ff)))
    r$res[, `:=`(depth = dp, C = C, n_est = r$clone$n_est, call_err = r$call_err)]
  }))
}))
x <- x[C <= 0.995 & method %in% c("emmax_snp", "ldscnr_consensus")]
x[, method := factor(c(emmax_snp = "EMMAX, single SNP", ldscnr_consensus = "LDscnR, LD units")[method],
                     levels = c("EMMAX, single SNP", "LDscnR, LD units"))]
x[, depth := factor(paste0(depth, "x"), levels = c("2x", "4x", "15x"))]
err <- unique(x[, .(depth, C, call_err)])[, .(genotype_error = mean(call_err)), by = depth]
pw <- x[, .(n_est = mean(n_est), power = mean(detected), n = .N), by = .(depth, C, method, pve)]
fd <- x[, .(n_est = mean(n_est), fdp = sum(n_false) / pmax(sum(n_regions), 1)), by = .(depth, C, method)]
fwrite(pw, "out/summary_depth_power.csv"); fwrite(fd, "out/summary_depth_fdr.csv")
print(err); print(dcast(pw[pve %in% c(0.1, 0.2, 0.3)], method + pve + depth ~ C, value.var = "power"))
print(dcast(fd, method + depth ~ C, value.var = "fdp"))

cols <- c("2x" = "#9ec5f4", "4x" = "#5c9be3", "15x" = "#174f94")
th <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
xs <- scale_x_continuous("Distinct clones among 300 isolates", breaks = round(sort(unique(pw[depth == "15x" & method == method[1] & pve == 0.2]$n_est))))
pA <- ggplot(pw[pve %in% c(0.2, 0.3)], aes(n_est, power, colour = depth, linetype = method, shape = method,
                                           group = interaction(depth, method))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_wrap(~ paste0("QTL explains ", pve * 100, "% of variance")) + xs +
  scale_y_continuous("Power (QTL detected)", limits = c(0, 1)) +
  scale_colour_manual("Depth of clone\nrepresentatives", values = cols) +
  scale_linetype_manual("Method", values = c("dashed", "solid")) +
  scale_shape_manual("Method", values = c(1, 16)) + th
pB <- ggplot(fd, aes(n_est, fdp, colour = depth, linetype = method, shape = method, group = interaction(depth, method))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + xs +
  scale_y_continuous("False discovery proportion", limits = c(0, NA)) +
  scale_colour_manual("Depth of clone\nrepresentatives", values = cols) +
  scale_linetype_manual("Method", values = c("dashed", "solid")) +
  scale_shape_manual("Method", values = c(1, 16)) + th
fig <- pA + pB + plot_layout(widths = c(2, 1), guides = "collect") + plot_annotation(tag_levels = "A")
ggsave("figures/gwas_depth.pdf", fig, width = 10, height = 3.4)
ggsave("figures/gwas_depth.png", fig, width = 10, height = 3.4, dpi = 200)
