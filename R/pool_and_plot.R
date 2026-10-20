## pool_and_plot.R -- pool out/res/*.res.rds into power / FDR summaries and the figure.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
## run from the repository root: Rscript R/pool_and_plot.R
files <- list.files("out/res", "\\.res\\.rds$", full.names = TRUE)
x <- rbindlist(lapply(files, function(f) {
  r <- readRDS(f); C <- as.numeric(sub("^C([0-9.]+)_r.*", "\\1", basename(f)))
  r$res[, `:=`(C = C, file = basename(f), n_true = r$clone$n_true, n_est = r$clone$n_est,
               ari = r$clone$ari, n_snp = r$n_snp)]
}))
## The Simes arm is computed in R/analyse_rep.R but left out of the summaries and figure.
METHODS <- c(emmax_snp = "EMMAX, single SNP", ldscnr_consensus = "LDscnR, LD units")
x <- x[method %in% names(METHODS)]
x[, method := factor(METHODS[method], levels = METHODS)]

clones <- unique(x[, .(C, file, n_true, n_est, ari, n_snp)])
cs <- clones[, .(reps = .N, n_true = mean(n_true), n_est = mean(n_est), n_est_min = min(n_est),
                 n_est_max = max(n_est), ari = mean(ari), n_snp = mean(n_snp)), by = C][order(C)]
pw <- x[, .(n_est = mean(n_est), power = mean(detected), se = sqrt(mean(detected) * (1 - mean(detected)) / .N),
            n = .N), by = .(C, method, pve)]
fd <- x[, .(n_est = mean(n_est), regions = sum(n_regions), false = sum(n_false),
            fdr = sum(n_false) / pmax(sum(n_regions), 1), false_per_scan = mean(n_false)), by = .(C, method)]
fwrite(cs, "out/summary_clones.csv"); fwrite(pw, "out/summary_power.csv"); fwrite(fd, "out/summary_fdr.csv")
print(cs); print(dcast(pw, method + pve ~ C, value.var = "power")); print(fd)

## ---- figure -------------------------------------------------------------------------
blues <- c("#9ec5f4", "#5c9be3", "#2a78d6", "#174f94")
cats  <- c("#2a78d6", "#eb6834")
th <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))
xs <- scale_x_continuous("Distinct clones among 300 isolates (inferred)", breaks = round(cs$n_est),
                         sec.axis = dup_axis(name = "Cloning rate", labels = cs$C))
pA <- ggplot(pw, aes(n_est, power, colour = factor(pve), group = factor(pve))) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = pmax(power - se, 0), ymax = pmin(power + se, 1)), width = 0, linewidth = 0.3) +
  facet_wrap(~method) +
  scale_x_continuous("Distinct clones among 300 isolates (inferred)", breaks = round(cs$n_est)) +
  scale_y_continuous("Power (QTL detected)", limits = c(0, 1)) +
  scale_colour_manual("Variance explained by QTL", values = blues,
                      labels = function(v) paste0(as.numeric(v) * 100, "%")) + th
pB <- ggplot(fd, aes(n_est, fdr, colour = method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + xs +
  scale_y_continuous("False discovery proportion", limits = c(0, NA)) +
  scale_colour_manual(NULL, values = cats) + th
fig <- pA / pB + plot_layout(heights = c(1.3, 1)) + plot_annotation(tag_levels = "A")
ggsave("figures/gwas_power_clonal.pdf", fig, width = 6, height = 5.6)
ggsave("figures/gwas_power_clonal.png", fig, width = 6, height = 5.6, dpi = 200)
