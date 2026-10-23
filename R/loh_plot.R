## loh_plot.R -- summarise out/loh/res/*.rds (R/loh_poc.R) into a table and figure.
## run from the repository root: Rscript R/loh_plot.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
files <- list.files("out/loh/res", "\\.rds$", full.names = TRUE)
lab <- function(f) sub("^C([0-9.]+)_.*", "\\1", basename(f))
truth <- rbindlist(lapply(files, function(f) readRDS(f)$truth[, C := lab(f)]))
calls <- rbindlist(lapply(files, function(f) readRDS(f)$calls[, C := lab(f)]))
t <- truth[testable == TRUE]
t[, het_bin := cut(n_famhet, c(0, 10, 20, 40, 80, 160, Inf),
                   labels = c("1-10", "11-20", "21-40", "41-80", "81-160", ">160"))]
rec <- t[, .(recall = mean(detected), n = .N, se = sqrt(mean(detected) * (1 - mean(detected)) / .N)),
         by = .(C, het_bin)]
prec <- calls[, .(precision = mean(correct), calls = .N), by = C]
tab <- merge(t[, .(testable_tracts = .N, recall = mean(detected),
                   recall_gt40 = mean(detected[n_famhet > 40])), by = C], prec, by = "C")
fwrite(tab, "out/loh/summary_loh.csv"); print(tab)
m <- t[detected == TRUE]
cat(sprintf("tract length: log-scale r = %.2f, median called/true = %.2f (n = %d)\n",
            cor(log(m$len), log(m$called_len)), median(m$called_len / m$len), nrow(m)))

cols <- c("0.99" = "#9ec5f4", "0.995" = "#2a78d6", "0.999" = "#174f94")
th <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA))
pA <- ggplot(rec, aes(het_bin, recall, colour = C, group = C)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.7) +
  scale_y_continuous("Tracts detected", limits = c(0, 1)) +
  scale_x_discrete("Heterozygous sites lost in the tract") +
  scale_colour_manual("Cloning rate", values = cols) + th
pB <- ggplot(m, aes(len / 1e3, called_len / 1e3, colour = C)) +
  geom_abline(linetype = "dashed", colour = "grey60") + geom_point(size = 0.8, alpha = 0.6, show.legend = FALSE) +
  scale_x_log10("True tract length (kb)") + scale_y_log10("Estimated tract length (kb)") +
  scale_colour_manual("Cloning rate", values = cols) + th
fig <- pA + pB + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave("figures/loh_poc.pdf", fig, width = 8, height = 3.2)
ggsave("figures/loh_poc.png", fig, width = 8, height = 3.2, dpi = 200)
