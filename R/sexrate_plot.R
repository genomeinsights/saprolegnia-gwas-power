## sexrate_plot.R -- which LD summaries track the rate of sex? Run from the repo root.
## x-axis: sexual births per deme per generation, N * (1 - C).
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
f <- list.files("out/sexrate/res", "\\.rds$", full.names = TRUE)
x <- rbindlist(lapply(f, function(ff) {
  m <- regmatches(basename(ff), regexec("^C([0-9.]+)_N([0-9]+)_r([0-9]+)", basename(ff)))[[1]]
  readRDS(ff)[, `:=`(C = as.numeric(m[2]), N = as.integer(m[3]), rep = as.integer(m[4]))]
}), fill = TRUE)
x[, sex_events := N * (1 - C)]
fwrite(x, "out/sexrate/summary_sexrate.csv")
long <- melt(x, id.vars = c("C", "N", "rep", "sex_events"),
             measure.vars = c("rbarD", "r2_between", "r2_10kb", "a_pred", "n_clones"))
## how well does each summary predict log(sex events) vs log(1 - C)?  (R2 of log-log fits)
fit <- long[is.finite(value) & value > 0, .(
  r2_sex_events = summary(lm(log(value) ~ log(sex_events)))$r.squared,
  r2_rate_of_sex = summary(lm(log(value) ~ log(1 - C)))$r.squared), by = variable]
print(fit)
lab <- c(rbarD = "rbarD (index of assoc.)", r2_between = "r2 between chrom.",
         r2_10kb = "r2 within 10 kb", a_pred = "LDscnR decay rate",
         n_clones = "Distinct clones")
long[, variable := factor(lab[as.character(variable)], levels = lab)]
cols <- c("500" = "#9ec5f4", "1000" = "#2a78d6", "2000" = "#174f94")
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
p <- ggplot(long[is.finite(value) & value > 0], aes(sex_events, value, colour = factor(N))) +
  geom_point(size = 1.4, alpha = 0.8) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  scale_x_log10("Sexual births per site per generation, N x (1 - cloning rate)") + scale_y_log10(NULL) +
  scale_colour_manual("Individuals\nper site (N)", values = cols) + th
ggsave("figures/sexrate_ld.pdf", p, width = 11, height = 2.8)
ggsave("figures/sexrate_ld.png", p, width = 11, height = 2.8, dpi = 200)
