## sexrate_loo.R -- leave-one-out prediction of the rate of sex (1 - cloning rate) from
## LD summaries: r-bar-d alone, plus LDscnR's decay rate, plus clone count (C <= 0.995). Run from the repo root.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
x <- fread("out/sexrate/summary_sexrate.csv")[is.finite(a_pred) & C <= 0.995]
x[, y := log10(1 - C)]
models <- list(
  "A. Index of association (rbarD)" = y ~ log(rbarD),
  "B. rbarD + LDscnR LD-decay rate" = y ~ log(rbarD) + log(a_pred),
  "C. rbarD + LDscnR decay + clone count" = y ~ log(rbarD) + log(a_pred) + log(n_clones))
loo <- rbindlist(lapply(names(models), function(m) {
  pr <- vapply(seq_len(nrow(x)), function(i) predict(lm(models[[m]], data = x[-i]), x[i]), 0)
  data.table(model = m, N = x$N, true = x$y, pred = pr)
}))
acc <- loo[, .(rmse = sqrt(mean((pred - true)^2)), within_2x = mean(abs(pred - true) <= log10(2)),
               r2 = cor(pred, true)^2), by = model]
print(acc); fwrite(acc, "out/sexrate/loo_accuracy.csv")
loo[, model := factor(model, levels = names(models))]
cols <- c("500" = "#9ec5f4", "1000" = "#2a78d6", "2000" = "#174f94")
p <- ggplot(loo, aes(10^true, 10^pred, colour = factor(N))) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_jitter(width = 0.03, height = 0, size = 1.6, alpha = 0.85) +
  facet_wrap(~ model) +
  scale_x_log10("True rate of sex (1 - cloning rate)") + scale_y_log10("Predicted (leave-one-out)") +
  scale_colour_manual("Individuals\nper site", values = cols) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/sexrate_loo.pdf", p, width = 9, height = 3.2)
ggsave("figures/sexrate_loo.png", p, width = 9, height = 3.2, dpi = 200)
