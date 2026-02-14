#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("inst", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

k_vals <- c(20, 40, 60, 80, 100)

make_base_df <- function() {
  df <- rbind(
    data.frame(Solver = "L1-fusion solver", Metric = "Time(s)", Method = "old", k = k_vals, Mean = 0.0006 * k_vals^2.05),
    data.frame(Solver = "L1-fusion solver", Metric = "Time(s)", Method = "new", k = k_vals, Mean = 0.0010 * k_vals^1.65),
    data.frame(Solver = "L1-fusion solver", Metric = "Memory (MB)", Method = "old", k = k_vals, Mean = 25 + 0.022 * k_vals^2),
    data.frame(Solver = "L1-fusion solver", Metric = "Memory (MB)", Method = "new", k = k_vals, Mean = 20 + 0.010 * k_vals^2),
    data.frame(Solver = "L2-fusion solver", Metric = "Time(s)", Method = "old", k = k_vals, Mean = 0.0003 * k_vals^2.25),
    data.frame(Solver = "L2-fusion solver", Metric = "Time(s)", Method = "new", k = k_vals, Mean = 0.0007 * k_vals^1.45),
    data.frame(Solver = "L2-fusion solver", Metric = "Memory (MB)", Method = "old", k = k_vals, Mean = 10 + 0.018 * k_vals^2.2),
    data.frame(Solver = "L2-fusion solver", Metric = "Memory (MB)", Method = "new", k = k_vals, Mean = 5 + 0.004 * k_vals^2)
  )
  df$SD <- ifelse(df$Metric == "Time(s)", 0.10 * df$Mean, 0.07 * df$Mean)
  df$Lower <- pmax(0, df$Mean - df$SD)
  df$Upper <- df$Mean + df$SD
  df$Method <- factor(df$Method, levels = c("old", "new"))
  df$Metric <- factor(df$Metric, levels = c("Time(s)", "Memory (MB)"))
  df$Solver <- factor(df$Solver, levels = c("L1-fusion solver", "L2-fusion solver"))
  df
}

df <- make_base_df()
cols <- c(old = "#B23A48", new = "#2A9D8F")

# Option 1: standard 2x2 scaling with ribbons.
p1 <- ggplot(df, aes(k, Mean, color = Method, fill = Method)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.16, linewidth = 0) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  facet_grid(Metric ~ Solver, scales = "free_y") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(title = "Option 1: Scaling Curves", x = "k", y = NULL, color = NULL, fill = NULL)

# Option 2: speedup and memory-ratio panels.
ratio_df <- do.call(rbind, lapply(split(df, list(df$Solver, df$Metric), drop = TRUE), function(d) {
  old <- d[d$Method == "old", ]
  new <- d[d$Method == "new", ]
  data.frame(
    Solver = old$Solver,
    Metric = old$Metric,
    k = old$k,
    Value = old$Mean / new$Mean
  )
}))
ratio_df$Panel <- ifelse(ratio_df$Metric == "Time(s)", "Speedup (old/new)", "Memory ratio (old/new)")
ratio_df$Panel <- factor(ratio_df$Panel, levels = c("Speedup (old/new)", "Memory ratio (old/new)"))

p2 <- ggplot(ratio_df, aes(k, Value, color = Solver)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  facet_wrap(~ Panel, scales = "free_y", nrow = 1) +
  labs(title = "Option 2: Relative Gain View", x = "k", y = NULL, color = NULL)

# Option 3: log-log time only with Big-O text.
time_df <- df[df$Metric == "Time(s)", ]
labels3 <- do.call(rbind, lapply(split(time_df, list(time_df$Solver, time_df$Method), drop = TRUE), function(d) d[which.max(d$k), ]))
labels3$BigO <- ifelse(labels3$Method == "old", "O(k^2)", "O(k^1.5~1.7)")

p3 <- ggplot(time_df, aes(k, Mean, color = Method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.9) +
  geom_text(data = labels3, aes(label = BigO), nudge_x = 4, hjust = 0, size = 2.8, show.legend = FALSE) +
  facet_wrap(~ Solver, scales = "free_y", nrow = 1) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.25))) +
  scale_y_log10() +
  scale_color_manual(values = cols) +
  labs(title = "Option 3: Log-Log Time With Big-O Tags", x = "k (log scale)", y = "Time(s) (log scale)", color = NULL)

# Option 4: endpoint comparison bars (largest k).
end_df <- df[df$k == max(df$k), c("Solver", "Metric", "Method", "Mean")]
p4 <- ggplot(end_df, aes(Method, Mean, fill = Method)) +
  geom_col(width = 0.66) +
  facet_grid(Metric ~ Solver, scales = "free_y") +
  scale_fill_manual(values = cols) +
  labs(title = "Option 4: Endpoint Comparison (largest k)", x = NULL, y = NULL, fill = NULL)

# Option 5: normalized to old=1 baseline.
norm_df <- do.call(rbind, lapply(split(df, list(df$Solver, df$Metric), drop = TRUE), function(d) {
  old <- d[d$Method == "old", c("k", "Mean")]
  d$Normalized <- d$Mean / old$Mean[match(d$k, old$k)]
  d
}))
p5 <- ggplot(norm_df, aes(k, Normalized, color = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  facet_grid(Metric ~ Solver, scales = "free_y") +
  scale_color_manual(values = cols) +
  labs(title = "Option 5: Normalized Curves (old = 1)", x = "k", y = "Normalized value", color = NULL)

ggsave(file.path(out_dir, "dummy_option_1_scaling.png"), p1, width = 6.8, height = 3.0, dpi = 180)
ggsave(file.path(out_dir, "dummy_option_2_relative_gain.png"), p2, width = 6.8, height = 3.0, dpi = 180)
ggsave(file.path(out_dir, "dummy_option_3_loglog_time.png"), p3, width = 6.8, height = 3.0, dpi = 180)
ggsave(file.path(out_dir, "dummy_option_4_endpoint_bars.png"), p4, width = 6.8, height = 3.0, dpi = 180)
ggsave(file.path(out_dir, "dummy_option_5_normalized.png"), p5, width = 6.8, height = 3.0, dpi = 180)

cat("Wrote 5 dummy figure options to inst/figures/\n")
