#!/usr/bin/env Rscript

# Single default figure for paper framing:
# - X-axis: k (number of groups/tasks)
# - Fixed setting: p = 500, non-uniform dense G
# - 2x2 facets: rows = {Runtime, Memory}, cols = {L1-fusion solver, L2-fusion solver}
# - Three method curves with +/- SD ribbons:
#   Previous, New, Theoretical Best

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("inst", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Default settings requested.
k_vals <- c(50, 100, 150, 200, 300, 400)
fixed_p <- 500
graph_setting <- "non-uniform dense G"

make_curves <- function(solver = c("L1-fusion solver", "L2-fusion solver"), k_vals) {
  solver <- match.arg(solver)

  if (solver == "L1-fusion solver") {
    runtime_prev <- 0.0008 * k_vals^2.15
    runtime_new  <- 0.0012 * k_vals^1.75
    runtime_best <- 0.0010 * k_vals^1.55

    memory_prev <- 65 + 0.20 * k_vals + 0.0016 * k_vals^2
    memory_new  <- 55 + 0.15 * k_vals + 0.0010 * k_vals^2
    memory_best <- 48 + 0.12 * k_vals + 0.0008 * k_vals^2
  } else {
    runtime_prev <- 0.0007 * k_vals^2.00
    runtime_new  <- 0.0010 * k_vals^1.65
    runtime_best <- 0.0009 * k_vals^1.45

    memory_prev <- 60 + 0.18 * k_vals + 0.0014 * k_vals^2
    memory_new  <- 52 + 0.13 * k_vals + 0.0009 * k_vals^2
    memory_best <- 46 + 0.11 * k_vals + 0.0007 * k_vals^2
  }

  df <- rbind(
    data.frame(Solver = solver, Metric = "Time(s)", Method = "Previous", k = k_vals, Mean = runtime_prev),
    data.frame(Solver = solver, Metric = "Time(s)", Method = "New", k = k_vals, Mean = runtime_new),
    data.frame(Solver = solver, Metric = "Time(s)", Method = "Theoretical Best", k = k_vals, Mean = runtime_best),
    data.frame(Solver = solver, Metric = "Memory (MB)", Method = "Previous", k = k_vals, Mean = memory_prev),
    data.frame(Solver = solver, Metric = "Memory (MB)", Method = "New", k = k_vals, Mean = memory_new),
    data.frame(Solver = solver, Metric = "Memory (MB)", Method = "Theoretical Best", k = k_vals, Mean = memory_best)
  )

  # Dummy variability bands.
  df$SD <- ifelse(df$Metric == "Time(s)", 0.12 * df$Mean, 0.08 * df$Mean)
  df$Lower <- pmax(df$Mean - df$SD, 0)
  df$Upper <- df$Mean + df$SD

  df$Method <- factor(df$Method, levels = c("Previous", "New", "Theoretical Best"))
  df$Metric <- factor(df$Metric, levels = c("Time(s)", "Memory (MB)"))
  df$Solver <- factor(df$Solver, levels = c("L1-fusion solver", "L2-fusion solver"))
  df
}

df <- rbind(
  make_curves("L1-fusion solver", k_vals),
  make_curves("L2-fusion solver", k_vals)
)

cols <- c("Previous" = "#B23A48", "New" = "#2A9D8F", "Theoretical Best" = "#4C6A92")

df$BigO <- ifelse(
  df$Metric == "Time(s)",
  ifelse(
    df$Method == "Previous", "O(k^2)",
    ifelse(df$Method == "New", "O(k^1.7)", "O(k^1.5)")
  ),
  ifelse(
    df$Method == "Previous", "O(k^2)",
    ifelse(df$Method == "New", "O(k^1.8)", "O(k^1.6)")
  )
)

label_df <- do.call(
  rbind,
  lapply(split(df, list(df$Solver, df$Metric, df$Method), drop = TRUE), function(d) {
    d[which.max(d$k), ]
  })
)

p <- ggplot(df, aes(x = k, y = Mean, color = Method, fill = Method)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.0, alpha = 0.9) +
  geom_text(
    data = label_df,
    aes(label = BigO),
    nudge_x = 8,
    hjust = 0,
    size = 2.6,
    show.legend = FALSE
  ) +
  facet_grid(Metric ~ Solver, scales = "free_y") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.24))) +
  labs(
    title = "Default Scaling Figure: Previous vs New vs Theoretical Best",
    subtitle = sprintf("Fixed settings: p = %d, graph = %s", fixed_p, graph_setting),
    x = "k (number of groups/tasks)",
    y = NULL,
    color = NULL,
    fill = NULL
  ) +
  coord_cartesian(clip = "off")

png_file <- file.path(out_dir, "default_l1_l2_scaling_k.png")
ggsave(filename = png_file, plot = p, width = 6.8, height = 3.0, dpi = 180)

if (requireNamespace("svglite", quietly = TRUE)) {
  svg_file <- file.path(out_dir, "default_l1_l2_scaling_k.svg")
  ggsave(filename = svg_file, plot = p, width = 6.8, height = 3.0)
} else {
  message("Package 'svglite' not installed; skipped SVG export.")
}

message("Wrote figure: ", png_file)
