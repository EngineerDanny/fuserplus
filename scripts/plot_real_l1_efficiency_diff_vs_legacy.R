#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

build_dir <- get_arg("--build_dir", "build")
dataset_summary_csv <- get_arg("--dataset_summary_csv", "data/processed/grouped_regression_datasets_summary.csv")
out_png <- get_arg("--out_png", "inst/figures/real_l1_efficiency_diff_vs_legacy_chain_approx_only.png")
out_points_csv <- get_arg("--out_points_csv", "build/real_l1_efficiency_diff_vs_legacy_chain_approx_only_points.csv")
out_fold_points_csv <- get_arg("--out_fold_points_csv", "build/real_l1_efficiency_diff_vs_legacy_chain_approx_only_fold_points.csv")
out_skipped_csv <- get_arg("--out_skipped_csv", "build/real_l1_efficiency_diff_vs_legacy_chain_approx_only_skipped.csv")
new_settings_only <- tolower(get_arg("--new_settings_only", "true")) %in% c("true", "t", "1", "yes", "y")
methods_arg <- tolower(get_arg("--methods", "chain_approx"))

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_points_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_fold_points_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_skipped_csv), recursive = TRUE, showWarnings = FALSE)

valid_regimes <- c("real_standard", "real_high_dim")
required_methods <- c(
  "Legacy Full-Edge",
  "Operator (Ours)",
  "Dense Sort (Ours)",
  "Chain Approx. (Ours)",
  "Featureless"
)
compare_methods <- c("Operator (Ours)", "Dense Sort (Ours)", "Chain Approx. (Ours)")
method_key_map <- c(
  operator = "Operator (Ours)",
  dense_sort = "Dense Sort (Ours)",
  chain_approx = "Chain Approx. (Ours)"
)
method_keys <- trimws(unlist(strsplit(methods_arg, ",", fixed = TRUE)))
method_keys <- method_keys[nzchar(method_keys)]
if (length(method_keys) > 0L) {
  if (!all(method_keys %in% names(method_key_map))) {
    stop("Unknown --methods value. Use comma-separated keys from: operator,dense_sort,chain_approx")
  }
  compare_methods <- unname(method_key_map[method_keys])
}

summary_files <- list.files(
  build_dir,
  pattern = "^real_l1_cv_fixed_budget.*_summary\\.csv$",
  full.names = TRUE
)
if (length(summary_files) == 0L) {
  stop("No summary CSV files found.")
}

dataset_slugs <- character(0)
regime_map <- data.frame(dataset = character(0), regime = character(0), stringsAsFactors = FALSE)
if (file.exists(dataset_summary_csv)) {
  ds <- read.csv(dataset_summary_csv, stringsAsFactors = FALSE)
  if (all(c("dataset", "regime") %in% names(ds))) {
    ds <- ds[!is.na(ds$dataset), c("dataset", "regime"), drop = FALSE]
    regime_map <- unique(ds)
    dataset_slugs <- sort(unique(regime_map$dataset), decreasing = TRUE)
  }
}

extract_dataset_slug <- function(path, known_slugs) {
  bn <- basename(path)
  if (length(known_slugs) > 0L) {
    hits <- known_slugs[vapply(known_slugs, function(s) grepl(s, bn, fixed = TRUE), logical(1))]
    if (length(hits) > 0L) return(hits[[which.max(nchar(hits))]])
  }
  x <- sub("^real_l1_cv_fixed_budget_dense_", "", bn)
  x <- sub("_summary\\.csv$", "", x)
  x <- sub("_(std_intercept|system_time|mem|nomem|rerun|updated_dense_sort|mem_conserve|light).*", "", x)
  x <- sub("_k[0-9]+$", "", x)
  x
}

is_new_settings_file <- function(path) {
  bn <- basename(path)
  # Keep the current naming style and exclude known legacy/experimental tags.
  if (!grepl("^real_l1_cv_fixed_budget_dense_.*_summary\\.csv$", bn)) return(FALSE)
  legacy_tags <- c(
    "std_intercept_system_time",
    "nomem",
    "updated_dense_sort",
    "mem_conserve",
    "light",
    "rerun"
  )
  !any(vapply(legacy_tags, function(tag) grepl(tag, bn, fixed = TRUE), logical(1)))
}

skipped <- list()
parsed <- list()

for (f in summary_files) {
  if (isTRUE(new_settings_only) && !is_new_settings_file(f)) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = basename(f), reason = "not_new_settings_file", stringsAsFactors = FALSE)
    next
  }

  x <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !all(c("method_label", "mean_test_rmse", "mean_elapsed_sec", "mean_memory_mb", "sd_test_rmse", "sd_elapsed_sec", "sd_memory_mb") %in% names(x))) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = basename(f), reason = "missing_columns", stringsAsFactors = FALSE)
    next
  }

  x$method_label <- as.character(x$method_label)
  if (!all(required_methods %in% x$method_label)) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = basename(f), reason = "missing_required_methods", stringsAsFactors = FALSE)
    next
  }

  dataset <- extract_dataset_slug(f, dataset_slugs)
  regime <- "unknown"
  if (nrow(regime_map) > 0L && dataset %in% regime_map$dataset) {
    regime <- as.character(regime_map$regime[match(dataset, regime_map$dataset)])
  }
  if (!(regime %in% valid_regimes)) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = basename(f), reason = "invalid_or_unknown_regime", stringsAsFactors = FALSE)
    next
  }

  x <- x[match(required_methods, x$method_label), , drop = FALSE]
  x$dataset <- dataset
  x$regime <- regime
  x$source_file <- basename(f)
  x$file_mtime <- as.numeric(file.info(f)$mtime)
  parsed[[length(parsed) + 1L]] <- x
}

if (length(parsed) == 0L) {
  stop("No valid summary files to plot.")
}

all_df <- do.call(rbind, parsed)
all_df <- all_df[order(all_df$dataset, -all_df$file_mtime), , drop = FALSE]
latest_sources <- all_df$source_file[!duplicated(all_df$dataset)]
use_df <- all_df[all_df$source_file %in% latest_sources, , drop = FALSE]

# Ensure complete method sets after deduplication.
use_df <- do.call(rbind, lapply(split(use_df, use_df$dataset), function(df) {
  if (all(required_methods %in% df$method_label)) df else NULL
}))
if (is.null(use_df) || nrow(use_df) == 0L) {
  stop("No datasets with complete required method set.")
}

points <- list()
fold_points <- list()

for (ds_name in unique(use_df$dataset)) {
  d <- use_df[use_df$dataset == ds_name, , drop = FALSE]
  d <- d[match(required_methods, d$method_label), , drop = FALSE]

  feat_rmse <- as.numeric(d$mean_test_rmse[d$method_label == "Featureless"][1L])
  non_feat_rmse <- as.numeric(d$mean_test_rmse[d$method_label != "Featureless"])
  if (is.finite(feat_rmse) && length(non_feat_rmse) > 0L && feat_rmse < min(non_feat_rmse, na.rm = TRUE)) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = unique(d$source_file), reason = "featureless_beats_all", stringsAsFactors = FALSE)
    next
  }

  legacy <- d[d$method_label == "Legacy Full-Edge", , drop = FALSE][1L, ]
  l_rmse <- as.numeric(legacy$mean_test_rmse)
  if (!is.finite(l_rmse) || l_rmse <= 0) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = unique(d$source_file), reason = "invalid_legacy_baseline", stringsAsFactors = FALSE)
    next
  }

  raw_file <- file.path(build_dir, sub("_summary\\.csv$", "_raw.csv", unique(d$source_file)))
  raw <- tryCatch(read.csv(raw_file, stringsAsFactors = FALSE), error = function(e) NULL)
  need_raw_cols <- c("fold", "method_label", "elapsed_sec", "memory_mb", "test_rmse")
  if (is.null(raw) || !all(need_raw_cols %in% names(raw))) {
    skipped[[length(skipped) + 1L]] <- data.frame(file = unique(d$source_file), reason = "missing_or_invalid_raw_file", stringsAsFactors = FALSE)
    next
  }

  raw <- raw[raw$method_label %in% c("Legacy Full-Edge", compare_methods), need_raw_cols, drop = FALSE]
  raw$fold <- as.character(raw$fold)

  for (m in compare_methods) {
    r <- d[d$method_label == m, , drop = FALSE][1L, ]
    rmse <- as.numeric(r$mean_test_rmse)
    if (!is.finite(rmse) || rmse <= 0) next

    pct_rmse_diff <- 100 * (rmse - l_rmse) / l_rmse
    rmse_dir <- ifelse(pct_rmse_diff <= 0, "better_or_equal_rmse", "worse_rmse")

    legacy_fold <- raw[raw$method_label == "Legacy Full-Edge", c("fold", "elapsed_sec", "memory_mb"), drop = FALSE]
    names(legacy_fold) <- c("fold", "legacy_elapsed_sec", "legacy_memory_mb")
    method_fold <- raw[raw$method_label == m, c("fold", "elapsed_sec", "memory_mb"), drop = FALSE]
    names(method_fold) <- c("fold", "method_elapsed_sec", "method_memory_mb")
    fold_merged <- merge(method_fold, legacy_fold, by = "fold", all = FALSE, sort = TRUE)

    if (nrow(fold_merged) == 0L) next
    ok <- is.finite(fold_merged$method_elapsed_sec) &
      is.finite(fold_merged$legacy_elapsed_sec) &
      is.finite(fold_merged$method_memory_mb) &
      is.finite(fold_merged$legacy_memory_mb) &
      fold_merged$method_elapsed_sec > 0 &
      fold_merged$legacy_elapsed_sec > 0 &
      fold_merged$method_memory_mb > 0 &
      fold_merged$legacy_memory_mb > 0
    fold_merged <- fold_merged[ok, , drop = FALSE]
    if (nrow(fold_merged) == 0L) next

    fold_merged$delta_log10_time <- log10(fold_merged$method_elapsed_sec / fold_merged$legacy_elapsed_sec)
    fold_merged$pct_time_diff <- 100 * (fold_merged$method_elapsed_sec / fold_merged$legacy_elapsed_sec - 1)
    fold_merged$delta_log10_memory <- log10(fold_merged$method_memory_mb / fold_merged$legacy_memory_mb)
    fold_merged$dataset <- ds_name
    fold_merged$regime <- as.character(r$regime)
    fold_merged$method_label <- m
    fold_merged$rmse_dir <- rmse_dir
    fold_merged$source_file <- as.character(r$source_file)

    fold_points[[length(fold_points) + 1L]] <- fold_merged[, c(
      "dataset", "regime", "method_label", "fold",
      "pct_time_diff", "delta_log10_time", "delta_log10_memory", "rmse_dir", "source_file"
    )]

    points[[length(points) + 1L]] <- data.frame(
      dataset = ds_name,
      regime = as.character(r$regime),
      method_label = m,
      pct_time_diff = mean(fold_merged$pct_time_diff, na.rm = TRUE),
      sd_pct_time_diff = sd(fold_merged$pct_time_diff, na.rm = TRUE),
      delta_log10_time = mean(fold_merged$delta_log10_time, na.rm = TRUE),
      delta_log10_memory = mean(fold_merged$delta_log10_memory, na.rm = TRUE),
      sd_delta_log10_time = sd(fold_merged$delta_log10_time, na.rm = TRUE),
      sd_delta_log10_memory = sd(fold_merged$delta_log10_memory, na.rm = TRUE),
      n_folds = nrow(fold_merged),
      pct_rmse_diff = pct_rmse_diff,
      abs_pct_rmse_diff = abs(pct_rmse_diff),
      rmse_dir = rmse_dir,
      source_file = as.character(r$source_file),
      stringsAsFactors = FALSE
    )
  }
}

plot_df <- if (length(points) > 0L) do.call(rbind, points) else NULL
if (is.null(plot_df) || nrow(plot_df) == 0L) {
  stop("No plot points after filtering.")
}

write.csv(plot_df, out_points_csv, row.names = FALSE)
if (length(fold_points) > 0L) {
  fold_df <- do.call(rbind, fold_points)
  write.csv(fold_df, out_fold_points_csv, row.names = FALSE)
}
if (length(skipped) > 0L) {
  skipped_df <- unique(do.call(rbind, skipped))
  write.csv(skipped_df, out_skipped_csv, row.names = FALSE)
}

plot_df$method_label <- factor(plot_df$method_label, levels = compare_methods)
plot_df$rmse_dir <- factor(plot_df$rmse_dir, levels = c("better_or_equal_rmse", "worse_rmse"))
plot_df$dataset <- factor(plot_df$dataset, levels = sort(unique(plot_df$dataset)))
fold_df <- if (length(fold_points) > 0L) do.call(rbind, fold_points) else NULL
if (!is.null(fold_df) && nrow(fold_df) > 0L) {
  fold_df$method_label <- factor(fold_df$method_label, levels = compare_methods)
  fold_df$rmse_dir <- factor(fold_df$rmse_dir, levels = c("better_or_equal_rmse", "worse_rmse"))
  fold_df$dataset <- factor(fold_df$dataset, levels = levels(plot_df$dataset))
}

p <- ggplot(plot_df, aes(x = delta_log10_time, y = delta_log10_memory)) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  {
    if (!is.null(fold_df) && nrow(fold_df) > 0L) {
      geom_point(
        data = fold_df,
        aes(x = delta_log10_time, y = delta_log10_memory, color = dataset),
        shape = 16,
        size = 1.1,
        alpha = 0.18,
        inherit.aes = FALSE
      )
    }
  } +
  geom_segment(
    aes(
      x = delta_log10_time - sd_delta_log10_time,
      xend = delta_log10_time + sd_delta_log10_time,
      y = delta_log10_memory,
      yend = delta_log10_memory,
      color = dataset
    ),
    linewidth = 0.5,
    alpha = 0.95
  ) +
  geom_segment(
    aes(
      x = delta_log10_time,
      xend = delta_log10_time,
      y = delta_log10_memory - sd_delta_log10_memory,
      yend = delta_log10_memory + sd_delta_log10_memory,
      color = dataset
    ),
    linewidth = 0.35,
    alpha = 0.9
  ) +
  geom_point(
    aes(fill = rmse_dir),
    shape = 21,
    color = "black",
    stroke = 0.2,
    size = 2.4,
    alpha = 0.9
  ) +
  ggrepel::geom_label_repel(
    aes(label = dataset, color = dataset),
    size = 2.7,
    label.size = 0.15,
    alpha = 0.85,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.alpha = 0.85,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      better_or_equal_rmse = "black",
      worse_rmse = "red3"
    )
  ) +
  facet_wrap(~method_label) +
  labs(
    title = "L1 Efficiency Difference vs Legacy Full-Edge",
    subtitle = "Faint dots: fold-level ratios. Solid dot + bars: mean +/- SD across folds",
    x = "log10(Time / Legacy Time)   (left is faster)",
    y = "log10(Memory / Legacy Memory)   (down is lower memory)",
    fill = "RMSE direction"
  ) +
  guides(color = "none")

ggsave(out_png, plot = p, width = 7, height = 5.5, dpi = 300)

cat("Saved plot: ", out_png, "\n", sep = "")
cat("Saved points: ", out_points_csv, "\n", sep = "")
if (file.exists(out_fold_points_csv)) {
  cat("Saved fold points: ", out_fold_points_csv, "\n", sep = "")
}
cat("Datasets plotted: ", length(unique(plot_df$dataset)), "\n", sep = "")
if (file.exists(out_skipped_csv)) {
  cat("Saved skipped list: ", out_skipped_csv, "\n", sep = "")
}
