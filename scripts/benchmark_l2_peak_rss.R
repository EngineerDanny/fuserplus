#!/usr/bin/env Rscript

# Peak-memory benchmark for old vs new L2.
# Uses `/usr/bin/time -l` to capture process-level max RSS.

DEFAULTS <- list(
  seed = 20260206L,
  k = 60L,
  p = 120L,
  n_group_train = 35L,
  n_group_test = 15L,
  sigma = 0.05,
  lambda = 1e-3,
  gamma = 1e-3,
  scaling = FALSE,
  g_structure = "sparse_chain", # sparse_chain | dense
  variants = "old,new",
  out_csv = "build/benchmark_l2_peak_rss.csv"
)

parse_args <- function(argv) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    key <- substring(key, 3L)
    if (i == length(argv) || startsWith(argv[[i + 1L]], "--")) {
      out[[key]] <- "TRUE"
      i <- i + 1L
    } else {
      out[[key]] <- argv[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

as_chr_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_first_num <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  m <- regexec(pattern, hit[1L], perl = TRUE)
  g <- regmatches(hit[1L], m)[[1L]]
  if (length(g) < 2L) return(NA_real_)
  as.numeric(g[2L])
}

run_variant <- function(variant, cfg) {
  tmp_csv <- tempfile(pattern = sprintf("l2_%s_", variant), fileext = ".csv")
  args <- c(
    "-l",
    "Rscript",
    "scripts/benchmark_l2_old_new.R",
    "--mode", "run",
    "--variant", variant,
    "--seed", as.character(cfg$seed),
    "--k", as.character(cfg$k),
    "--p", as.character(cfg$p),
    "--n_group_train", as.character(cfg$n_group_train),
    "--n_group_test", as.character(cfg$n_group_test),
    "--sigma", as.character(cfg$sigma),
    "--lambda", as.character(cfg$lambda),
    "--gamma", as.character(cfg$gamma),
    "--scaling", ifelse(cfg$scaling, "true", "false"),
    "--g_structure", cfg$g_structure,
    "--output_csv", tmp_csv
  )

  lines <- system2("/usr/bin/time", args = args, stdout = TRUE, stderr = TRUE)

  max_rss <- parse_first_num(lines, "^\\s*([0-9]+)\\s+maximum resident set size\\s*$")
  peak_foot <- parse_first_num(lines, "^\\s*([0-9]+)\\s+peak memory footprint\\s*$")
  elapsed <- parse_first_num(lines, "Elapsed time \\(seconds\\)=([0-9.]+)")
  test_rmse <- parse_first_num(lines, "Test RMSE=([0-9.eE+-]+)")

  data.frame(
    variant = variant,
    seed = cfg$seed,
    k = cfg$k,
    p = cfg$p,
    n_group_train = cfg$n_group_train,
    n_group_test = cfg$n_group_test,
    sigma = cfg$sigma,
    lambda = cfg$lambda,
    gamma = cfg$gamma,
    scaling = cfg$scaling,
    g_structure = cfg$g_structure,
    elapsed_sec = elapsed,
    test_rmse = test_rmse,
    max_rss_bytes = max_rss,
    max_rss_mb = max_rss / (1024^2),
    peak_footprint_bytes = peak_foot,
    peak_footprint_mb = peak_foot / (1024^2),
    stringsAsFactors = FALSE
  )
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
cfg <- list(
  seed = as.integer(argv$seed %||% DEFAULTS$seed),
  k = as.integer(argv$k %||% DEFAULTS$k),
  p = as.integer(argv$p %||% DEFAULTS$p),
  n_group_train = as.integer(argv$n_group_train %||% DEFAULTS$n_group_train),
  n_group_test = as.integer(argv$n_group_test %||% DEFAULTS$n_group_test),
  sigma = as.numeric(argv$sigma %||% DEFAULTS$sigma),
  lambda = as.numeric(argv$lambda %||% DEFAULTS$lambda),
  gamma = as.numeric(argv$gamma %||% DEFAULTS$gamma),
  scaling = as_bool(argv$scaling %||% as.character(DEFAULTS$scaling), default = FALSE),
  g_structure = as.character(argv$g_structure %||% DEFAULTS$g_structure),
  variants = as_chr_vec(as.character(argv$variants %||% DEFAULTS$variants)),
  out_csv = as.character(argv$out_csv %||% DEFAULTS$out_csv)
)

if (!all(cfg$variants %in% c("old", "new"))) {
  stop("variants must be a comma-separated subset of: old,new")
}

rows <- lapply(cfg$variants, run_variant, cfg = cfg)
out <- do.call(rbind, rows)

dir.create(dirname(cfg$out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, cfg$out_csv, row.names = FALSE)

cat("Peak RSS benchmark summary:\n")
print(out[, c("variant", "elapsed_sec", "test_rmse", "max_rss_mb", "peak_footprint_mb")], row.names = FALSE)
cat(sprintf("\nSaved: %s\n", cfg$out_csv))
