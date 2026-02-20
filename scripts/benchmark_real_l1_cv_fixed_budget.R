#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    library(RSpectra)
  }
  if (requireNamespace("irlba", quietly = TRUE)) {
    library(irlba)
  }
})

# Load local L1 solver implementations.
source("R/l1_fusion_new_utils.R")
source("R/l1_fusion.R")
source("R/l1_fusion_operator_new.R")
source("R/l1_fusion_dfs_chain.R")
source("R/l1_fusion_chain_specialized.R")
source("R/l1_fusion_dense_sort.R")
source("R/l1_fusion_new.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

DEFAULTS <- list(
  data_rds = "data/processed/communities_crime/grouped_regression.rds",
  l1_methods = "legacy_full_edge,operator,dense_sort,chain_approx,featureless",
  l1_budget = 1300L,
  folds = 5L,
  seed = 20260217L,
  min_group_size = 5L,
  k_keep = NA,
  g_structure = "dense", # dense | sparse_chain
  lambda = 1e-3,
  gamma = 1e-3,
  mu = 1e-4,
  tol = 1e-5,
  intercept = TRUE,
  scaling = FALSE,
  standardize = TRUE,
  conserve_memory = TRUE,
  measure_memory = TRUE,
  edge_block = 256L,
  c_flag_old = FALSE,
  time_log_scale = FALSE,
  out_prefix = "real_l1_cv_fixed_budget"
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
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  vals <- as.numeric(parts)
  if (any(is.na(vals))) stop("Could not parse numeric vector.")
  vals
}

as_int <- function(x) {
  out <- as.integer(x)
  if (is.na(out)) stop("Could not parse integer value.")
  out
}

as_chr_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

check_choices <- function(values, allowed, arg_name) {
  bad <- setdiff(values, allowed)
  if (length(bad)) {
    stop(sprintf(
      "Invalid %s value(s): %s. Allowed: %s",
      arg_name, paste(bad, collapse = ", "), paste(allowed, collapse = ", ")
    ))
  }
}

predict_from_beta <- function(beta, X, groups) {
  group_levels <- sort(unique(groups))
  cols <- if (!is.null(colnames(beta)) && all(as.character(group_levels) %in% colnames(beta))) {
    match(as.character(groups), colnames(beta))
  } else {
    as.integer(factor(groups, levels = group_levels))
  }

  p <- ncol(X)
  has_intercept <- nrow(beta) == (p + 1L)
  beta_sub <- beta[, cols, drop = FALSE]
  if (has_intercept) {
    rowSums(X * t(beta_sub[seq_len(p), , drop = FALSE])) + as.numeric(beta_sub[p + 1L, ])
  } else {
    rowSums(X * t(beta_sub))
  }
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

stratified_folds <- function(groups, nfold = 5L, seed = 1L) {
  set.seed(seed)
  folds <- integer(length(groups))
  split_idx <- split(seq_along(groups), groups)
  for (ii in split_idx) {
    perm <- sample(ii)
    tag <- rep(seq_len(nfold), length.out = length(perm))
    folds[perm] <- tag
  }
  folds
}

prepare_real_data <- function(data_rds, min_group_size = 5L, k_keep = NA_integer_) {
  d0 <- readRDS(data_rds)
  if (!all(c("X", "y", "groups") %in% names(d0))) {
    stop("RDS must contain at least X, y, groups.")
  }

  X <- as.matrix(d0$X)
  y <- as.numeric(d0$y)
  groups <- as.integer(d0$groups)
  if (nrow(X) != length(y) || nrow(X) != length(groups)) {
    stop("RDS has inconsistent X/y/groups dimensions.")
  }
  if (!all(is.finite(y))) {
    stop("Response y contains non-finite values.")
  }
  if (stats::sd(y) <= 0) {
    stop("Response y is constant; choose a dataset with response variation.")
  }

  group_counts <- sort(table(groups), decreasing = TRUE)
  keep_levels <- as.integer(names(group_counts[group_counts >= min_group_size]))
  if (length(keep_levels) < 2L) {
    stop("Need at least 2 groups after min_group_size filtering.")
  }
  if (!is.na(k_keep) && is.finite(k_keep) && k_keep > 0L && length(keep_levels) > k_keep) {
    keep_levels <- keep_levels[seq_len(k_keep)]
  }

  keep_rows <- groups %in% keep_levels
  X <- X[keep_rows, , drop = FALSE]
  y <- y[keep_rows]
  groups_orig <- groups[keep_rows]
  groups <- as.integer(factor(groups_orig, levels = keep_levels))
  group_levels_orig <- sort(unique(groups_orig))

  list(
    X = X,
    y = y,
    groups = groups,
    group_levels_orig = group_levels_orig,
    G_dense = d0$G_dense %||% NULL,
    G_sparse_chain = d0$G_sparse_chain %||% NULL
  )
}

build_chain_graph <- function(k) {
  G <- matrix(0, k, k)
  if (k > 1L) {
    for (i in 1:(k - 1L)) {
      G[i, i + 1L] <- 1
      G[i + 1L, i] <- 1
    }
  }
  G
}

subset_graph_l1 <- function(data_obj, g_structure) {
  k <- length(unique(data_obj$groups))
  lev <- data_obj$group_levels_orig
  if (g_structure == "dense") {
    if (!is.null(data_obj$G_dense)) {
      G <- as.matrix(data_obj$G_dense[lev, lev, drop = FALSE])
      diag(G) <- 0
      return(G)
    }
    G <- matrix(1, k, k)
    diag(G) <- 0
    return(G)
  }

  if (g_structure == "sparse_chain") {
    if (!is.null(data_obj$G_sparse_chain)) {
      G <- as.matrix(data_obj$G_sparse_chain[lev, lev, drop = FALSE])
      diag(G) <- 0
      return(G)
    }
    return(build_chain_graph(k))
  }

  stop("Unknown g_structure: ", g_structure)
}

fit_l1 <- function(method, d, cfg) {
  if (method == "legacy_full_edge") {
    return(fusedLassoProximal(
      X = d$X_train,
      Y = d$y_train,
      groups = d$groups_train,
      lambda = cfg$lambda,
      gamma = cfg$gamma,
      G = d$G_l1,
      mu = cfg$mu,
      tol = cfg$tol,
      num.it = cfg$l1_budget,
      c.flag = cfg$c_flag_old,
      intercept = cfg$intercept,
      conserve.memory = cfg$conserve_memory,
      scaling = cfg$scaling
    ))
  }

  solver <- switch(
    method,
    operator = "operator",
    dense_sort = "dense_sort",
    chain_approx = "chain_specialized",
    stop("Unknown L1 method: ", method)
  )

  fusedLassoProximalNew(
    X = d$X_train,
    Y = d$y_train,
    groups = d$groups_train,
    lambda = cfg$lambda,
    gamma = cfg$gamma,
    G = d$G_l1,
    mu = cfg$mu,
    tol = cfg$tol,
    num.it = cfg$l1_budget,
    c.flag = FALSE,
    intercept = cfg$intercept,
    conserve.memory = cfg$conserve_memory,
    scaling = cfg$scaling,
    edge.block = cfg$edge_block,
    solver = solver
  )
}

fit_one <- function(method, d, cfg) {
  run_method_timed_worker(method = method, d = d, cfg = cfg)
}
method_label <- function(method) {
  switch(
    method,
    legacy_full_edge = "Legacy Full-Edge",
    operator = "Operator (Ours)",
    dense_sort = "Dense Sort (Ours)",
    chain_approx = "Chain Approx. (Ours)",
    featureless = "Featureless",
    method
  )
}

worker_cfg_from_parent <- function(cfg) {
  list(
    lambda = cfg$lambda,
    gamma = cfg$gamma,
    mu = cfg$mu,
    tol = cfg$tol,
    l1_budget = cfg$l1_budget,
    c_flag_old = cfg$c_flag_old,
    intercept = cfg$intercept,
    scaling = cfg$scaling,
    conserve_memory = cfg$conserve_memory,
    edge_block = cfg$edge_block
  )
}

parse_time_l_elapsed_sec <- function(lines) {
  if (!length(lines)) return(NA_real_)
  cand <- grep("\\breal\\b", lines, value = TRUE)
  if (!length(cand)) return(NA_real_)
  for (ln in cand) {
    val <- suppressWarnings(as.numeric(sub("^\\s*([0-9]*\\.?[0-9]+)\\s+real.*$", "\\1", ln)))
    if (is.finite(val)) return(val)
  }
  NA_real_
}

parse_time_l_peak_rss_mb <- function(lines) {
  if (!length(lines)) return(NA_real_)
  cand <- grep("maximum resident set size", lines, ignore.case = TRUE, value = TRUE)
  if (!length(cand)) return(NA_real_)

  for (ln in rev(cand)) {
    if (grepl(":", ln, fixed = TRUE)) {
      val <- suppressWarnings(as.numeric(trimws(sub(".*:\\s*([0-9.]+).*", "\\1", ln))))
    } else {
      val <- suppressWarnings(as.numeric(trimws(sub("^\\s*([0-9.]+).*", "\\1", ln))))
    }
    if (is.finite(val)) {
      ln_low <- tolower(ln)
      if (grepl("kbytes|kb", ln_low)) {
        return(val / 1024) # kB -> MB
      }
      # macOS /usr/bin/time -l reports bytes for this field.
      return(val / (1024^2)) # bytes -> MB
    }
  }
  NA_real_
}

run_method_timed_worker <- function(method, d, cfg) {
  payload_path <- tempfile(pattern = "peakrss_payload_", fileext = ".rds")
  out_path <- tempfile(pattern = "peakrss_out_", fileext = ".rds")
  stdout_path <- tempfile(pattern = "peakrss_stdout_", fileext = ".log")
  stderr_path <- tempfile(pattern = "peakrss_stderr_", fileext = ".log")
  on.exit(unlink(c(payload_path, out_path, stdout_path, stderr_path), force = TRUE), add = TRUE)

  payload <- list(
    method = method,
    d = d,
    cfg = worker_cfg_from_parent(cfg)
  )
  saveRDS(payload, payload_path)

  run_args <- c(
    "-l",
    "Rscript",
    "scripts/benchmark_real_l1_cv_fixed_budget.R",
    "--worker_payload", payload_path,
    "--worker_out", out_path
  )
  run_status <- tryCatch(
    system2(
      "/usr/bin/time",
      run_args,
      stdout = stdout_path,
      stderr = stderr_path
    ),
    error = function(e) 1L
  )
  if (is.null(run_status)) run_status <- 0L
  stderr_lines <- if (file.exists(stderr_path)) readLines(stderr_path, warn = FALSE) else character(0)

  out <- if (file.exists(out_path)) readRDS(out_path) else list(ok = FALSE, error = "Missing worker output.")
  if (run_status != 0L || !isTRUE(out$ok)) {
    stop(sprintf(
      "Timed worker failed for method '%s' (status=%s): %s",
      method, as.character(run_status), as.character(out$error %||% "unknown")
    ))
  }

  elapsed_sec <- parse_time_l_elapsed_sec(stderr_lines)
  if (!is.finite(elapsed_sec)) {
    stop(sprintf("Could not parse elapsed wall time from /usr/bin/time output for method '%s'.", method))
  }
  memory_mb <- parse_time_l_peak_rss_mb(stderr_lines)
  if (!is.finite(memory_mb)) {
    warning(sprintf("Could not parse peak RSS from /usr/bin/time output for method '%s'. Setting memory=0.", method))
    memory_mb <- 0
  }

  if (!isTRUE(cfg$measure_memory)) {
    memory_mb <- 0
  }

  list(
    elapsed_sec = as.numeric(elapsed_sec),
    elapsed_ms = as.numeric(elapsed_sec * 1000),
    memory_mb = as.numeric(memory_mb),
    test_rmse = as.numeric(out$test_rmse),
    warning = as.character(out$warning %||% "none")
  )
}

run_fit_worker <- function(argv) {
  payload_path <- argv$worker_payload
  out_path <- argv$worker_out %||% tempfile(pattern = "peakrss_worker_", fileext = ".rds")
  payload <- readRDS(payload_path)
  method <- payload$method
  d <- payload$d
  cfg <- payload$cfg

  warning_msg <- "none"
  ok <- TRUE
  err <- ""

  tryCatch({
    test_rmse_val <- NA_real_
    if (method == "featureless") {
      ybar <- mean(d$y_train)
      test_rmse_val <- rmse(d$y_test, rep(ybar, length(d$y_test)))
    } else {
      fit <- withCallingHandlers(
        {
          fit_l1(method, d, cfg)
        },
        warning = function(w) {
          warning_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      )
      if (is.null(fit)) {
        stop("Fit failed to return coefficients for method: ", method)
      }
      yhat <- predict_from_beta(fit, d$X_test, d$groups_test)
      test_rmse_val <- rmse(d$y_test, yhat)
    }
    out <- list(
      ok = TRUE,
      method = method,
      test_rmse = test_rmse_val,
      warning = warning_msg
    )
    saveRDS(out, out_path)
  }, error = function(e) {
    ok <<- FALSE
    err <<- conditionMessage(e)
    saveRDS(list(ok = FALSE, method = method, error = err), out_path)
  })

  if (!ok) {
    stop("Worker failed: ", err)
  }
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(key) if (key %in% names(argv)) argv[[key]] else DEFAULTS[[key]]

if (!is.null(argv$worker_payload)) {
  run_fit_worker(argv)
  quit(save = "no", status = 0L)
}

cfg <- list(
  data_rds = as.character(get_opt("data_rds")),
  l1_methods = as_chr_vec(as.character(get_opt("l1_methods"))),
  l1_budget = as_int(get_opt("l1_budget")),
  folds = as_int(get_opt("folds")),
  seed = as_int(get_opt("seed")),
  min_group_size = as_int(get_opt("min_group_size")),
  k_keep = as.integer(get_opt("k_keep")),
  g_structure = as.character(get_opt("g_structure")),
  lambda = as.numeric(get_opt("lambda")),
  gamma = as.numeric(get_opt("gamma")),
  mu = as.numeric(get_opt("mu")),
  tol = as.numeric(get_opt("tol")),
  intercept = as_bool(get_opt("intercept"), TRUE),
  scaling = as_bool(get_opt("scaling"), FALSE),
  standardize = as_bool(get_opt("standardize"), TRUE),
  conserve_memory = as_bool(get_opt("conserve_memory"), FALSE),
  measure_memory = as_bool(get_opt("measure_memory"), TRUE),
  edge_block = as_int(get_opt("edge_block")),
  c_flag_old = as_bool(get_opt("c_flag_old"), FALSE),
  time_log_scale = as_bool(get_opt("time_log_scale"), FALSE),
  out_prefix = as.character(get_opt("out_prefix"))
)

check_choices(
  cfg$l1_methods,
  c("legacy_full_edge", "operator", "dense_sort", "chain_approx", "featureless"),
  "l1_methods"
)
check_choices(cfg$g_structure, c("dense", "sparse_chain"), "g_structure")
if (cfg$l1_budget <= 0L) stop("l1_budget must be > 0.")
if (cfg$folds < 2L) stop("folds must be >= 2.")
if (cfg$min_group_size < cfg$folds) {
  stop("min_group_size must be >= folds to ensure each group can appear in each fold.")
}

if (!dir.exists("build")) dir.create("build", recursive = TRUE)

d0 <- prepare_real_data(
  data_rds = cfg$data_rds,
  min_group_size = cfg$min_group_size,
  k_keep = cfg$k_keep
)
G_l1 <- subset_graph_l1(d0, g_structure = cfg$g_structure)

cat(sprintf("Real data: %s\n", cfg$data_rds))
cat(sprintf("n=%d, p=%d, k=%d after filtering\n", nrow(d0$X), ncol(d0$X), length(unique(d0$groups))))
cat(sprintf("L1 methods: %s\n", paste(cfg$l1_methods, collapse = ", ")))
cat(sprintf("Graph structure: %s\n", cfg$g_structure))
cat(sprintf("CV folds: %d (stratified by group)\n", cfg$folds))
cat(sprintf("Fixed L1 budget: %d\n", cfg$l1_budget))
cat(sprintf("Feature standardization: %s\n", ifelse(cfg$standardize, "TRUE (train-fold stats)", "FALSE")))
cat("Timing source: /usr/bin/time -l (worker process wall time)\n")
cat(sprintf(
  "Memory measurement pass: %s\n",
  ifelse(cfg$measure_memory, "TRUE (peak RSS via /usr/bin/time -l in worker process)", "FALSE")
))

rows <- list()
row_id <- 0L

fold_assign <- stratified_folds(d0$groups, nfold = cfg$folds, seed = cfg$seed)
if (length(unique(fold_assign)) < cfg$folds) {
  stop("Failed to assign all requested folds; reduce folds or increase group sizes.")
}

for (fold_id in seq_len(cfg$folds)) {
  test_idx <- which(fold_assign == fold_id)
  train_idx <- which(fold_assign != fold_id)
  if (!length(test_idx) || !length(train_idx)) {
    stop(sprintf("Fold %d has empty train/test split.", fold_id))
  }

  d <- list(
    X_train = d0$X[train_idx, , drop = FALSE],
    y_train = d0$y[train_idx],
    groups_train = d0$groups[train_idx],
    X_test = d0$X[test_idx, , drop = FALSE],
    y_test = d0$y[test_idx],
    groups_test = d0$groups[test_idx],
    G_l1 = G_l1
  )

  if (cfg$standardize) {
    x_center <- colMeans(d$X_train)
    x_scale <- apply(d$X_train, 2L, stats::sd)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
    d$X_train <- scale(d$X_train, center = x_center, scale = x_scale)
    d$X_test <- scale(d$X_test, center = x_center, scale = x_scale)
  }

  for (method in cfg$l1_methods) {
    cat(sprintf("fold=%d method=%s budget=%d ... ", fold_id, method, cfg$l1_budget))
    res <- fit_one(method = method, d = d, cfg = cfg)
    cat(sprintf(
      "time=%.3fms memory=%.3fMB test_rmse=%.6f\n",
      res$elapsed_ms, res$memory_mb, res$test_rmse
    ))

    row_id <- row_id + 1L
    rows[[row_id]] <- data.frame(
      fold = fold_id,
      method = method,
      budget = cfg$l1_budget,
      elapsed_sec = res$elapsed_sec,
      elapsed_ms = res$elapsed_ms,
      memory_mb = res$memory_mb,
      test_rmse = res$test_rmse,
      warning = res$warning,
      stringsAsFactors = FALSE
    )
  }
}

raw_df <- do.call(rbind, rows)
raw_df$method_label <- vapply(raw_df$method, method_label, character(1))

method_levels <- cfg$l1_methods
method_labels <- vapply(method_levels, method_label, character(1))
raw_df$method <- factor(raw_df$method, levels = method_levels)
raw_df$method_label <- factor(raw_df$method_label, levels = method_labels)

mean_df <- aggregate(
  cbind(test_rmse, elapsed_sec, elapsed_ms, memory_mb) ~ method + method_label,
  data = raw_df,
  FUN = mean
)
names(mean_df)[names(mean_df) == "test_rmse"] <- "mean_test_rmse"
names(mean_df)[names(mean_df) == "elapsed_sec"] <- "mean_elapsed_sec"
names(mean_df)[names(mean_df) == "elapsed_ms"] <- "mean_elapsed_ms"
names(mean_df)[names(mean_df) == "memory_mb"] <- "mean_memory_mb"

sd_test <- aggregate(test_rmse ~ method + method_label, data = raw_df, FUN = sd)
sd_time <- aggregate(elapsed_sec ~ method + method_label, data = raw_df, FUN = sd)
sd_time_ms <- aggregate(elapsed_ms ~ method + method_label, data = raw_df, FUN = sd)
sd_mem <- aggregate(memory_mb ~ method + method_label, data = raw_df, FUN = sd)
names(sd_test)[names(sd_test) == "test_rmse"] <- "sd_test_rmse"
names(sd_time)[names(sd_time) == "elapsed_sec"] <- "sd_elapsed_sec"
names(sd_time_ms)[names(sd_time_ms) == "elapsed_ms"] <- "sd_elapsed_ms"
names(sd_mem)[names(sd_mem) == "memory_mb"] <- "sd_memory_mb"

summary_df <- merge(mean_df, sd_test, by = c("method", "method_label"), all.x = TRUE)
summary_df <- merge(summary_df, sd_time, by = c("method", "method_label"), all.x = TRUE)
summary_df <- merge(summary_df, sd_time_ms, by = c("method", "method_label"), all.x = TRUE)
summary_df <- merge(summary_df, sd_mem, by = c("method", "method_label"), all.x = TRUE)
summary_df$sd_test_rmse[is.na(summary_df$sd_test_rmse)] <- 0
summary_df$sd_elapsed_sec[is.na(summary_df$sd_elapsed_sec)] <- 0
summary_df$sd_elapsed_ms[is.na(summary_df$sd_elapsed_ms)] <- 0
summary_df$sd_memory_mb[is.na(summary_df$sd_memory_mb)] <- 0
raw_df$warning_flag <- raw_df$warning != "none"
warn_counts <- aggregate(warning_flag ~ method + method_label, data = raw_df, FUN = sum)
names(warn_counts)[names(warn_counts) == "warning_flag"] <- "n_warnings"

warn_text <- aggregate(
  warning ~ method + method_label,
  data = raw_df,
  FUN = function(x) {
    vals <- unique(x[x != "none"])
    if (!length(vals)) return("none")
    paste(vals, collapse = " | ")
  }
)
names(warn_text)[names(warn_text) == "warning"] <- "warning_types"

summary_df <- merge(summary_df, warn_counts, by = c("method", "method_label"), all.x = TRUE)
summary_df <- merge(summary_df, warn_text, by = c("method", "method_label"), all.x = TRUE)
summary_df$n_warnings[is.na(summary_df$n_warnings)] <- 0L
summary_df$warning_types[is.na(summary_df$warning_types)] <- "none"
summary_df <- summary_df[order(match(as.character(summary_df$method), method_levels)), ]

raw_csv <- sprintf("build/%s_raw.csv", cfg$out_prefix)
summary_csv <- sprintf("build/%s_summary.csv", cfg$out_prefix)
write.csv(raw_df, raw_csv, row.names = FALSE)
write.csv(summary_df, summary_csv, row.names = FALSE)

cat("\nSaved files:\n")
cat(sprintf("- %s\n", raw_csv))
cat(sprintf("- %s\n", summary_csv))
