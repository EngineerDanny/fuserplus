#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fuserplus)
  library(ggplot2)
})

# Load local development solver variants.
source("R/l1_fusion_new_utils.R")
source("R/l1_fusion_operator_new.R")
source("R/l1_fusion_operator_ws_new.R")
source("R/l1_fusion_dfs_chain.R")
source("R/l1_fusion_chain_specialized.R")
source("R/l1_fusion_new.R")

DEFAULTS <- list(
  data_rds = "data/processed/communities_crime_l2.rds",
  methods = "old_l1,operator,operator_ws,chain_specialized",
  budgets = "20,40,80,120,180,260",
  reps = 1L,
  seed = 20260214L,
  val_frac = 0.20,
  lambda = 1e-3,
  gamma = 1e-3,
  tol = 1e-5,
  intercept = FALSE,
  scaling = FALSE,
  conserve_memory = FALSE,
  edge_block = 256L,
  out_csv = "inst/figures/real_l1_budget_sweep.csv",
  out_png = "inst/figures/real_l1_loss_vs_time.png"
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

as_num_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  vals <- as.numeric(parts)
  if (any(is.na(vals))) stop("Could not parse numeric vector.")
  vals
}

as_chr_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

predict_from_beta <- function(beta, X, groups) {
  group_levels <- sort(unique(groups))
  yhat <- numeric(nrow(X))
  if (is.null(colnames(beta))) {
    colnames(beta) <- as.character(group_levels)
  }
  for (g in group_levels) {
    idx <- which(groups == g)
    g_col <- match(as.character(g), colnames(beta))
    if (is.na(g_col)) {
      g_col <- which(group_levels == g)
    }
    yhat[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta[, g_col, drop = FALSE])
  }
  yhat
}

mse <- function(y, yhat) mean((y - yhat)^2)

stratified_split <- function(groups, val_frac, seed) {
  set.seed(seed)
  val_idx <- integer(0L)
  split_idx <- split(seq_along(groups), groups)
  for (idx in split_idx) {
    n <- length(idx)
    if (n <= 2L) next
    n_val <- floor(n * val_frac)
    n_val <- max(1L, min(n_val, n - 2L))
    val_idx <- c(val_idx, sample(idx, n_val))
  }
  val_idx <- sort(unique(val_idx))
  train_idx <- setdiff(seq_along(groups), val_idx)
  list(train = train_idx, val = val_idx)
}

fit_one <- function(method, d, budget, cfg) {
  warnings_seen <- character(0)
  fit <- NULL
  gc(reset = TRUE)
  tm <- system.time({
    fit <- withCallingHandlers({
      if (method == "old_l1") {
        fusedLassoProximal(
          X = d$X_train,
          Y = d$y_train,
          groups = d$groups_train,
          lambda = cfg$lambda,
          gamma = cfg$gamma,
          G = d$G,
          tol = cfg$tol,
          num.it = budget,
          c.flag = FALSE,
          intercept = cfg$intercept,
          conserve.memory = cfg$conserve_memory,
          scaling = cfg$scaling
        )
      } else if (method %in% c("operator", "operator_ws", "chain_specialized")) {
        solver <- method
        fusedLassoProximalNew(
          X = d$X_train,
          Y = d$y_train,
          groups = d$groups_train,
          lambda = cfg$lambda,
          gamma = cfg$gamma,
          G = d$G,
          tol = cfg$tol,
          num.it = budget,
          c.flag = FALSE,
          intercept = cfg$intercept,
          conserve.memory = cfg$conserve_memory,
          scaling = cfg$scaling,
          edge.block = cfg$edge_block,
          solver = solver
        )
      } else {
        stop("Unknown method: ", method)
      }
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  })

  yhat <- predict_from_beta(fit, d$X_val, d$groups_val)
  list(
    elapsed_sec = unname(tm[["elapsed"]]),
    val_mse = mse(d$y_val, yhat),
    iterations = if (method == "old_l1") {
      suppressWarnings(as.integer(fusedLassoProximalIterationsTaken()))
    } else {
      suppressWarnings(as.integer(fusedLassoProximalNewIterationsTaken()))
    },
    active_edges = if (method == "old_l1") {
      choose(length(unique(d$groups_train)), 2L)
    } else {
      suppressWarnings(as.integer(fusedLassoProximalNewActiveEdges()))
    },
    warnings = if (length(warnings_seen)) paste(unique(warnings_seen), collapse = " | ") else "none"
  )
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(key) if (key %in% names(argv)) argv[[key]] else DEFAULTS[[key]]

cfg <- list(
  data_rds = as.character(get_opt("data_rds")),
  methods = as_chr_vec(as.character(get_opt("methods"))),
  budgets = as.integer(as_num_vec(as.character(get_opt("budgets")))),
  reps = as.integer(get_opt("reps")),
  seed = as.integer(get_opt("seed")),
  val_frac = as.numeric(get_opt("val_frac")),
  lambda = as.numeric(get_opt("lambda")),
  gamma = as.numeric(get_opt("gamma")),
  tol = as.numeric(get_opt("tol")),
  intercept = as_bool(get_opt("intercept"), FALSE),
  scaling = as_bool(get_opt("scaling"), FALSE),
  conserve_memory = as_bool(get_opt("conserve_memory"), FALSE),
  edge_block = as.integer(get_opt("edge_block")),
  out_csv = as.character(get_opt("out_csv")),
  out_png = as.character(get_opt("out_png"))
)

valid_methods <- c("old_l1", "operator", "operator_ws", "chain_specialized")
if (length(cfg$methods) == 0L) cfg$methods <- valid_methods
if (any(!cfg$methods %in% valid_methods)) {
  stop("Invalid method. Allowed: ", paste(valid_methods, collapse = ", "))
}

d0 <- readRDS(cfg$data_rds)
if (!all(c("X", "y", "groups", "G_dense") %in% names(d0))) {
  stop("RDS must contain at least X, y, groups, G_dense.")
}

group_counts <- table(d0$groups)
keep_levels_orig <- as.integer(names(group_counts[group_counts >= 2L]))
if (length(keep_levels_orig) < length(group_counts)) {
  dropped <- length(group_counts) - length(keep_levels_orig)
  cat(sprintf("Dropping %d singleton group(s) for legacy old_l1 compatibility.\n", dropped))
}
keep_rows <- d0$groups %in% keep_levels_orig
X_use <- d0$X[keep_rows, , drop = FALSE]
y_use <- d0$y[keep_rows]
groups_use_orig <- d0$groups[keep_rows]
group_levels_orig <- sort(unique(groups_use_orig))
groups_use <- as.integer(factor(groups_use_orig, levels = group_levels_orig))
G_use <- d0$G_dense[group_levels_orig, group_levels_orig, drop = FALSE]

idx <- stratified_split(groups_use, cfg$val_frac, cfg$seed)
d <- list(
  X_train = X_use[idx$train, , drop = FALSE],
  y_train = y_use[idx$train],
  groups_train = groups_use[idx$train],
  X_val = X_use[idx$val, , drop = FALSE],
  y_val = y_use[idx$val],
  groups_val = groups_use[idx$val],
  G = G_use
)

cat(sprintf("Data split: train=%d val=%d p=%d k=%d\n",
            nrow(d$X_train), nrow(d$X_val), ncol(d$X_train), length(unique(d$groups_train))))
cat(sprintf("Methods: %s\n", paste(cfg$methods, collapse = ", ")))
cat(sprintf("Budgets: %s\n", paste(cfg$budgets, collapse = ", ")))

rows <- list()
row_id <- 0L
for (rep_id in seq_len(cfg$reps)) {
  for (method in cfg$methods) {
    for (budget in cfg$budgets) {
      cat(sprintf("Running rep=%d method=%s budget=%d ... ", rep_id, method, budget))
      res <- fit_one(method, d, budget, cfg)
      cat(sprintf("time=%.3fs val_mse=%.6f\n", res$elapsed_sec, res$val_mse))

      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        rep = rep_id,
        method = method,
        budget = budget,
        elapsed_sec = res$elapsed_sec,
        val_mse = res$val_mse,
        iterations = res$iterations,
        active_edges = res$active_edges,
        warnings = res$warnings,
        stringsAsFactors = FALSE
      )
    }
  }
}

res_df <- do.call(rbind, rows)
mean_df <- aggregate(
  cbind(elapsed_sec, val_mse) ~ method + budget,
  data = res_df,
  FUN = mean
)
mean_df <- mean_df[order(mean_df$method, mean_df$budget), ]

dir.create(dirname(cfg$out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, cfg$out_csv, row.names = FALSE)

method_labels <- c(
  old_l1 = "old_l1",
  operator = "operator",
  operator_ws = "operator_ws",
  chain_specialized = "chain_specialized"
)
mean_df$method_label <- method_labels[mean_df$method]

p <- ggplot(mean_df, aes(x = elapsed_sec, y = val_mse, color = method_label, shape = method_label)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  labs(
    title = "Real Data L1 Budget Sweep (Communities & Crime)",
    subtitle = sprintf("Validation MSE vs wall time | lambda=%.1e gamma=%.1e", cfg$lambda, cfg$gamma),
    x = "Time (seconds)",
    y = "Validation MSE",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 12)

ggsave(cfg$out_png, p, width = 8.0, height = 5.0, dpi = 160)

cat(sprintf("\nSaved raw results: %s\n", cfg$out_csv))
cat(sprintf("Saved figure: %s\n", cfg$out_png))

summary_cols <- mean_df[, c("method", "budget", "elapsed_sec", "val_mse")]
summary_cols <- summary_cols[order(summary_cols$method, summary_cols$budget), ]
cat("\nMean results by method x budget:\n")
print(summary_cols, row.names = FALSE)
