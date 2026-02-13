#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fuserplus))
source("R/l1_fusion.R")
source("R/l1_fusion_new_utils.R")
source("R/l1_fusion_operator_new.R")
source("R/l1_fusion_dense_sort.R")
source("R/l1_fusion_dfs_chain.R")
source("R/l1_fusion_chain_specialized.R")
source("R/l1_fusion_new.R")

DEFAULTS <- list(
  mode = "sweep",                          # sweep | run
  method = NULL,                           # required in mode=run
  methods = "old_l1,operator,dense_sort_scaffold,dfs_chain,chain_specialized",
  k = NULL,                                # required in mode=run
  k_values = "40,80,120,160",
  reps = 2L,
  seed = 20260214L,
  p = 300L,
  n_group_train = 20L,
  n_group_test = 10L,
  sigma = 0.1,
  lambda = 1e-3,
  gamma = 2.0,
  tol = 1e-20,
  num_it = 700L,
  intercept = FALSE,
  scaling = FALSE,
  conserve_memory = FALSE,
  c_flag_old = FALSE,
  c_flag_new = FALSE,
  graph_mode = "dense_uniform",            # dense_uniform | dense_nonuniform | sparse
  g_offdiag = 1.0,
  g_diag = 1.0,
  weight_min = 0.2,
  weight_max = 2.0,
  sparse_prob = 0.1,
  sparse_ensure_chain = TRUE,
  require_uniform_weights = TRUE,
  fallback_dense = TRUE,
  chain_use_mst = TRUE,
  chain_start = 1L,
  chain_min_weight = 1e-8,
  data_model = "shared_plus_group_noise",  # shared_plus_group_noise | sparse_mixed
  beta_base_sd = 0.04,
  beta_group_sd = 0.006,
  sparse_shared_prob = 0.02,
  sparse_group_prob_scale = 0.02,
  sparse_group_mean = 1.0,
  sparse_group_sd = 0.25,
  sparse_shared_mean = -1.0,
  sparse_shared_sd = 0.25,
  print_spec = TRUE
)

parse_args <- function(argv) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
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
  if (is.null(x) || nchar(x) == 0L) return(numeric(0))
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) {
    stop("Could not parse numeric vector.")
  }
  vals
}

as_chr_vec <- function(x) {
  if (is.null(x) || nchar(x) == 0L) return(character(0))
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts[nzchar(parts)]
}

predict_from_beta <- function(beta, X, groups) {
  group_levels <- sort(unique(groups))
  if (is.null(colnames(beta))) colnames(beta) <- as.character(group_levels)
  yhat <- numeric(nrow(X))
  for (g in group_levels) {
    idx <- which(groups == g)
    g_col <- match(as.character(g), colnames(beta))
    if (is.na(g_col)) g_col <- which(group_levels == g)
    yhat[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta[, g_col, drop = FALSE])
  }
  yhat
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

build_G <- function(k, cfg) {
  mode <- cfg$graph_mode
  diag_val <- cfg$g_diag

  if (mode == "dense_uniform") {
    G <- matrix(cfg$g_offdiag, k, k)
    diag(G) <- diag_val
    return(G)
  }

  if (mode == "dense_nonuniform") {
    if (cfg$weight_min <= 0 || cfg$weight_max < cfg$weight_min) {
      stop("For dense_nonuniform, require 0 < weight_min <= weight_max.")
    }
    U <- matrix(runif(k * k, min = cfg$weight_min, max = cfg$weight_max), k, k)
    G <- (U + t(U)) / 2
    diag(G) <- diag_val
    return(G)
  }

  if (mode == "sparse") {
    if (cfg$sparse_prob < 0 || cfg$sparse_prob > 1) {
      stop("For sparse graph_mode, sparse_prob must be in [0,1].")
    }
    if (cfg$weight_min <= 0 || cfg$weight_max < cfg$weight_min) {
      stop("For sparse graph_mode, require 0 < weight_min <= weight_max.")
    }

    G <- matrix(0, k, k)
    if (k > 1L) {
      for (i in seq_len(k - 1L)) {
        for (j in (i + 1L):k) {
          if (runif(1) <= cfg$sparse_prob) {
            w <- runif(1, min = cfg$weight_min, max = cfg$weight_max)
            G[i, j] <- w
            G[j, i] <- w
          }
        }
      }
      if (isTRUE(cfg$sparse_ensure_chain)) {
        for (i in seq_len(k - 1L)) {
          if (G[i, i + 1L] <= 0) {
            w <- runif(1, min = cfg$weight_min, max = cfg$weight_max)
            G[i, i + 1L] <- w
            G[i + 1L, i] <- w
          }
        }
      }
    }
    diag(G) <- diag_val
    return(G)
  }

  stop("Unknown graph_mode: ", mode)
}

generate_beta <- function(k, p, data_model, cfg) {
  if (data_model == "shared_plus_group_noise") {
    base_beta <- matrix(rnorm(p, sd = cfg$beta_base_sd), p, 1L)
    return(base_beta %*% matrix(1, 1, k) + matrix(rnorm(p * k, sd = cfg$beta_group_sd), p, k))
  }

  if (data_model == "sparse_mixed") {
    beta_true <- matrix(0, nrow = p, ncol = k)
    nz.group <- rbinom(p * k, 1, cfg$sparse_group_prob_scale / k)
    nz.shared <- rbinom(p, 1, cfg$sparse_shared_prob)
    beta_true[which(nz.group == 1)] <- rnorm(sum(nz.group), cfg$sparse_group_mean, cfg$sparse_group_sd)
    beta_true[which(nz.shared == 1), ] <- rnorm(sum(nz.shared), cfg$sparse_shared_mean, cfg$sparse_shared_sd)
    return(beta_true)
  }

  stop("Unknown data_model: ", data_model)
}

generate_data <- function(seed, k, p, n_group_train, n_group_test, sigma, cfg) {
  set.seed(seed)
  groups_train <- rep(seq_len(k), each = n_group_train)
  groups_test <- rep(seq_len(k), each = n_group_test)
  G <- build_G(k, cfg = cfg)
  beta_true <- generate_beta(k, p, cfg$data_model, cfg)

  X_train <- matrix(rnorm(length(groups_train) * p), nrow = length(groups_train), ncol = p)
  X_test <- matrix(rnorm(length(groups_test) * p), nrow = length(groups_test), ncol = p)

  gen_y <- function(X, groups) {
    y <- numeric(nrow(X))
    for (g in seq_len(k)) {
      idx <- which(groups == g)
      y[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta_true[, g, drop = FALSE]) +
        rnorm(length(idx), sd = sigma)
    }
    y
  }

  list(
    X_train = X_train,
    y_train = gen_y(X_train, groups_train),
    groups_train = groups_train,
    X_test = X_test,
    y_test = gen_y(X_test, groups_test),
    groups_test = groups_test,
    G = G
  )
}

fit_once <- function(method, d, cfg) {
  warnings_seen <- character(0)
  fit <- NULL
  gc(reset = TRUE)
  t <- system.time({
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
          num.it = cfg$num_it,
          c.flag = cfg$c_flag_old,
          intercept = cfg$intercept,
          conserve.memory = cfg$conserve_memory,
          scaling = cfg$scaling
        )
      } else if (method %in% c("operator", "dense_sort_scaffold", "dfs_chain", "chain_specialized")) {
        solver <- switch(
          method,
          operator = "operator",
          dense_sort_scaffold = "dense_sort",
          dfs_chain = "dfs_chain",
          chain_specialized = "chain_specialized"
        )
        fusedLassoProximalNew(
          X = d$X_train,
          Y = d$y_train,
          groups = d$groups_train,
          lambda = cfg$lambda,
          gamma = cfg$gamma,
          G = d$G,
          tol = cfg$tol,
          num.it = cfg$num_it,
          c.flag = cfg$c_flag_new,
          intercept = cfg$intercept,
          conserve.memory = cfg$conserve_memory,
          scaling = cfg$scaling,
          solver = solver,
          chain.use.mst = cfg$chain_use_mst,
          chain.start = cfg$chain_start,
          chain.min.weight = cfg$chain_min_weight,
          require_uniform_weights = cfg$require_uniform_weights,
          fallback = cfg$fallback_dense
        )
      } else {
        stop("Unknown method: ", method)
      }
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  })

  iter <- NA_integer_
  if (method == "old_l1" && exists("fusedLassoProximalIterationsTaken", mode = "function")) {
    iter <- suppressWarnings(as.integer(fusedLassoProximalIterationsTaken()))
  }
  if (method != "old_l1") {
    iter <- suppressWarnings(as.integer(fusedLassoProximalNewIterationsTaken()))
  }

  active_edges <- if (method == "old_l1") {
    choose(length(unique(d$groups_train)), 2L)
  } else {
    suppressWarnings(as.integer(fusedLassoProximalNewActiveEdges()))
  }

  yhat <- predict_from_beta(fit, d$X_test, d$groups_test)
  gc_after <- gc()
  memory_mb <- sum(gc_after[, 7], na.rm = TRUE)
  list(
    method = method,
    elapsed = unname(t[["elapsed"]]),
    memory_mb = memory_mb,
    test_rmse = rmse(d$y_test, yhat),
    iterations = iter,
    active_edges = active_edges,
    warnings = if (length(warnings_seen)) paste(unique(warnings_seen), collapse = " | ") else "none"
  )
}

print_result <- function(x, k, p, rep_id) {
  cat(sprintf(
    paste0(
      "RESULT method=%s k=%d p=%d rep=%d ",
      "Elapsed time (seconds)=%.6f Memory (MB)=%.2f Test RMSE=%.6f ",
      "iterations=%s active_edges=%s warnings=%s\n"
    ),
    x$method, k, p, rep_id, x$elapsed, x$memory_mb, x$test_rmse,
    as.character(x$iterations), as.character(x$active_edges), shQuote(x$warnings)
  ))
}

print_spec <- function(cfg, methods, k_vals) {
  cat("Synthetic Data Spec:\n")
  cat(sprintf("  data_model=%s\n", cfg$data_model))
  if (cfg$data_model == "shared_plus_group_noise") {
    cat(sprintf("  beta[:, g] = base_beta + noise_g, base_sd=%.6f, group_sd=%.6f\n", cfg$beta_base_sd, cfg$beta_group_sd))
  } else {
    cat(sprintf(
      "  sparse_mixed: shared_prob=%.6f, group_prob_scale=%.6f, shared=(mean=%.3f,sd=%.3f), group=(mean=%.3f,sd=%.3f)\n",
      cfg$sparse_shared_prob, cfg$sparse_group_prob_scale,
      cfg$sparse_shared_mean, cfg$sparse_shared_sd,
      cfg$sparse_group_mean, cfg$sparse_group_sd
    ))
  }
  cat(sprintf("  X ~ N(0,1), y = X beta_g + eps, eps_sd=%.6f\n", cfg$sigma))
  if (cfg$graph_mode == "dense_uniform") {
    cat(sprintf("  G mode=dense_uniform: offdiag=%.6f, diag=%.6f\n", cfg$g_offdiag, cfg$g_diag))
  } else if (cfg$graph_mode == "dense_nonuniform") {
    cat(sprintf("  G mode=dense_nonuniform: weights~U(%.3f, %.3f), diag=%.6f\n",
                cfg$weight_min, cfg$weight_max, cfg$g_diag))
  } else if (cfg$graph_mode == "sparse") {
    cat(sprintf("  G mode=sparse: edge_prob=%.3f, weights~U(%.3f, %.3f), ensure_chain=%s, diag=%.6f\n",
                cfg$sparse_prob, cfg$weight_min, cfg$weight_max,
                ifelse(cfg$sparse_ensure_chain, "TRUE", "FALSE"), cfg$g_diag))
  }
  cat(sprintf("  methods=%s\n", paste(methods, collapse = ",")))
  cat(sprintf("  k_values=%s\n", paste(k_vals, collapse = ",")))
  cat(sprintf("  p=%d, n_group_train=%d, n_group_test=%d, reps=%d\n\n",
              cfg$p, cfg$n_group_train, cfg$n_group_test, cfg$reps))
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
get_arg <- function(x, key) if (key %in% names(x)) x[[key]] else NULL
get_opt <- function(key) {
  cli <- get_arg(argv, key)
  if (!is.null(cli)) return(cli)
  DEFAULTS[[key]]
}

cfg <- list(
  mode = as.character(get_opt("mode")),
  method = as.character(get_opt("method")),
  methods = as_chr_vec(as.character(get_opt("methods"))),
  k_run = get_opt("k"),
  k_values = as.integer(as_num_vec(get_opt("k_values"))),
  reps = as.integer(get_opt("reps")),
  seed = as.integer(get_opt("seed")),
  p = as.integer(get_opt("p")),
  n_group_train = as.integer(get_opt("n_group_train")),
  n_group_test = as.integer(get_opt("n_group_test")),
  sigma = as.numeric(get_opt("sigma")),
  lambda = as.numeric(get_opt("lambda")),
  gamma = as.numeric(get_opt("gamma")),
  tol = as.numeric(get_opt("tol")),
  num_it = as.integer(get_opt("num_it")),
  intercept = as_bool(get_opt("intercept"), default = FALSE),
  scaling = as_bool(get_opt("scaling"), default = FALSE),
  conserve_memory = as_bool(get_opt("conserve_memory"), default = FALSE),
  c_flag_old = as_bool(get_opt("c_flag_old"), default = FALSE),
  c_flag_new = as_bool(get_opt("c_flag_new"), default = FALSE),
  graph_mode = as.character(get_opt("graph_mode")),
  g_offdiag = as.numeric(get_opt("g_offdiag")),
  g_diag = as.numeric(get_opt("g_diag")),
  weight_min = as.numeric(get_opt("weight_min")),
  weight_max = as.numeric(get_opt("weight_max")),
  sparse_prob = as.numeric(get_opt("sparse_prob")),
  sparse_ensure_chain = as_bool(get_opt("sparse_ensure_chain"), default = TRUE),
  require_uniform_weights = as_bool(get_opt("require_uniform_weights"), default = TRUE),
  fallback_dense = as_bool(get_opt("fallback_dense"), default = TRUE),
  chain_use_mst = as_bool(get_opt("chain_use_mst"), default = TRUE),
  chain_start = as.integer(get_opt("chain_start")),
  chain_min_weight = as.numeric(get_opt("chain_min_weight")),
  data_model = as.character(get_opt("data_model")),
  beta_base_sd = as.numeric(get_opt("beta_base_sd")),
  beta_group_sd = as.numeric(get_opt("beta_group_sd")),
  sparse_shared_prob = as.numeric(get_opt("sparse_shared_prob")),
  sparse_group_prob_scale = as.numeric(get_opt("sparse_group_prob_scale")),
  sparse_group_mean = as.numeric(get_opt("sparse_group_mean")),
  sparse_group_sd = as.numeric(get_opt("sparse_group_sd")),
  sparse_shared_mean = as.numeric(get_opt("sparse_shared_mean")),
  sparse_shared_sd = as.numeric(get_opt("sparse_shared_sd")),
  print_spec = as_bool(get_opt("print_spec"), default = TRUE)
)

valid_methods <- c("old_l1", "operator", "dense_sort_scaffold", "dfs_chain", "chain_specialized")
if (!length(cfg$methods)) cfg$methods <- valid_methods
if (any(!cfg$methods %in% valid_methods)) {
  stop("Invalid method in --methods. Allowed: ", paste(valid_methods, collapse = ", "))
}

if (cfg$print_spec) {
  spec_k <- if (cfg$mode == "run" && !is.null(cfg$k_run) && cfg$k_run != "") as.integer(cfg$k_run) else cfg$k_values
  print_spec(cfg, cfg$methods, spec_k)
}

if (cfg$mode == "run") {
  if (is.null(cfg$method) || cfg$method == "NULL" || cfg$method == "") {
    stop("--method is required in mode=run")
  }
  if (!(cfg$method %in% valid_methods)) {
    stop("--method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  if (is.null(cfg$k_run) || cfg$k_run == "NULL" || cfg$k_run == "") {
    stop("--k is required in mode=run")
  }
  k <- as.integer(cfg$k_run)
  d <- generate_data(cfg$seed, k, cfg$p, cfg$n_group_train, cfg$n_group_test, cfg$sigma, cfg)
  res <- fit_once(cfg$method, d, cfg)
  print_result(res, k = k, p = cfg$p, rep_id = 1L)
  quit(save = "no", status = 0)
}

if (cfg$mode != "sweep") {
  stop("Unknown mode: ", cfg$mode)
}
if (!length(cfg$k_values)) {
  stop("k_values must be non-empty in mode=sweep")
}

rows <- list()
idx <- 1L
for (k in cfg$k_values) {
  for (m in cfg$methods) {
    for (r in seq_len(cfg$reps)) {
      d <- generate_data(
        seed = cfg$seed + 1000L * r + k,
        k = k,
        p = cfg$p,
        n_group_train = cfg$n_group_train,
        n_group_test = cfg$n_group_test,
        sigma = cfg$sigma,
        cfg = cfg
      )
      res <- fit_once(m, d, cfg)
      print_result(res, k = k, p = cfg$p, rep_id = r)
      rows[[idx]] <- data.frame(
        method = res$method,
        k = k,
        p = cfg$p,
        rep = r,
        elapsed = res$elapsed,
        memory_mb = res$memory_mb,
        test_rmse = res$test_rmse,
        iterations = res$iterations,
        active_edges = res$active_edges,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
}

run_df <- do.call(rbind, rows)
mean_df <- aggregate(cbind(elapsed, memory_mb, test_rmse, iterations, active_edges) ~ method + k + p, run_df, mean)
cat("\nMean benchmark metrics by method and k:\n")
summary_df <- mean_df[, c("method", "k", "elapsed", "memory_mb", "test_rmse", "iterations", "active_edges")]
names(summary_df) <- c(
  "Method", "k", "Elapsed time (seconds)", "Memory (MB)", "Test RMSE",
  "Iterations", "Active Edges"
)
print(summary_df, row.names = FALSE)

expo_rows <- list()
idx <- 1L
for (m in unique(run_df$method)) {
  sub <- mean_df[mean_df$method == m, , drop = FALSE]
  if (nrow(sub) >= 2L && all(sub$elapsed > 0) && all(sub$k > 0)) {
    fit_t <- lm(log(elapsed) ~ log(k), data = sub)
    fit_rmse <- lm(log(test_rmse) ~ log(k), data = sub)
    expo_rows[[idx]] <- data.frame(
      method = m,
      time_exponent_vs_k = unname(coef(fit_t)[["log(k)"]]),
      rmse_exponent_vs_k = unname(coef(fit_rmse)[["log(k)"]]),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
if (length(expo_rows)) {
  expo_df <- do.call(rbind, expo_rows)
  cat("\nEstimated log-log exponents vs k:\n")
  print(expo_df, row.names = FALSE)
}
