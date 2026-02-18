#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fuserplus)
  library(atime)
})

source("R/l1_fusion_new_utils.R")
source("R/l1_fusion_operator_new.R")
source("R/l1_fusion_dfs_chain.R")
source("R/l1_fusion_chain_specialized.R")

DEFAULTS <- list(
  data_mode = "real",  # real | synthetic
  data_rds = "data/processed/world_bank_wdi/grouped_regression.rds",
  min_group_size = 2L,
  seed = 20260206L,
  k_values = "10,20,30,40,50,80,100",
  p = 120L,
  n_group_train = 35L,
  sigma = 0.05,
  lambda = 1e-3,
  gamma = 1e-3,
  mu = 1e-4,
  tol = 1e-4,
  num_it = 1200L,
  scaling = FALSE,
  intercept = FALSE,
  conserve_memory = FALSE,
  edge_block = 256L,
  c_flag_old = FALSE,
  g_structure = "sparse_chain",    # dense | sparse_chain
  times = 3L,
  seconds_limit = Inf,
  out_prefix = "atime_l1_variants"
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

as_choice <- function(x, choices, arg_name) {
  if (!(x %in% choices)) {
    stop(sprintf("%s must be one of: %s", arg_name, paste(choices, collapse = ", ")))
  }
  x
}

as_num_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  vals <- as.integer(parts)
  if (any(is.na(vals))) stop("Could not parse --k_values")
  sort(unique(vals))
}

build_G <- function(k, structure, dense_diag = 0) {
  if (structure == "dense") {
    G <- matrix(1, k, k)
    diag(G) <- dense_diag
    return(G)
  }
  if (structure == "sparse_chain") {
    G <- matrix(0, k, k)
    if (k > 1L) {
      for (i in 1:(k - 1L)) {
        G[i, i + 1L] <- 1
        G[i + 1L, i] <- 1
      }
    }
    return(G)
  }
  stop("Unknown g_structure: ", structure)
}

prepare_real_pool <- function(data_rds, min_group_size = 2L) {
  d0 <- readRDS(data_rds)
  if (!all(c("X", "y", "groups") %in% names(d0))) {
    stop("Real-data RDS must contain at least X, y, groups.")
  }
  X <- as.matrix(d0$X)
  y <- as.numeric(d0$y)
  groups <- as.integer(d0$groups)
  if (length(groups) != nrow(X) || length(y) != nrow(X)) {
    stop("RDS has inconsistent X/y/groups dimensions.")
  }
  if (stats::sd(y) <= 0) {
    stop("Real-data RDS has constant y; choose a dataset with response variation.")
  }
  counts <- sort(table(groups), decreasing = TRUE)
  counts <- counts[counts >= min_group_size]
  if (!length(counts)) {
    stop("No groups satisfy min_group_size in real dataset.")
  }

  list(
    X = X,
    y = y,
    groups = groups,
    group_order = as.integer(names(counts)),
    G_dense = d0$G_dense %||% NULL,
    G_sparse_chain = d0$G_sparse_chain %||% NULL
  )
}

slice_real_data <- function(k, pool, g_structure) {
  if (k > length(pool$group_order)) {
    stop(sprintf("Requested k=%d but only %d groups available after filtering.", k, length(pool$group_order)))
  }
  keep_levels <- pool$group_order[seq_len(k)]
  keep_rows <- pool$groups %in% keep_levels

  X <- pool$X[keep_rows, , drop = FALSE]
  y <- pool$y[keep_rows]
  if (stats::sd(y) <= 0) {
    stop(sprintf("Selected real-data slice for k=%d has constant y; use another dataset or k range.", k))
  }
  groups_orig <- pool$groups[keep_rows]
  groups <- as.integer(factor(groups_orig, levels = keep_levels))

  if (g_structure == "dense") {
    if (!is.null(pool$G_dense)) {
      G <- as.matrix(pool$G_dense[keep_levels, keep_levels, drop = FALSE])
      diag(G) <- 0
    } else {
      G <- build_G(k, "dense", dense_diag = 0)
    }
  } else if (g_structure == "sparse_chain") {
    if (!is.null(pool$G_sparse_chain)) {
      G <- as.matrix(pool$G_sparse_chain[keep_levels, keep_levels, drop = FALSE])
      diag(G) <- 0
    } else {
      G <- build_G(k, "sparse_chain", dense_diag = 0)
    }
  } else {
    stop("Unknown g_structure: ", g_structure)
  }

  list(X = X, y = y, groups = groups, G = G)
}

summarize_atime <- function(obj) {
  m <- as.data.frame(obj$measurements)
  out <- data.frame(
    k = m$N,
    Method = m$expr.name,
    `Elapsed time (seconds)` = m$median,
    `Memory (MB)` = m$kilobytes / 1024,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  method_levels <- c(
    "Legacy Full-Edge",
    "Operator (Ours)",
    "Chain Approx. (Ours)"
  )
  out$Method <- factor(out$Method, levels = method_levels)
  out <- out[order(out$k, out$Method), , drop = FALSE]
  out$Method <- as.character(out$Method)
  rownames(out) <- NULL
  out
}

print_wide_tables <- function(df, header) {
  methods <- unique(df$Method)
  k_vals <- sort(unique(df$k))

  elapsed <- data.frame(k = k_vals)
  memory <- data.frame(k = k_vals)
  for (m in methods) {
    idx <- df$Method == m
    elapsed[[m]] <- df$`Elapsed time (seconds)`[idx][match(k_vals, df$k[idx])]
    memory[[m]] <- df$`Memory (MB)`[idx][match(k_vals, df$k[idx])]
  }

  cat("\n", header, "\n", sep = "")
  cat("Elapsed time (seconds):\n")
  print(elapsed, row.names = FALSE)
  cat("\nMemory (MB):\n")
  print(memory, row.names = FALSE)
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(key) if (key %in% names(argv)) argv[[key]] else DEFAULTS[[key]]

`%||%` <- function(x, y) if (is.null(x)) y else x

data_mode <- as_choice(as.character(get_opt("data_mode")), c("real", "synthetic"), "--data_mode")
data_rds <- as.character(get_opt("data_rds"))
min_group_size <- as.integer(get_opt("min_group_size"))
seed <- as.integer(get_opt("seed"))
k_values <- as_num_vec(get_opt("k_values"))
p <- as.integer(get_opt("p"))
n_group_train <- as.integer(get_opt("n_group_train"))
sigma <- as.numeric(get_opt("sigma"))
lambda <- as.numeric(get_opt("lambda"))
gamma <- as.numeric(get_opt("gamma"))
mu <- as.numeric(get_opt("mu"))
tol <- as.numeric(get_opt("tol"))
num_it <- as.integer(get_opt("num_it"))
scaling <- as_bool(get_opt("scaling"), FALSE)
intercept <- as_bool(get_opt("intercept"), FALSE)
conserve_memory <- as_bool(get_opt("conserve_memory"), FALSE)
edge_block <- as.integer(get_opt("edge_block"))
c_flag_old <- as_bool(get_opt("c_flag_old"), FALSE)
g_structure <- as.character(get_opt("g_structure"))
times <- as.integer(get_opt("times"))
seconds_limit <- as.numeric(get_opt("seconds_limit"))
out_prefix <- as.character(get_opt("out_prefix"))

if (!dir.exists("inst/figures")) dir.create("inst/figures", recursive = TRUE)
if (!dir.exists("build")) dir.create("build", recursive = TRUE)

real_pool <- NULL
if (data_mode == "real") {
  real_pool <- prepare_real_pool(data_rds = data_rds, min_group_size = min_group_size)
  max_k <- length(real_pool$group_order)
  if (max(k_values) > max_k) {
    stop(sprintf("Requested max k=%d exceeds available k=%d from %s.", max(k_values), max_k, data_rds))
  }
  cat(sprintf("Using real data: %s\n", data_rds))
  cat(sprintf("Available groups (>= %d rows): %d\n", min_group_size, max_k))
} else {
  cat("Using synthetic data mode.\n")
}

# -----------------------------
# L1 atime benchmark (old_l1 vs operator vs chain_specialized)
# -----------------------------
l1_obj <- atime::atime(
  N = k_values,
  times = times,
  seconds.limit = seconds_limit,
  N.env.parent = environment(),
  setup = {
    k <- N
    if (data_mode == "real") {
      d <- slice_real_data(k = k, pool = real_pool, g_structure = g_structure)
      X <- d$X
      y <- d$y
      groups <- d$groups
      G <- d$G
    } else {
      set.seed(seed + 5000L + N)
      groups <- rep(seq_len(k), each = n_group_train)
      n <- length(groups)
      X <- matrix(rnorm(n * p), nrow = n, ncol = p)

      beta <- matrix(0, nrow = p, ncol = k)
      nonzero <- rbinom(p * k, 1, 0.02 / max(1L, k))
      if (sum(nonzero) > 0L) {
        beta[which(nonzero == 1L)] <- rnorm(sum(nonzero), 0.8, 0.25)
      }
      shared <- rbinom(p, 1, 0.02)
      if (sum(shared) > 0L) {
        beta[which(shared == 1L), ] <- rnorm(sum(shared), -0.8, 0.25)
      }

      y <- numeric(n)
      for (g in seq_len(k)) {
        idx <- which(groups == g)
        y[idx] <- X[idx, , drop = FALSE] %*% beta[, g] + rnorm(length(idx), sd = sigma)
      }
      G <- build_G(k, g_structure, dense_diag = 0)
    }
  },
  expr.list = setNames(
    list(
      quote(
        fusedLassoProximal(
          X, y, groups,
          lambda = lambda, gamma = gamma, G = G,
          mu = mu, tol = tol, num.it = num_it,
          c.flag = c_flag_old,
          intercept = intercept,
          conserve.memory = conserve_memory,
          scaling = scaling
        )
      ),
      quote(
        fusedLassoProximalNewOperator(
          X, y, groups,
          lambda = lambda, gamma = gamma, G = G,
          mu = mu, tol = tol, num.it = num_it,
          c.flag = FALSE,
          intercept = intercept,
          conserve.memory = conserve_memory,
          scaling = scaling,
          edge.block = edge_block
        )
      ),
      quote(
        fusedLassoProximalChainSpecialized(
          X, y, groups,
          lambda = lambda, gamma = gamma, G = G,
          mu = mu, tol = tol, num.it = num_it,
          c.flag = FALSE,
          intercept = intercept,
          conserve.memory = conserve_memory,
          scaling = scaling,
          edge.block = edge_block
        )
      )
    ),
    c(
      "Legacy Full-Edge",
      "Operator (Ours)",
      "Chain Approx. (Ours)"
    )
  )
)

l1_df <- summarize_atime(l1_obj)
l1_summary_file <- sprintf("build/%s_summary.csv", out_prefix)
write.csv(l1_df, l1_summary_file, row.names = FALSE)
l1_rds_file <- sprintf("build/%s_atime_obj.rds", out_prefix)
saveRDS(l1_obj, l1_rds_file)

l1_plot <- plot(l1_obj) +
  ggplot2::labs(x = "k (number of subsets)")
l1_file <- sprintf("inst/figures/%s_atime.png", out_prefix)
ggplot2::ggsave(l1_file, l1_plot, width = 5.0, height = 3.4, dpi = 500)

print_wide_tables(l1_df, "L1-fusion solver: atime (Legacy Full-Edge vs Operator vs Chain Approx.)")

cat("\nSaved files:\n")
cat(sprintf("- %s\n", l1_summary_file))
cat(sprintf("- %s\n", l1_rds_file))
cat(sprintf("- %s\n", l1_file))
