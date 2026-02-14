#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fuserplus)
  library(atime)
  library(grid)
})

source("R/l1_fusion_new_utils.R")
source("R/l1_fusion_operator_new.R")
source("R/l1_fusion_new.R")
source("R/l2_fusion_new.R")

DEFAULTS <- list(
  seed = 20260214L,
  k_values = "20,40,60,80",
  p = 120L,
  n_group_train = 30L,
  sigma = 0.05,
  lambda = 1e-3,
  gamma = 1e-3,
  tol = 1e-4,
  num_it = 400L,
  scaling = FALSE,
  g_structure = "dense",           # dense | sparse_chain
  times = 3L,
  seconds_limit = 2.0,
  edge_block = 256L,
  c_flag_old = FALSE,
  out_prefix = "atime_l1_l2"
)

parse_args <- function(argv) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    key <- substring(key, 3)
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
  vals <- as.integer(parts)
  if (any(is.na(vals))) stop("Could not parse --k_values")
  sort(unique(vals))
}

build_G <- function(k, structure) {
  if (structure == "dense") {
    return(matrix(1, k, k))
  }
  if (structure == "sparse_chain") {
    G <- matrix(0, k, k)
    if (k > 1) {
      for (i in 1:(k - 1)) {
        G[i, i + 1] <- 1
        G[i + 1, i] <- 1
      }
    }
    return(G)
  }
  stop("Unknown g_structure: ", structure)
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
  out[order(out$k, out$Method), ]
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

seed <- as.integer(get_opt("seed"))
k_values <- as_num_vec(get_opt("k_values"))
p <- as.integer(get_opt("p"))
n_group_train <- as.integer(get_opt("n_group_train"))
sigma <- as.numeric(get_opt("sigma"))
lambda <- as.numeric(get_opt("lambda"))
gamma <- as.numeric(get_opt("gamma"))
tol <- as.numeric(get_opt("tol"))
num_it <- as.integer(get_opt("num_it"))
scaling <- as_bool(get_opt("scaling"), FALSE)
g_structure <- as.character(get_opt("g_structure"))
times <- as.integer(get_opt("times"))
seconds_limit <- as.numeric(get_opt("seconds_limit"))
edge_block <- as.integer(get_opt("edge_block"))
c_flag_old <- as_bool(get_opt("c_flag_old"), TRUE)
out_prefix <- as.character(get_opt("out_prefix"))

if (!dir.exists("inst/figures")) dir.create("inst/figures", recursive = TRUE)

# -----------------------------
# L1 atime benchmark
# -----------------------------
l1_obj <- atime::atime(
  N = k_values,
  times = times,
  seconds.limit = seconds_limit,
  N.env.parent = environment(),
  setup = {
    set.seed(seed + N)
    k <- N
    groups <- rep(seq_len(k), each = n_group_train)
    n <- length(groups)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    beta <- matrix(rnorm(p * k, sd = 0.15), p, k)
    y <- numeric(n)
    for (g in seq_len(k)) {
      idx <- which(groups == g)
      y[idx] <- X[idx, , drop = FALSE] %*% beta[, g] + rnorm(length(idx), sd = sigma)
    }
    G <- build_G(k, g_structure)
  },
  expr.list = list(
    old = quote(
      fusedLassoProximal(
        X, y, groups,
        lambda = lambda, gamma = gamma, G = G,
        tol = tol, num.it = num_it,
        intercept = FALSE, scaling = scaling,
        c.flag = c_flag_old
      )
    ),
    new = quote(
      fusedLassoProximalNewOperator(
        X, y, groups,
        lambda = lambda, gamma = gamma, G = G,
        tol = tol, num.it = num_it,
        intercept = FALSE, scaling = scaling,
        edge.block = edge_block
      )
    )
  )
)

l1_df <- summarize_atime(l1_obj)
write.csv(l1_df, sprintf("scripts/%s_l1_summary.csv", out_prefix), row.names = FALSE)
l1_plot <- plot(l1_obj) + ggplot2::ggtitle("A. L1-fusion solver")

# -----------------------------
# L2 atime benchmark
# -----------------------------
l2_obj <- atime::atime(
  N = k_values,
  times = times,
  seconds.limit = seconds_limit,
  N.env.parent = environment(),
  setup = {
    set.seed(seed + 10000 + N)
    k <- N
    groups <- rep(seq_len(k), each = n_group_train)
    n <- length(groups)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    beta <- matrix(rnorm(p * k, sd = 0.15), p, k)
    y <- numeric(n)
    for (g in seq_len(k)) {
      idx <- which(groups == g)
      y[idx] <- X[idx, , drop = FALSE] %*% beta[, g] + rnorm(length(idx), sd = sigma)
    }
    G <- build_G(k, g_structure)
  },
  expr.list = list(
    old = quote(
      fusedL2DescentGLMNet(
        X, y, groups,
        lambda = lambda, gamma = gamma, G = G, scaling = scaling
      )
    ),
    new = quote(
      fusedL2DescentGLMNetNew(
        X, y, groups,
        lambda = lambda, gamma = gamma, G = G, scaling = scaling
      )
    )
  )
)

l2_df <- summarize_atime(l2_obj)
write.csv(l2_df, sprintf("scripts/%s_l2_summary.csv", out_prefix), row.names = FALSE)
l2_plot <- plot(l2_obj) + ggplot2::ggtitle("B. L2-fusion solver")

combined_file <- sprintf("inst/figures/%s_combined_atime.png", out_prefix)
png(filename = combined_file, width = 1224, height = 630, res = 180)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
print(l1_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(l2_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Keep only combined atime figure for this prefix.
unlink(sprintf("inst/figures/%s_l1_atime.png", out_prefix), force = TRUE)
unlink(sprintf("inst/figures/%s_l2_atime.png", out_prefix), force = TRUE)

print_wide_tables(l1_df, "L1-fusion solver: atime old vs new")
print_wide_tables(l2_df, "L2-fusion solver: atime old vs new")

cat("\nSaved files:\n")
cat(sprintf("- scripts/%s_l1_summary.csv\n", out_prefix))
cat(sprintf("- scripts/%s_l2_summary.csv\n", out_prefix))
cat(sprintf("- inst/figures/%s_combined_atime.png\n", out_prefix))
