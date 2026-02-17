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

as_num_vec <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  vals <- as.integer(parts)
  if (any(is.na(vals))) stop("Could not parse --k_values")
  sort(unique(vals))
}

build_G <- function(k, structure) {
  if (structure == "dense") {
    G <- matrix(1, k, k)
    diag(G) <- 0
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
    "Legacy Full-Edge (old_l1)",
    "Full-Edge Operator",
    "Chain-Specialized (Approx)"
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

# -----------------------------
# L1 atime benchmark (old_l1 vs operator vs chain_specialized)
# -----------------------------
l1_obj <- atime::atime(
  N = k_values,
  times = times,
  seconds.limit = seconds_limit,
  N.env.parent = environment(),
  setup = {
    set.seed(seed + 5000L + N)
    k <- N
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
    G <- build_G(k, g_structure)
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
      "Legacy Full-Edge (old_l1)",
      "Full-Edge Operator",
      "Chain-Specialized (Approx)"
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

print_wide_tables(l1_df, "L1-fusion solver: atime (old_l1 vs operator vs chain_specialized)")

cat("\nSaved files:\n")
cat(sprintf("- %s\n", l1_summary_file))
cat(sprintf("- %s\n", l1_rds_file))
cat(sprintf("- %s\n", l1_file))
