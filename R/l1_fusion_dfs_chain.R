# DFS-chain approximation for L1 fused regression.
# Approximates a general graph by a DFS-induced chain and solves
# the resulting problem with the operator-based solver.

.dfs_chain_max_spanning_tree <- function(G) {
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("Parameter 'G' must be a square matrix.")
  }

  k <- nrow(G)
  if (k <= 1L) {
    return(matrix(0, nrow = k, ncol = k))
  }

  W <- G
  W[is.na(W)] <- 0
  diag(W) <- 0
  W <- 0.5 * (W + t(W))

  in_tree <- rep(FALSE, k)
  in_tree[1L] <- TRUE
  tree <- matrix(0, nrow = k, ncol = k)

  for (step in 2L:k) {
    best_w <- -Inf
    best_u <- NA_integer_
    best_v <- NA_integer_
    u_set <- which(in_tree)
    v_set <- which(!in_tree)

    for (u in u_set) {
      w_row <- W[u, v_set]
      idx <- which.max(w_row)
      w_val <- w_row[idx]
      if (length(w_val) && is.finite(w_val) && (w_val > best_w)) {
        best_w <- w_val
        best_u <- u
        best_v <- v_set[idx]
      }
    }

    # If graph is disconnected/non-positive, connect arbitrarily to keep a chainable tree.
    if (is.na(best_u) || is.na(best_v)) {
      best_u <- u_set[1L]
      best_v <- v_set[1L]
    }

    tree[best_u, best_v] <- 1
    tree[best_v, best_u] <- 1
    in_tree[best_v] <- TRUE
  }

  tree
}

.dfs_chain_order_from_adj <- function(adj, start = 1L) {
  k <- nrow(adj)
  start <- as.integer(start)
  if (!is.finite(start) || start < 1L || start > k) {
    start <- 1L
  }

  visited <- rep(FALSE, k)
  order <- integer(0L)
  stack <- c(start)

  while (length(stack)) {
    u <- stack[length(stack)]
    stack <- stack[-length(stack)]
    if (visited[u]) {
      next
    }
    visited[u] <- TRUE
    order <- c(order, u)

    nbrs <- which(adj[u, ] > 0)
    nbrs <- nbrs[!visited[nbrs]]
    if (length(nbrs)) {
      # Push in reverse so traversal is deterministic and ascending.
      nbrs <- sort(nbrs, decreasing = TRUE)
      stack <- c(stack, nbrs)
    }
  }

  if (length(order) < k) {
    order <- c(order, which(!visited))
  }
  order
}

.dfs_chain_graph <- function(G, use_mst = TRUE, start = 1L, min_weight = 1e-8) {
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("Parameter 'G' must be a square matrix.")
  }
  k <- nrow(G)
  if (k <= 1L) {
    return(list(order = seq_len(k), G_chain = G))
  }

  if (!is.numeric(min_weight) || length(min_weight) != 1L || !is.finite(min_weight) || min_weight <= 0) {
    stop("Parameter 'min_weight' must be a positive scalar.")
  }

  base_adj <- if (isTRUE(use_mst)) {
    .dfs_chain_max_spanning_tree(G)
  } else {
    A <- (0.5 * (G + t(G)))
    diag(A) <- 0
    (A > 0) * 1
  }

  ord <- .dfs_chain_order_from_adj(base_adj, start = start)

  G_chain <- matrix(0, nrow = k, ncol = k)
  for (i in seq_len(k - 1L)) {
    u <- ord[i]
    v <- ord[i + 1L]
    w <- G[u, v]
    if (!is.finite(w) || w <= 0) {
      w <- min_weight
    }
    G_chain[u, v] <- w
    G_chain[v, u] <- w
  }

  list(order = ord, G_chain = G_chain)
}

#' DFS-chain approximation for L1 fused regression.
#'
#' @description
#' Approximates the fusion graph by a DFS-induced chain, then solves with
#' `fusedLassoProximalNewOperator()`. This is a one-shot approximation method
#' for general graphs.
#'
#' @param X matrix of covariates (n by p)
#' @param Y response vector (length n)
#' @param groups group indicators (length n)
#' @param lambda sparsity hyperparameter (scalar)
#' @param gamma fusion hyperparameter (scalar)
#' @param G pairwise fusion weight matrix (k by k)
#' @param mu smoothness parameter for proximal updates
#' @param tol optimization tolerance
#' @param num.it maximum iterations
#' @param lam.max optional precomputed maximal eigenvalue
#' @param c.flag kept for API parity; ignored in this implementation
#' @param intercept whether to include per-group intercept terms
#' @param penalty.factors per-feature sparsity weights
#' @param conserve.memory whether to avoid materializing XX/XY lists
#' @param scaling whether to scale losses by group size
#' @param edge.block number of edges processed per block
#' @param diagnostics whether to record per-iteration diagnostics
#' @param trace_every record diagnostics every `trace_every` iterations
#' @param chain.use.mst if TRUE, build DFS order on a maximum spanning tree
#' @param chain.start starting node for DFS order
#' @param chain.min.weight fallback positive weight for chain edges
#'
#' @return Coefficient matrix from operator solver on DFS-chain approximation.
#' @export
fusedLassoProximalDFSChainApprox <- function(
  X,
  Y,
  groups,
  lambda,
  gamma,
  G,
  mu = 1e-04,
  tol = 1e-06,
  num.it = 1000,
  lam.max = NULL,
  c.flag = FALSE,
  intercept = TRUE,
  penalty.factors = NULL,
  conserve.memory = NULL,
  scaling = TRUE,
  edge.block = 256L,
  diagnostics = FALSE,
  trace_every = 1L,
  chain.use.mst = TRUE,
  chain.start = 1L,
  chain.min.weight = 1e-8
) {
  if (c.flag) {
    warning("fusedLassoProximalDFSChainApprox does not use c.flag; running R implementation.")
  }

  chain <- .dfs_chain_graph(
    G = G,
    use_mst = chain.use.mst,
    start = chain.start,
    min_weight = chain.min.weight
  )

  fit <- fusedLassoProximalNewOperator(
    X = X,
    Y = Y,
    groups = groups,
    lambda = lambda,
    gamma = gamma,
    G = chain$G_chain,
    mu = mu,
    tol = tol,
    num.it = num.it,
    lam.max = lam.max,
    c.flag = FALSE,
    intercept = intercept,
    penalty.factors = penalty.factors,
    conserve.memory = conserve.memory,
    scaling = scaling,
    edge.block = edge.block,
    diagnostics = diagnostics,
    trace_every = trace_every
  )

  attr(fit, "dfs_chain_order") <- chain$order
  fit
}
