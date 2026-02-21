# Chain-specialized L1 fused solver.
# Uses a DFS-induced chain graph and a vectorized chain fusion operator.

.chain_fusion_operator_term <- function(B.fusion, edge.w, gamma, mu, return_stats = FALSE) {
  k <- ncol(B.fusion)
  if (gamma <= 0 || k <= 1L || length(edge.w) == 0L) {
    zero <- matrix(0, nrow(B.fusion), k)
    if (!return_stats) {
      return(zero)
    }
    return(list(delta = zero, linear = 0, prox_sq = 0, prox_inf = 0))
  }

  D <- B.fusion[, seq_len(k - 1L), drop = FALSE] - B.fusion[, 2:k, drop = FALSE]
  Z <- gamma * sweep(D, 2L, edge.w, `*`)
  P <- pmax(pmin(Z / mu, 1), -1)
  C <- gamma * sweep(P, 2L, edge.w, `*`)

  out <- matrix(0, nrow(B.fusion), k)
  out[, seq_len(k - 1L)] <- out[, seq_len(k - 1L), drop = FALSE] + C
  out[, 2:k] <- out[, 2:k, drop = FALSE] - C

  if (!return_stats) {
    return(out)
  }

  list(
    delta = out,
    linear = sum(Z * P),
    prox_sq = sum(P^2),
    prox_inf = max(abs(P))
  )
}

.chain_max_node_weight <- function(edge.w) {
  k <- length(edge.w) + 1L
  if (k <= 1L || length(edge.w) == 0L) {
    return(0)
  }

  deg <- numeric(k)
  deg[1L] <- edge.w[1L]
  deg[k] <- edge.w[k - 1L]
  if (k > 2L) {
    deg[2:(k - 1L)] <- edge.w[1:(k - 2L)] + edge.w[2:(k - 1L)]
  }
  max(deg)
}

.run_chain_iterations <- function(
    prep,
    lambda,
    gamma,
    mu,
    tol,
    num.it,
    lam.max,
    intercept,
    edge.w,
    trace_state = NULL) {
  if (is.null(lam.max)) {
    lam.max <- .compute_lam_max(prep$XX, prep$X.list)
  }

  max.node.weight <- .chain_max_node_weight(edge.w)
  L.U <- lam.max + (lambda^2 + 2 * gamma^2 * max.node.weight) / mu
  L.U.inv <- 1 / L.U

  state <- .init_operator_state(prep$p, prep$k, intercept)
  B.old <- state$B.old
  W <- state$W
  weighted.delta.f <- state$weighted.delta.f

  iter <- 0L
  converged <- FALSE

  for (i in seq_len(num.it)) {
    iter <- i
    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      trace_state$iter_global <- trace_state$iter_global + 1L
    }

    if (!is.null(prep$penalty.factors)) {
      B.sparsity <- B.old[seq_len(prep$p), , drop = FALSE] *
        matrix(prep$penalty.factors, nrow = prep$p, ncol = prep$k)
    } else {
      B.sparsity <- B.old[seq_len(prep$p), , drop = FALSE]
    }

    if (intercept) {
      B.sparsity <- rbind(B.sparsity, 0)
      B.fusion <- rbind(B.old[seq_len(prep$p), , drop = FALSE], 0)
    } else {
      B.fusion <- B.old
    }

    delta.lik <- .compute_delta_lik(
      XX = prep$XX,
      XY = prep$XY,
      X.list = prep$X.list,
      Y = prep$Y,
      groups = prep$groups,
      group.names = prep$group.names,
      B.old = B.old,
      samp.sizes = prep$samp.sizes
    )

    z.sparsity <- lambda * B.sparsity
    A.sparsity <- .clip_unit(z.sparsity / mu)
    sparsity.linear <- sum(z.sparsity * A.sparsity)
    sparsity.sq <- sum(A.sparsity^2)
    sparsity.inf <- max(abs(A.sparsity))

    fusion.stats <- .chain_fusion_operator_term(
      B.fusion = B.fusion,
      edge.w = edge.w,
      gamma = gamma,
      mu = mu,
      return_stats = TRUE
    )

    delta.reg <- lambda * A.sparsity + fusion.stats$delta
    delta.f <- delta.lik + delta.reg
    B.new <- W - L.U.inv * delta.f

    weighted.delta.f <- weighted.delta.f + L.U.inv * 0.5 * i * delta.f
    Z <- -weighted.delta.f
    W <- (i * B.new + 2 * Z) / (i + 2)

    improvement <- sum(abs(B.old - B.new))

    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      if (i == 1L || (i %% trace_state$trace_every) == 0L || i == num.it) {
        objective <- .objective_l1_fused(
          B = B.new,
          prep = prep,
          lambda = lambda,
          gamma = gamma,
          intercept = intercept
        )
        .append_trace_row(trace_state, list(
          stage = "chain_specialized",
          iter_local = i,
          iter_global = trace_state$iter_global,
          objective = objective,
          dual_linear = sparsity.linear + fusion.stats$linear,
          dual_quadratic = 0.5 * mu * (sparsity.sq + fusion.stats$prox_sq),
          dual_surrogate_reg = (sparsity.linear + fusion.stats$linear) -
            (0.5 * mu * (sparsity.sq + fusion.stats$prox_sq)),
          dual_feas_inf = max(sparsity.inf, fusion.stats$prox_inf),
          update_l1 = improvement,
          kkt_inf = max(abs(delta.f)),
          kkt_fro = sqrt(mean(delta.f^2))
        ))
      }
    }

    if (improvement < tol * prep$p) {
      converged <- TRUE
      break
    }

    B.old <- B.new
  }

  list(
    B = B.new,
    iterations = iter,
    converged = converged,
    active_edges = length(edge.w),
    total_edges = length(edge.w)
  )
}

.reorder_prep_for_chain <- function(prep, order_idx) {
  out <- prep
  out$group.names <- prep$group.names[order_idx]
  out$samp.sizes <- prep$samp.sizes[order_idx]

  if (!is.null(prep$XX)) {
    out$XX <- prep$XX[order_idx]
  }
  if (!is.null(prep$XY)) {
    out$XY <- prep$XY[order_idx]
  }
  if (!is.null(prep$X.list)) {
    out$X.list <- prep$X.list[order_idx]
  }
  if (!is.null(prep$Y.list)) {
    out$Y.list <- prep$Y.list[order_idx]
  }
  if (!is.null(prep$y2)) {
    out$y2 <- prep$y2[order_idx]
  }

  k <- prep$k
  edge.u <- seq_len(k - 1L)
  edge.v <- edge.u + 1L
  out$edge.u <- edge.u
  out$edge.v <- edge.v
  out
}

#' Chain-specialized DFS approximation for fused L1 regression.
#'
#' @description
#' Builds a DFS-induced chain from `G`, reorders groups along the chain, and
#' runs a chain-specialized vectorized solver for the fusion operator.
#'
#' @inheritParams fusedLassoProximalDFSChainApprox
#'
#' @return Coefficient matrix (p by k, or (p+1) by k if intercept=TRUE).
#' @export
fusedLassoProximalChainSpecialized <- function(
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
    chain.min.weight = 1e-8) {
  .validate_l1new_common(lambda, gamma, mu, edge.block)
  if (c.flag) {
    warning("fusedLassoProximalChainSpecialized does not use c.flag; running R implementation.")
  }
  if (is.null(conserve.memory)) {
    conserve.memory <- (ncol(X) >= 10000)
  }
  if (!is.logical(conserve.memory) || length(conserve.memory) != 1L) {
    stop("Parameter 'conserve.memory' must be TRUE/FALSE or NULL.")
  }

  chain <- .dfs_chain_graph(
    G = G,
    use_mst = chain.use.mst,
    start = chain.start,
    min_weight = chain.min.weight
  )
  ord <- chain$order
  k <- nrow(G)
  if (k <= 1L) {
    ord <- seq_len(k)
  }

  # Edge weights on adjacent chain positions in chain order.
  edge.w <- numeric(max(0L, k - 1L))
  if (k > 1L) {
    for (i in seq_len(k - 1L)) {
      w <- G[ord[i], ord[i + 1L]]
      if (!is.finite(w) || w <= 0) {
        w <- chain.min.weight
      }
      edge.w[i] <- w
    }
  }

  prep <- .prepare_l1new_inputs(
    X = X,
    Y = Y,
    groups = groups,
    G = G,
    intercept = intercept,
    penalty.factors = penalty.factors,
    conserve.memory = conserve.memory,
    scaling = scaling
  )
  prep.chain <- .reorder_prep_for_chain(prep, ord)
  prep.chain$edge.w <- edge.w

  trace.state <- .init_trace_state(diagnostics = diagnostics, trace_every = trace_every)

  fit.raw <- .run_chain_iterations(
    prep = prep.chain,
    lambda = lambda,
    gamma = gamma,
    mu = mu,
    tol = tol,
    num.it = num.it,
    lam.max = lam.max,
    intercept = intercept,
    edge.w = edge.w,
    trace_state = trace.state
  )

  B.chain <- .finalize_l1new_fit(
    fit.raw = fit.raw,
    prep = prep.chain,
    lambda = lambda,
    intercept = intercept,
    trace_state = trace.state
  )

  # Reorder columns back to original group order for compatibility.
  inv.ord <- order(ord)
  B.out <- B.chain[, inv.ord, drop = FALSE]
  colnames(B.out) <- prep$group.names
  B.out
}
