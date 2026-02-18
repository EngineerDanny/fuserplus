# Dense-sort L1 fusion solver for complete (dense) group graphs.
# The fusion operator is evaluated row-wise using sorting, avoiding explicit
# edge loops over choose(k, 2) pairs.

.dense_sort_graph_check <- function(G, tol = 1e-12, require_uniform_weights = TRUE) {
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    return(list(ok = FALSE, reason = "G must be a square matrix."))
  }
  k <- nrow(G)
  if (k <= 1L) {
    return(list(ok = TRUE, reason = "Trivial graph.", k = k, dense = TRUE, uniform = TRUE, weight = 0))
  }
  if (max(abs(G - t(G))) > tol) {
    return(list(ok = FALSE, reason = "G must be symmetric for dense-sort solver."))
  }

  off <- G[row(G) != col(G)]
  if (any(off <= tol)) {
    return(list(ok = FALSE, reason = "Dense-sort solver requires strictly positive off-diagonal fusion weights."))
  }

  dense <- all(abs(off) > tol)
  if (!dense) {
    return(list(
      ok = FALSE,
      reason = "Dense-sort solver currently supports complete graphs only (all off-diagonal entries non-zero).",
      k = k,
      dense = FALSE
    ))
  }

  w_ref <- off[1L]
  uniform <- max(abs(off - w_ref)) <= tol
  if (require_uniform_weights && !uniform) {
    return(list(
      ok = FALSE,
      reason = "Dense-sort solver currently requires uniform off-diagonal fusion weights.",
      k = k,
      dense = TRUE,
      uniform = FALSE
    ))
  }

  list(
    ok = TRUE,
    reason = "Supported dense complete graph.",
    k = k,
    dense = TRUE,
    uniform = uniform,
    weight = as.numeric(w_ref)
  )
}

.dense_sort_pairwise_abs_sum <- function(x) {
  # Computes sum_{u<v} |x_u - x_v| in O(k log k) via sorting.
  xs <- sort(as.numeric(x), method = "quick")
  k <- length(xs)
  if (k <= 1L) {
    return(0)
  }
  idx <- seq_len(k)
  as.numeric(sum((2 * idx - k - 1) * xs))
}

.dense_sort_fusion_operator_term <- function(
  B.fusion,
  gamma.eff,
  mu,
  return_stats = FALSE,
  workspace = NULL
) {
  p.use <- nrow(B.fusion)
  k <- ncol(B.fusion)
  use_ws <- !is.null(workspace) &&
    is.environment(workspace) &&
    !is.null(workspace$delta) &&
    is.matrix(workspace$delta) &&
    identical(dim(workspace$delta), c(p.use, k))

  if (use_ws) {
    delta <- workspace$delta
    delta[] <- 0
    if (is.null(workspace$xs) || length(workspace$xs) != k) {
      workspace$xs <- numeric(k)
    }
    if (is.null(workspace$pref) || length(workspace$pref) != (k + 1L)) {
      workspace$pref <- numeric(k + 1L)
    }
    if (is.null(workspace$l.sat) || length(workspace$l.sat) != k) {
      workspace$l.sat <- integer(k)
    }
    if (is.null(workspace$u.uns) || length(workspace$u.uns) != k) {
      workspace$u.uns <- integer(k)
    }
    if (is.null(workspace$tmp.num) || length(workspace$tmp.num) != k) {
      workspace$tmp.num <- numeric(k)
    }
    if (is.null(workspace$s.sorted) || length(workspace$s.sorted) != k) {
      workspace$s.sorted <- numeric(k)
    }
  } else {
    delta <- matrix(0, nrow = p.use, ncol = k)
  }

  if (gamma.eff <= 0 || k <= 1L) {
    if (!return_stats) {
      return(delta)
    }
    return(list(delta = delta, linear = 0, prox_sq = 0, prox_inf = 0))
  }

  cval <- gamma.eff / mu
  tau <- mu / gamma.eff
  i.seq <- seq_len(k)
  k.vec <- rep.int(k, k)

  if (use_ws) {
    xs <- workspace$xs
    pref <- workspace$pref
    l.sat <- workspace$l.sat
    u.uns <- workspace$u.uns
    tmp.num <- workspace$tmp.num
    s.sorted <- workspace$s.sorted
  } else {
    xs <- numeric(k)
    pref <- numeric(k + 1L)
    l.sat <- integer(k)
    u.uns <- integer(k)
    tmp.num <- numeric(k)
    s.sorted <- numeric(k)
  }

  for (j in seq_len(p.use)) {
    x <- B.fusion[j, ]
    ord <- order(x, method = "radix")
    xs[] <- x[ord]
    pref[1L] <- 0
    pref[2:(k + 1L)] <- cumsum(xs)

    # left saturated count: # {j < i : xs[j] <= xs[i] - tau}
    l.sat[] <- pmin(
      findInterval(xs - tau, xs, rightmost.closed = TRUE),
      i.seq - 1L
    )

    # right unsaturated upper index: max j with xs[j] < xs[i] + tau
    u.uns[] <- pmax(
      i.seq,
      findInterval(xs + tau, xs, left.open = TRUE)
    )

    n.left.uns <- (i.seq - 1L) - l.sat
    sum.left.uns <- pref[i.seq] - pref[l.sat + 1L]
    left.uns <- cval * (n.left.uns * xs - sum.left.uns)
    left.uns[n.left.uns <= 0L] <- 0

    n.right.uns <- u.uns - i.seq
    sum.right.uns <- pref[u.uns + 1L] - pref[i.seq + 1L]
    right.uns <- cval * (sum.right.uns - n.right.uns * xs)
    right.uns[n.right.uns <= 0L] <- 0

    right.sat <- k.vec - u.uns
    tmp.num[] <- l.sat + left.uns - right.uns - right.sat
    s.sorted[] <- tmp.num

    row.delta <- gamma.eff * s.sorted
    delta[j, ord] <- row.delta
  }

  if (!return_stats) {
    return(delta)
  }

  list(
    delta = delta,
    linear = NA_real_,
    prox_sq = NA_real_,
    prox_inf = 1
  )
}

.objective_l1_fused_dense_complete <- function(B, prep, lambda, gamma.eff, intercept) {
  B.pen <- if (intercept) B[seq_len(prep$p), , drop = FALSE] else B

  sparsity.w <- matrix(prep$penalty.factors, nrow = prep$p, ncol = prep$k)
  sparsity.pen <- lambda * sum(abs(B.pen) * sparsity.w)

  fusion.pen <- 0
  if (gamma.eff > 0 && prep$k > 1L) {
    fusion.pen <- gamma.eff * sum(vapply(seq_len(nrow(B.pen)), function(j) {
      .dense_sort_pairwise_abs_sum(B.pen[j, ])
    }, numeric(1L)))
  }

  loss <- 0
  if (!is.null(prep$XX)) {
    for (g in seq_len(prep$k)) {
      b <- B[, g, drop = FALSE]
      qf <- drop(t(b) %*% prep$XX[[g]] %*% b)
      lin <- drop(t(b) %*% prep$XY[[g]])
      loss <- loss + (qf - 2 * lin + prep$y2[g]) / (2 * prep$samp.sizes[g])
    }
  } else {
    for (g in seq_len(prep$k)) {
      resid <- prep$Y.list[[g]] - prep$X.list[[g]] %*% B[, g, drop = FALSE]
      loss <- loss + sum(resid^2) / (2 * prep$samp.sizes[g])
    }
  }

  loss + sparsity.pen + fusion.pen
}

.genFusedLassoProximalDenseSortCore <- function(
  prep,
  lambda,
  gamma,
  edge.weight,
  mu,
  tol,
  num.it,
  lam.max,
  intercept,
  trace_state = NULL
) {
  if (is.null(lam.max)) {
    lam.max <- .compute_lam_max(prep$XX, prep$X.list)
  }

  max.node.weight <- (prep$k - 1L) * edge.weight
  L.U <- lam.max + (lambda^2 + 2 * gamma^2 * max.node.weight) / mu
  L.U.inv <- 1 / L.U
  gamma.eff <- gamma * edge.weight

  state <- .init_operator_state(prep$p, prep$k, intercept)
  B.old <- state$B.old
  W <- state$W
  weighted.delta.f <- state$weighted.delta.f
  penalty.mat <- NULL
  if (!is.null(prep$penalty.factors)) {
    penalty.mat <- matrix(prep$penalty.factors, nrow = prep$p, ncol = prep$k)
  }
  fusion.workspace <- new.env(parent = emptyenv())
  fusion.workspace$delta <- matrix(0, nrow(B.old), ncol(B.old))
  fusion.workspace$xs <- numeric(prep$k)
  fusion.workspace$pref <- numeric(prep$k + 1L)
  fusion.workspace$l.sat <- integer(prep$k)
  fusion.workspace$u.uns <- integer(prep$k)
  fusion.workspace$tmp.num <- numeric(prep$k)
  fusion.workspace$s.sorted <- numeric(prep$k)

  iter <- 0L
  converged <- FALSE

  for (i in seq_len(num.it)) {
    iter <- i
    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      trace_state$iter_global <- trace_state$iter_global + 1L
    }

    if (!is.null(penalty.mat)) {
      B.sparsity <- B.old[seq_len(prep$p), , drop = FALSE] * penalty.mat
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

    need.fusion.stats <- !is.null(trace_state) && isTRUE(trace_state$enabled)
    fusion.out <- .dense_sort_fusion_operator_term(
      B.fusion = B.fusion,
      gamma.eff = gamma.eff,
      mu = mu,
      return_stats = need.fusion.stats,
      workspace = fusion.workspace
    )
    if (need.fusion.stats) {
      fusion.delta <- fusion.out$delta
    } else {
      fusion.delta <- fusion.out
    }

    delta.reg <- lambda * A.sparsity + fusion.delta
    delta.f <- delta.lik + delta.reg
    B.new <- W - L.U.inv * delta.f

    weighted.delta.f <- weighted.delta.f + L.U.inv * 0.5 * i * delta.f
    Z <- -weighted.delta.f
    W <- (i * B.new + 2 * Z) / (i + 2)

    improvement <- sum(abs(B.old - B.new))

    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      if (i == 1L || (i %% trace_state$trace_every) == 0L || i == num.it) {
        objective <- .objective_l1_fused_dense_complete(
          B = B.new,
          prep = prep,
          lambda = lambda,
          gamma.eff = gamma.eff,
          intercept = intercept
        )
          dual.linear <- sparsity.linear + fusion.out$linear
          dual.quad <- 0.5 * mu * (sparsity.sq + fusion.out$prox_sq)
          .append_trace_row(trace_state, list(
          stage = "dense_sort",
          iter_local = i,
          iter_global = trace_state$iter_global,
          objective = objective,
          dual_linear = dual.linear,
          dual_quadratic = dual.quad,
          dual_surrogate_reg = dual.linear - dual.quad,
            dual_feas_inf = max(sparsity.inf, fusion.out$prox_inf, na.rm = TRUE),
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
    active_edges = choose(prep$k, 2L),
    total_edges = choose(prep$k, 2L)
  )
}

#' Dense-sort L1 fusion solver for complete group graphs.
#'
#' @description
#' Uses a complete-graph dense-sort fusion operator to avoid explicit edge loops
#' in the fusion update. Falls back to `fusedLassoProximalNewOperator()` only
#' when graph assumptions are not met and `fallback = TRUE`.
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
#' @param c.flag Kept for API parity; ignored in this implementation
#' @param intercept whether to include per-group intercept terms
#' @param penalty.factors per-feature sparsity weights
#' @param conserve.memory whether to avoid materializing XX/XY lists
#' @param scaling whether to scale losses by group size
#' @param edge.block Kept for API parity; not used by dense-sort core
#' @param diagnostics whether to record per-iteration objective/KKT diagnostics
#' @param trace_every record diagnostics every `trace_every` iterations
#' @param require_uniform_weights if TRUE, require equal off-diagonal weights in `G`
#' @param graph_tol numerical tolerance for graph checks
#' @param fallback if TRUE, use operator solver when graph assumptions are not met
#'
#' @return Coefficient matrix (p by k, or (p+1) by k if intercept=TRUE)
#' @export
fusedLassoProximalDenseSortScaffold <- function(
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
  require_uniform_weights = TRUE,
  graph_tol = 1e-12,
  fallback = TRUE
) {
  .validate_l1new_common(lambda, gamma, mu, edge.block)
  if (c.flag) {
    warning("fusedLassoProximalDenseSortScaffold does not use c.flag; running R implementation.")
  }
  if (is.null(conserve.memory)) {
    conserve.memory <- (ncol(X) >= 10000)
  }
  if (!is.logical(conserve.memory) || length(conserve.memory) != 1L) {
    stop("Parameter 'conserve.memory' must be TRUE/FALSE or NULL.")
  }

  graph.info <- .dense_sort_graph_check(
    G = G,
    tol = graph_tol,
    require_uniform_weights = require_uniform_weights
  )
  if (!isTRUE(graph.info$ok)) {
    if (!isTRUE(fallback)) {
      stop(graph.info$reason)
    }
    warning(paste("Dense-sort assumptions not met; falling back to fusedLassoProximalNewOperator().", graph.info$reason))
    return(fusedLassoProximalNewOperator(
      X = X,
      Y = Y,
      groups = groups,
      lambda = lambda,
      gamma = gamma,
      G = G,
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
    ))
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
  trace.state <- .init_trace_state(diagnostics = diagnostics, trace_every = trace_every)

  fit.raw <- .genFusedLassoProximalDenseSortCore(
    prep = prep,
    lambda = lambda,
    gamma = gamma,
    edge.weight = graph.info$weight,
    mu = mu,
    tol = tol,
    num.it = num.it,
    lam.max = lam.max,
    intercept = intercept,
    trace_state = trace.state
  )

  .finalize_l1new_fit(
    fit.raw = fit.raw,
    prep = prep,
    lambda = lambda,
    intercept = intercept,
    trace_state = trace.state
  )
}
