# Backward-compatible wrapper over the split L1-new solvers.
# Supports operator-only, operator+working-set, and chain/dense special cases
# through a single `solver` entry point.

#' Fused L1 solver (compatibility wrapper).
#'
#' @param ... Arguments forwarded to the selected solver.
#' @param solver Solver variant. One of `"operator"`, `"operator_ws"`,
#'   `"dfs_chain"`, `"chain_specialized"`, `"dense_sort"`.
#'   If NULL, dispatch is controlled by `working_set` for backward compatibility.
#' @param working_set Deprecated dispatch switch for backward compatibility.
#'   If `solver` is NULL, `TRUE -> "operator_ws"` and `FALSE -> "operator"`.
#' @param ws_init_edges Initial number of active edges for working-set mode.
#' @param ws_add_edges Max number of edges added per outer working-set round.
#' @param ws_max_outer Max number of outer working-set rounds.
#' @param ws_inner_it Inner-iteration budget per outer round.
#' @param ws_violation_tol Violation threshold for edge activation.
#' @param ws_final_full If TRUE, run a final full-edge refinement pass.
#' @param ws_final_it Iteration budget for final full-edge refinement.
#' @param screening Edge screening strategy in working-set mode: "none" or "grad_zero".
#' @param screen_margin Screening margin in [0,1] for `screening = "grad_zero"`.
#' @param screen_max_drop_frac Maximum fraction of edges allowed to be screened.
#' @param screen_min_keep Minimum number of edges to keep after screening.
#' @param chain.use.mst If TRUE, build chain from maximum spanning tree (chain modes).
#' @param chain.start DFS start node (chain modes).
#' @param chain.min.weight Positive fallback edge weight (chain modes).
#' @param require_uniform_weights If TRUE, dense-sort requires uniform off-diagonal
#'   weights in `G`.
#' @param graph_tol Numerical tolerance for dense-sort graph checks.
#' @param fallback If TRUE, dense-sort falls back to operator when assumptions fail.
#'
#' @return Coefficient matrix.
#' @export
fusedLassoProximalNew <- function(
  ...,
  solver = NULL,
  working_set = FALSE,
  ws_init_edges = 256L,
  ws_add_edges = 256L,
  ws_max_outer = 8L,
  ws_inner_it = NULL,
  ws_violation_tol = 1e-06,
  ws_final_full = TRUE,
  ws_final_it = NULL,
  screening = c("none", "grad_zero"),
  screen_margin = 0,
  screen_max_drop_frac = 1,
  screen_min_keep = 0L,
  chain.use.mst = TRUE,
  chain.start = 1L,
  chain.min.weight = 1e-8,
  require_uniform_weights = TRUE,
  graph_tol = 1e-12,
  fallback = TRUE
) {
  if (!is.logical(working_set) || length(working_set) != 1L) {
    stop("Parameter 'working_set' must be TRUE/FALSE.")
  }
  screening <- match.arg(screening)
  if (is.null(solver)) {
    solver <- if (working_set) "operator_ws" else "operator"
  }
  solver <- match.arg(
    solver,
    choices = c("operator", "operator_ws", "dfs_chain", "chain_specialized", "dense_sort")
  )

  if (solver == "operator") {
    fusedLassoProximalNewOperator(...)
  } else if (solver == "operator_ws") {
    fusedLassoProximalNewWorkingSet(
      ...,
      ws_init_edges = ws_init_edges,
      ws_add_edges = ws_add_edges,
      ws_max_outer = ws_max_outer,
      ws_inner_it = ws_inner_it,
      ws_violation_tol = ws_violation_tol,
      ws_final_full = ws_final_full,
      ws_final_it = ws_final_it,
      screening = screening,
      screen_margin = screen_margin,
      screen_max_drop_frac = screen_max_drop_frac,
      screen_min_keep = screen_min_keep
    )
  } else if (solver == "dfs_chain") {
    fusedLassoProximalDFSChainApprox(
      ...,
      chain.use.mst = chain.use.mst,
      chain.start = chain.start,
      chain.min.weight = chain.min.weight
    )
  } else if (solver == "chain_specialized") {
    fusedLassoProximalChainSpecialized(
      ...,
      chain.use.mst = chain.use.mst,
      chain.start = chain.start,
      chain.min.weight = chain.min.weight
    )
  } else if (solver == "dense_sort") {
    fusedLassoProximalDenseSortScaffold(
      ...,
      require_uniform_weights = require_uniform_weights,
      graph_tol = graph_tol,
      fallback = fallback
    )
  } else {
    stop("Unknown solver: ", solver)
  }
}
