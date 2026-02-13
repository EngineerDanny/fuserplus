# fuserplus

Fused lasso for high-dimensional regression over groups.

`fuserplus` includes:
- L1 fusion regression (`fusedLassoProximal`)
- L2 fusion regression (`fusedL2DescentGLMNet`)
- development solver variants for L1 benchmarking: `operator`, `operator_ws`, `dense_sort`, `dfs_chain`, `chain_specialized`

## Installation

```r
install.packages("remotes")
remotes::install_github("EngineerDanny/fuserplus")
```

## Quick Start (Package API)

```r
library(fuserplus)
set.seed(123)

k <- 4
p <- 100
n.group <- 15
sigma <- 0.05
groups <- rep(1:k, each = n.group)

beta <- matrix(0, p, k)
beta[1:5, ] <- 1

X <- matrix(rnorm(length(groups) * p), nrow = length(groups), ncol = p)
y <- numeric(length(groups))
for (g in 1:k) {
  idx <- which(groups == g)
  y[idx] <- X[idx, ] %*% beta[, g] + rnorm(length(idx), 0, sigma)
}

G <- matrix(1, k, k)

beta.l1 <- fusedLassoProximal(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-3, G = G,
  tol = 1e-5, num.it = 2000, intercept = FALSE
)

beta.l2 <- fusedL2DescentGLMNet(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-3, G = G
)
```

## How To Test Solver Variants

The solver variants are currently benchmark-oriented development entry points.
Use the benchmark scripts below to compare speed/accuracy across methods.

### 1) Compare core L1 methods

```bash
Rscript scripts/benchmark_l1_dense_sort_vs_operator.R \
  --mode sweep \
  --methods old_l1,operator,dense_sort_scaffold,dfs_chain,chain_specialized \
  --k_values 40,80,120,160 \
  --reps 2 \
  --p 300 \
  --n_group_train 20 \
  --n_group_test 10 \
  --lambda 1e-3 \
  --gamma 2.0 \
  --num_it 700 \
  --graph_mode dense_uniform
```

Graph modes:
- `dense_uniform`
- `dense_nonuniform`
- `sparse`

### 2) Compare operator vs working-set (`operator_ws`)

```bash
Rscript scripts/benchmark_l1_operator_ws.R \
  --mode baseline \
  --k 120 \
  --p 100 \
  --n_group_train 10 \
  --n_group_test 5 \
  --lambda 1e-3 \
  --gamma 1e-3 \
  --num_it 1200 \
  --g_structure dense
```

## Run Package Tests

```bash
R -q -e "testthat::test_dir('tests/testthat')"
```

## Notes On Solver Assumptions

- `dense_sort` is intended for dense complete graphs with near-uniform off-diagonal weights.
- `dfs_chain` and `chain_specialized` are chain-based approximations for general graphs.
- `operator` / `operator_ws` are edge-explicit methods and are the general reference paths.

For detailed modeling examples, see `vignettes/subgroup_fusion.Rmd`.
