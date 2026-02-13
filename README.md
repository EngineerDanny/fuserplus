# fuserplus
Fused lasso for high-dimensional regression over groups. This package implements the model described in [Dondelinger et al. (2016)](https://arxiv.org/abs/1611.00953).

`fuserplus` provides L1/L2 fusion estimators and multiple L1 solver backends:
`operator`, `operator_ws`, `dfs_chain`, `chain_specialized`, and `dense_sort`.

## Installation

```r
library('devtools')
install_github('EngineerDanny/fuserplus')
```
## Example

See also the included vignette.

```r
library(fuserplus)
set.seed(123)

# Generate simple heterogeneous dataset
k = 4 # number of groups
p = 100 # number of covariates
n.group = 15 # number of samples per group
sigma = 0.05 # observation noise sd
groups = rep(1:k, each=n.group) # group indicators

# sparse linear coefficients
beta = matrix(0, p, k)
nonzero.ind = rbinom(p*k, 1, 0.025/k) # Independent coefficients
nonzero.shared = rbinom(p, 1, 0.025) # shared coefficients
beta[which(nonzero.ind==1)] = rnorm(sum(nonzero.ind), 1, 0.25) 
beta[which(nonzero.shared==1),] = rnorm(sum(nonzero.shared), -1, 0.25)

X = lapply(1:k, function(k.i) matrix(rnorm(n.group*p),n.group, p)) # covariates 
y = sapply(1:k, function(k.i) X[[k.i]] %*% beta[,k.i] + rnorm(n.group, 0, sigma)) # response
X = do.call('rbind', X)

# Pairwise Fusion strength hyperparameters (tau(k,k'))
# Same for all pairs in this example
G = matrix(1, k, k) 

# Use L1 fusion to estimate betas (with near-optimal sparsity and 
# information sharing among groups)
beta.l1 = fusedLassoProximal(X, y, groups, lambda=0.001, tol=9e-5, 
                             gamma=0.001, G, intercept=FALSE,
                             num.it=2000) 

# Use L2 fusion to estimate betas (with near-optimal information sharing among groups)
# Note: fusedL2DescentGLMNet now accepts the original X/y/groups and performs
# the block-diagonal transformation internally.
beta.l2 = fusedL2DescentGLMNet(X, y, groups,
                               lambda=0.001,
                               G=G,
                               gamma=0.001)
```
