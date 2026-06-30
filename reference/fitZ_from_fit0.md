# Estimate covariate effects with measurement parameters fixed (two-step EM)

Fixes `mPhi` at `fit0$mPhi` and estimates multinomial logit coefficients
`mGamma` (Q x (T-1)) via an EM algorithm with a BFGS M-step.

## Usage

``` r
fitZ_from_fit0(
  fit0,
  data,
  Y.names,
  Zp.names,
  tol = 1e-06,
  maxIter = 200L,
  incomplete = FALSE,
  include.intercept = TRUE,
  rebase = "C1",
  starting_val = NULL,
  verbose = FALSE
)
```

## Arguments

- fit0:

  Output of `lca_step1()$fit0`.

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- Zp.names:

  Character vector of covariate column names.

- tol:

  Convergence tolerance. Default `1e-6`.

- maxIter:

  Maximum EM iterations. Default `200`.

- incomplete:

  Logical. FIML for partially missing indicators. See the `Missing Data`
  section of `vignette("tseLCA", package = "tseLCA")`. Default `FALSE`.

- include.intercept:

  Logical. Prepend intercept to covariate design matrix. Default `TRUE`.

- rebase:

  Character or integer. Reference class for the multinomial logit
  parameterization (e.g. `"C1"`, `"C2"`, or an integer). Default `"C1"`.
  Must match the `rebase` used in
  [`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md)
  so class column ordering is consistent.

- starting_val:

  Optional Q x (T-1) starting value matrix for `mGamma`.

- verbose:

  Logical. Print convergence messages. Default `FALSE`.

## Value

A list with the following elements:

- `mGamma`:

  Q x (T-1) numeric matrix of multinomial logit coefficients, where Q is
  the number of columns in the covariate design matrix (including
  intercept if `include.intercept = TRUE`). Rows are named by covariate,
  columns by non-reference class (e.g. `"C2"`, `"C3"`).

- `mPhi`:

  Expanded item parameter matrix (items x classes), fixed at `fit0$mPhi`
  throughout estimation.

- `vOmega`:

  Length-T vector of marginal class proportions implied by the final
  `mGamma`, computed as column means of the fitted class probability
  matrix.

- `LLKSeries`:

  Single-column matrix of observed-data log-likelihoods, one row per EM
  iteration. Useful for diagnosing convergence.

- `converged`:

  Logical. `TRUE` if the EM loop exited before `maxIter` iterations or
  if the final log-likelihood change was below `tol`.

- `n_obs`:

  Integer. Number of observations used in estimation after listwise
  deletion on covariates.

## Examples

``` r
# \donttest{
d  <- generate_data(200, "high", "covariate", seed = 1)
s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)

# Estimate two-step gamma with mPhi fixed at Step-1 values
fZ <- fitZ_from_fit0(
  fit0     = s1$fit0,
  data     = d,
  Y.names  = paste0("Y", 1:6),
  Zp.names = "Zp",
  verbose  = TRUE
)
#> fitZ EM converged in 9 iterations.
fZ$mGamma   # Q x (T-1) coefficient matrix
#>                  C2         C3
#> Intercept  1.988800 -3.1317130
#> Zp        -1.017498  0.9190021
fZ$converged
#> [1] TRUE
# }
```
