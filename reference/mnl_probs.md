# Compute multinomial logistic class probabilities given covariates

Evaluates P(X = t \| Zp) for each observation using a multinomial logit
with one or more covariates and class-specific intercepts and slopes.

## Usage

``` r
mnl_probs(Zp, params)
```

## Arguments

- Zp:

  Numeric vector of length n, or numeric matrix of dimension n x P,
  where P is the number of covariates. A vector is treated as a single
  covariate (P = 1).

- params:

  List with elements `$b0` (length-T intercepts, reference = 0) and `$b`
  (length-T slopes when P = 1, or P x T slope matrix when P \> 1,
  reference class = 1). See `bk2018_params$covariate_params`.

## Value

An n x T matrix of class probabilities (rows sum to 1).

## Examples

``` r
# Single covariate: class membership probabilities for Zp = 1..5
mnl_probs(1:5, bk2018_params$covariate_params)
#>           [,1]       [,2]       [,3]
#> [1,] 0.2037951 0.78188368 0.01432126
#> [2,] 0.3842557 0.54234331 0.07340104
#> [3,] 0.4905617 0.25471421 0.25472410
#> [4,] 0.3842487 0.07339685 0.54235449
#> [5,] 0.2037890 0.01432027 0.78189074

# Multiple covariates (n = 5, P = 2)
Zp_mat <- matrix(rnorm(10), nrow = 5, ncol = 2)
params2 <- list(b0 = c(0, 0.5, -0.5), b = matrix(rnorm(6), nrow = 2, ncol = 3))
mnl_probs(Zp_mat, params2)
#>            [,1]       [,2]       [,3]
#> [1,] 0.28898001 0.51907007 0.19194992
#> [2,] 0.42139659 0.06281686 0.51578655
#> [3,] 0.25569690 0.59988291 0.14442019
#> [4,] 0.34009877 0.06106763 0.59883360
#> [5,] 0.06023329 0.92395210 0.01581461
```
