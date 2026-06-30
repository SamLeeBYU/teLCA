# Estimate two-step covariate model via multilevLCA (optional reference path)

Calls
[`multilevLCA::multiLCA`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
with `fixedpars = 1` and `Z = Zp.names` to fit the two-step covariate
model. This is the original multilevLCA approach and is used when
`get.twostep.vcov = TRUE` in
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
to obtain multilevLCA's corrected standard errors for the two-step gamma
estimates.

## Usage

``` r
fitZ_from_multiLCA(
  data,
  Y.names,
  n_classes,
  Zp.names,
  maxIter.measurement,
  measurement.tol,
  covariate.tol,
  iter.measurement,
  R2.threshold,
  incomplete = FALSE,
  rebase = "C1",
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- n_classes:

  Integer. Number of latent classes.

- Zp.names:

  Character vector of covariate column names.

- maxIter.measurement:

  Maximum EM iterations.

- measurement.tol:

  Convergence tolerance.

- covariate.tol:

  NR tolerance for the covariate model.

- iter.measurement:

  Number of random restarts.

- R2.threshold:

  Entropy R\\^2\\ restart threshold.

- incomplete:

  Logical. FIML for partially missing indicators. See the `Missing Data`
  section of `vignette("tseLCA", package = "tseLCA")`. Default `FALSE`.

- rebase:

  Character or integer. Reference class for column naming of `$mGamma`.
  Must match the `rebase` used in
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
  so coefficient labels are consistent. Default `"C1"`.

- verbose:

  Logical.

## Value

A list with the following elements:

- `mGamma`:

  Q x (T-1) numeric matrix of multinomial logit coefficients. Rows are
  named by covariate (including `"Intercept"`), columns by non-reference
  class (e.g. `"C2"`, `"C3"`).

- `mPhi`:

  Item parameter matrix (items x classes) from the fixed-parameter
  multilevLCA fit.

- `vOmega`:

  Length-T vector of marginal class proportions, computed as the average
  of the fitted class probability matrix (`vPi_avg` in multilevLCA
  output).

- `LLKSeries`:

  Matrix of observed-data log-likelihoods across EM iterations, passed
  through directly from the multilevLCA fit.

- `raw_fit`:

  The full
  [`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
  output object, including `$Varmat_cor` (corrected variance matrix) and
  `$SEs_cor_gamma` (corrected standard errors for `mGamma`) if
  available.

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Two-step estimation via multiLCA (fixedpars = 1)
fZ_ml <- fitZ_from_multiLCA(
  data                = d,
  Y.names             = paste0("Y", 1:6),
  n_classes           = 3,
  Zp.names            = "Zp",
  maxIter.measurement = 5000L,
  measurement.tol     = 1e-8,
  covariate.tol       = 1e-6,
  iter.measurement    = 10L,
  R2.threshold        = 0.70
)
fZ_ml$mGamma           # two-step estimates
#>                  C2         C3
#> Intercept  1.990672 -3.1319910
#> Zp        -1.018352  0.9190157
fZ_ml$raw_fit$Varmat_cor   # multilevLCA corrected vcov
#>             [,1]         [,2]         [,3]         [,4]
#> [1,]  0.28553970 -0.119040291  0.051831951 -0.016507545
#> [2,] -0.11904029  0.065539163 -0.009799449  0.005370232
#> [3,]  0.05183195 -0.009799449  0.613479393 -0.157124127
#> [4,] -0.01650754  0.005370232 -0.157124127  0.043473866
# }
```
