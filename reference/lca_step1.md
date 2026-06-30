# Fit the LCA measurement model (Step 1)

Estimates the latent class measurement model with multilevLCA and
optionally, fixes `mPhi` and estimates covariate effects (two-step
initialization) with
[`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).

## Usage

``` r
lca_step1(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-08,
  covariate.tol = 1e-06,
  iter.measurement = 10L,
  R2.threshold = 0.7,
  use.two.step = TRUE,
  estimate.one.step = TRUE,
  incomplete = FALSE,
  maxIter.fitZ = 200L,
  include.intercept = TRUE,
  rebase = "C1",
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame containing at minimum the indicator columns.

- Y.names:

  Character vector of item column names.

- n_classes:

  Integer. Number of latent classes.

- Zp.names:

  Character vector of covariate column names, or `NULL`.

- maxIter.measurement:

  Maximum EM iterations before giving up on convergence. Default
  `5000L`.

- measurement.tol:

  Convergence tolerance. Default `1e-8`.

- covariate.tol:

  Convergence tolerance for the `fitZ` M-step. Default `1e-6`.

- iter.measurement:

  Number of random restarts when entropy R\\^2\\ is low. Default `10`.

- R2.threshold:

  Entropy R\\^2\\ below which restarts are triggered. Default `0.7`.

- use.two.step:

  Logical. If `TRUE`, also estimate `fitZ` with
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md)
  if `Zp.names` is applied. Default `TRUE`.

- estimate.one.step:

  Logical. If `FALSE`, skip the unconditional EM and only compute
  `fitZ`. Default `TRUE`.

- incomplete:

  Logical. FIML for partially missing indicators. See the `Missing Data`
  section of `vignette("tseLCA", package = "tseLCA")`. Default `FALSE`.

- maxIter.fitZ:

  Maximum EM iterations for
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).
  Default `200`.

- include.intercept:

  Logical. Prepend intercept to covariate design matrix. Default `TRUE`.

- rebase:

  Character or integer specifying the reference latent class. Use
  `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
  measurement model is permuted so this class becomes column 1, making
  it the reference for all downstream multinomial logit
  parameterizations.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A list with `$fit0`
([`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
measurement model) and `$fitZ` (two-step covariate model from
[`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md),
or `NULL`).

## Examples

``` r
# \donttest{
d <- generate_data(200, "high", "covariate", seed = 1)

# Measurement model only
s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
s1$fit0$vPi    # estimated class prevalences
#>                
#> P(C1) 0.3495138
#> P(C2) 0.2915216
#> P(C3) 0.3589645
s1$fit0$mPhi   # item-response probabilities
#>                C1         C2         C3
#> P(Y1|C) 0.8702096 0.79456644 0.12317767
#> P(Y2|C) 0.9016604 0.88525528 0.10247858
#> P(Y3|C) 0.8743309 0.87570434 0.06720021
#> P(Y4|C) 0.8565891 0.09127798 0.06686104
#> P(Y5|C) 0.8909744 0.09780804 0.02807791
#> P(Y6|C) 0.8206322 0.13853263 0.09135284

# With two-step covariate initialization
s1z <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3,
                 Zp.names = "Zp", use.two.step = TRUE, verbose = TRUE)
#> fitZ EM converged in 9 iterations.
s1z$fitZ$mGamma   # two-step gamma estimates
#>                  C2         C3
#> Intercept  1.988800 -3.1317130
#> Zp        -1.017498  0.9190021
# }
```
