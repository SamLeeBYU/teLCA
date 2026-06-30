# Extract the variance-covariance matrix from a tseLCA model object

For measurement models, returns the BHHH variance-covariance matrix in
the unconstrained log-ratio parameterization (NOT the probability
scale). Row and column names identify each parameter as `log(pi_t/pi_1)`
(class prevalences) or `log(P(Y=k|C_t)/P(Y=0|C_t))` (item-response
probabilities). An attribute `"parameterization"` is attached to remind
the user of the scale.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
vcov(object, boundary.tol = 0.01, ...)

# S3 method for class 'tseLCA_covariate'
vcov(object, which = c("three_step", "two_step"), ...)

# S3 method for class 'tseLCA_distal'
vcov(object, ...)

# S3 method for class 'tseLCA_both'
vcov(
  object,
  which = c("covariate", "distal", "both"),
  step = c("three_step", "two_step"),
  ...
)
```

## Arguments

- object:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- boundary.tol:

  Scalar. Parameters within this tolerance of 0 or 1 are treated as
  fixed. Default `1e-2`.

- ...:

  Further arguments (currently unused).

- which:

  Character. `"three_step"` (default) or `"two_step"` for covariate
  models; `"covariate"`, `"distal"`, or `"both"` for both models.

- step:

  Character. For `tseLCA_both`: `"three_step"` (default) or
  `"two_step"`.

## Value

A named square matrix in the unconstrained log-ratio parameterization.
Row/column names identify each parameter as `log(pi_t/pi_1)` or
`log(P(Y=k|C_t)/P(Y=0|C_t))`. An attribute `"parameterization"` is
attached as a reminder. Returns `NULL` invisibly if `fit0$mU` is not
available. For structural models, returns the Step-3 vcov matrix; the
two-step vcov is only available when `get.twostep.vcov = TRUE`.

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
V <- vcov(fit_m)
# Names show log-ratio parameterization:
rownames(V)
#>  [1] "log(pi_C2/pi_C1)"           "log(pi_C3/pi_C1)"          
#>  [3] "log(P(Y1=1|C1)/P(Y1=0|C1))" "log(P(Y2=1|C1)/P(Y2=0|C1))"
#>  [5] "log(P(Y3=1|C1)/P(Y3=0|C1))" "log(P(Y4=1|C1)/P(Y4=0|C1))"
#>  [7] "log(P(Y5=1|C1)/P(Y5=0|C1))" "log(P(Y6=1|C1)/P(Y6=0|C1))"
#>  [9] "log(P(Y1=1|C2)/P(Y1=0|C2))" "log(P(Y2=1|C2)/P(Y2=0|C2))"
#> [11] "log(P(Y3=1|C2)/P(Y3=0|C2))" "log(P(Y4=1|C2)/P(Y4=0|C2))"
#> [13] "log(P(Y5=1|C2)/P(Y5=0|C2))" "log(P(Y6=1|C2)/P(Y6=0|C2))"
#> [15] "log(P(Y1=1|C3)/P(Y1=0|C3))" "log(P(Y2=1|C3)/P(Y2=0|C3))"
#> [17] "log(P(Y3=1|C3)/P(Y3=0|C3))" "log(P(Y4=1|C3)/P(Y4=0|C3))"
#> [19] "log(P(Y5=1|C3)/P(Y5=0|C3))" "log(P(Y6=1|C3)/P(Y6=0|C3))"
attr(V, "parameterization")
#> [1] "log-ratio (unconstrained); NOT probabilities"
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
vcov(fit)   # Q*(T-1) x Q*(T-1) vcov matrix with named rows/cols
#>              Intercept:C2        Zp:C2 Intercept:C3        Zp:C3
#> Intercept:C2  0.391644881 -0.173583653  0.001643327 -0.002599746
#> Zp:C2        -0.173583653  0.090099886  0.016315352 -0.002300251
#> Intercept:C3  0.001643327  0.016315352  0.517169347 -0.130664301
#> Zp:C3        -0.002599746 -0.002300251 -0.130664301  0.035941355
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
vcov(fit)   # T x T vcov matrix with mu_C1..mu_CT row/col names
#>               mu_C1         mu_C2         mu_C3
#> mu_C1  1.365920e-02  5.055316e-05 -8.311167e-05
#> mu_C2  5.055316e-05  1.301282e-02 -9.122807e-04
#> mu_C3 -8.311167e-05 -9.122807e-04  2.342768e-02
# }
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 0.5)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", Zo.name = "Zo",
                  use.simple.cov = TRUE)
vcov(fit, which = "covariate")
#>              Intercept:C2        Zp:C2 Intercept:C3        Zp:C3
#> Intercept:C2  0.391644881 -0.173583653  0.001643327 -0.002599746
#> Zp:C2        -0.173583653  0.090099886  0.016315352 -0.002300251
#> Intercept:C3  0.001643327  0.016315352  0.517169347 -0.130664301
#> Zp:C3        -0.002599746 -0.002300251 -0.130664301  0.035941355
vcov(fit, which = "distal")
#>               mu_C1         mu_C2         mu_C3
#> mu_C1  3.772208e-03  3.058231e-05 -1.253885e-05
#> mu_C2  3.058231e-05  5.147141e-03 -6.123296e-06
#> mu_C3 -1.253885e-05 -6.123296e-06  3.590917e-03
# }
```
