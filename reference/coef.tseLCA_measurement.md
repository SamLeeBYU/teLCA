# Extract coefficients from a tseLCA model object

Extract coefficients from a tseLCA model object

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
coef(object, ...)

# S3 method for class 'tseLCA_covariate'
coef(object, which = c("three_step", "two_step"), ...)

# S3 method for class 'tseLCA_distal'
coef(object, ...)

# S3 method for class 'tseLCA_both'
coef(
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

- ...:

  Further arguments (currently unused).

- which:

  Character. For covariate and both models: `"three_step"` (default) or
  `"two_step"`. For both models also accepts `"covariate"`, `"distal"`,
  or `"both"`.

- step:

  Character. For `tseLCA_both`: `"three_step"` (default) or
  `"two_step"`.

## Value

The coefficient matrix (covariate models), named numeric vector (distal
models), or a named list of both (measurement or both models).

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
coef(fit_m)   # returns list with $prevalences and $item_probs
#> $prevalences
#>                
#> P(C1) 0.2755108
#> P(C2) 0.3497761
#> P(C3) 0.3747131
#> 
#> $item_probs
#>                C1         C2         C3
#> P(Y1|C) 0.8596137 0.93397796 0.04399101
#> P(Y2|C) 0.8167472 0.89269738 0.16741645
#> P(Y3|C) 0.9999900 0.77243231 0.06487957
#> P(Y4|C) 0.8302998 0.05619882 0.13766913
#> P(Y5|C) 0.8222312 0.03677733 0.05498234
#> P(Y6|C) 0.7656717 0.20241699 0.07538715
#> 
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
coef(fit)                      # three-step estimates
#>                  C2         C3
#> Intercept  2.233384 -3.2742157
#> Zp        -1.156988  0.9400712
coef(fit, which = "two_step")  # two-step starting values
#>                  C2         C3
#> Intercept  1.988800 -3.1317130
#> Zp        -1.017498  0.9190021
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
coef(fit)   # named vector of class means
#>       mu_C1       mu_C2       mu_C3 
#> -0.82227322  1.09461008  0.04916818 
# }
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", Zo.name = "Zo",
                  use.simple.cov = TRUE)
coef(fit, which = "covariate")
#>                  C2         C3
#> Intercept  2.233384 -3.2742157
#> Zp        -1.156988  0.9400712
coef(fit, which = "distal")
#>       mu_C1       mu_C2       mu_C3 
#> -1.04363124  0.00515846  1.14210878 
coef(fit, which = "both")
#> $covariate
#>                  C2         C3
#> Intercept  2.233384 -3.2742157
#> Zp        -1.156988  0.9400712
#> 
#> $distal
#>       mu_C1       mu_C2       mu_C3 
#> -1.04363124  0.00515846  1.14210878 
#> 
# }
```
