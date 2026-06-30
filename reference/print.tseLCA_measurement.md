# Print a tseLCA model object

Compact one-line or table summary printed to the console.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
print(x, ...)

# S3 method for class 'tseLCA_covariate'
print(x, digits = 4, ...)

# S3 method for class 'tseLCA_distal'
print(x, digits = 4, ...)

# S3 method for class 'tseLCA_both'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- ...:

  Further arguments.

- digits:

  Integer. Number of decimal places for coefficient tables.

## Value

Invisibly returns `x`.

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
print(fit_m)
#> tseLCA -- measurement model
#>   Classes: 3   Log-lik: -295.0095   AIC: 630.02   BIC: 682.12
#>   Entropy R²: 0.8683
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
print(fit)
#> tseLCA -- three-step covariate model
#>   Classes: 3   Estimator: ML   Log-lik: -548.6403   AIC: 1177.28   BIC: 1309.21
#>   Entropy R² (covariate-adjusted): 0.8589
#> 
#> Covariate coefficients (three-step):
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.2334    0.6258  3.5688 < 0.001 ***
#> Zp:C2         -1.1570    0.3002 -3.8545 < 0.001 ***
#> Intercept:C3  -3.2742    0.7191 -4.5529 < 0.001 ***
#> Zp:C3          0.9401    0.1896  4.9587 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
# \donttest{
d   <- generate_data(200, "high", "distal", seed = 2)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zo.name = "Zo", use.simple.cov = TRUE)
print(fit)
#> tseLCA -- three-step distal outcome model
#>   Classes: 3   Estimator: ML   Family: gaussian
#>   Log-lik: -892.7558   AIC: 1831.51   BIC: 1907.37
#> 
#> Distal outcome means by class:
#>              Estimate Std.Error z.value     p.value
#> mu_C1 (mean)  -0.8223    0.1169 -7.0356 < 0.001 ***
#> mu_C2 (mean)   1.0946    0.1141  9.5956 < 0.001 ***
#> mu_C3 (mean)   0.0492    0.1531  0.3212 0.7480     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
