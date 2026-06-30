# Summarize a tseLCA model object

Verbose summary including model fit, class prevalences, item-response
probabilities, and coefficient tables with standard errors and p-values.

## Usage

``` r
# S3 method for class 'tseLCA_measurement'
summary(object, ...)

# S3 method for class 'tseLCA_covariate'
summary(object, digits = 4, ...)

# S3 method for class 'tseLCA_distal'
summary(object, digits = 4, ...)

# S3 method for class 'tseLCA_both'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  A `tseLCA` object returned by
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- ...:

  Further arguments (currently unused).

- digits:

  Integer. Number of decimal places for coefficient tables.

## Value

Invisibly returns `object`.

## Examples

``` r
d    <- generate_data(100, "high", "covariate", seed = 1)
fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
summary(fit_m)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -295.0095
#> AIC            : 630.0190
#> BIC            : 682.1224
#> Entropy R²     : 0.8683
#> 
#> Class prevalences:
#>             
#> P(C1) 0.2755
#> P(C2) 0.3498
#> P(C3) 0.3747
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8596 0.9340 0.0440
#> P(Y2|C) 0.8167 0.8927 0.1674
#> P(Y3|C) 1.0000 0.7724 0.0649
#> P(Y4|C) 0.8303 0.0562 0.1377
#> P(Y5|C) 0.8222 0.0368 0.0550
#> P(Y6|C) 0.7657 0.2024 0.0754
# \donttest{
d   <- generate_data(200, "high", "covariate", seed = 1)
fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
                  Zp.names = "Zp", use.simple.cov = TRUE)
summary(fit)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -548.6403
#> AIC            : 1177.2805
#> BIC            : 1309.2132
#> Entropy R²     : 0.8589  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.9888 -3.1317
#> Zp        -1.0175  0.9190
#> 
#> Three-step estimates:
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
summary(fit)
#> -- tseLCA Three-Step Distal Outcome Model -------------------
#> Latent classes : 3
#> Estimator      : ML
#> Family         : gaussian
#> Log-likelihood : -892.7558
#> AIC            : 1831.5116
#> BIC            : 1907.3729
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error z.value     p.value
#> mu_C1 (mean)  -0.8223    0.1169 -7.0356 < 0.001 ***
#> mu_C2 (mean)   1.0946    0.1141  9.5956 < 0.001 ***
#> mu_C3 (mean)   0.0492    0.1531  0.3212 0.7480     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
