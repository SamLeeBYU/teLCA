# tseLCA Workflow

``` r

library(tseLCA)
```

## Overview

`tseLCA` implements the BCH and ML bias-adjusted three-step estimators
for latent class analysis (LCA) with covariates and distal outcomes,
following the methodological framework for both BCH and Vermunt’s ML
approaches from Bakk, Tekle & Vermunt (2013). `tseLCA` also builds on
top of the two-step LCA estimation procedure outlined by Bakk & Kuha
(2018), and using the R package `multilevLCA` for efficient measurement
model estimation from Lyrvall et al. (2025). `tseLCA` provides analytic
sandwich variance estimation that propagates measurement uncertainty
through the classification-error correction in the final step.

The three-step approach separates the model into:

1.  **Step 1** — Estimate the LCA measurement model (class-conditional
    item probabilities and class prevalences).
2.  **Step 2** — Assign posterior class probabilities and compute the
    misclassification matrix.
3.  **Step 3** — Estimate the structural model (covariate effects or
    distal outcome means) using the bias-adjusted weights.

------------------------------------------------------------------------

## Synthetic data

The built-in data-generating process replicates the design of Bakk &
Kuha (2018). Each dataset has six binary indicators
($`Y_1, \ldots, Y_6`$) drawn from a three-class LCA, plus either a
covariate $`Z_p \sim \text{Uniform}\{1,\ldots,5\}`$ predicting class
membership, or a continuous distal outcome $`Z_o`$ predicted by class
membership.

``` r

# High separation: P(Y_h = 1 | class) = 0.9 / 0.1
d <- generate_data(
  n = 500,
  separation = "high",
  scenario = "covariate",
  seed = 1
)
head(d)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp
#> 1  1  1  1  0  0  0 2  1
#> 2  0  0  0  0  0  0 3  4
#> 3  1  0  1  0  0  0 2  1
#> 4  1  1  0  1  1  1 1  2
#> 5  0  0  0  0  0  0 3  5
#> 6  1  1  1  1  1  1 1  3
```

``` r

# Low separation: P(Y_h = 1 | class) = 0.7 / 0.3
# Zp and X are identical to 'd' because seed = 1
d.low <- generate_data(
  n = 500,
  separation = "low",
  scenario = "covariate",
  seed = 1
)
head(d.low)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp
#> 1  1  1  1  0  0  0 2  1
#> 2  1  0  0  1  0  1 3  4
#> 3  0  0  1  0  0  0 2  1
#> 4  1  1  0  0  1  1 1  2
#> 5  0  0  1  1  1  0 3  5
#> 6  1  1  1  1  1  1 1  3
```

------------------------------------------------------------------------

## Step 1: Measurement model

[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
with no `Zp.names` or `Zo.name` fits the measurement model only,
returning a `tseLCA_measurement` object. Internally this calls
[`multilevLCA::multiLCA()`](https://rdrr.io/pkg/multilevLCA/man/multiLCA.html)
with random restarts when entropy $`R^2`$ is low.

``` r

d.measurement <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  measurement.tol = 1e-8
)
summary(d.measurement)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -1455.5052
#> AIC            : 2951.0104
#> BIC            : 3035.3025
#> Entropy R²     : 0.8780
#> 
#> Class prevalences:
#>             
#> P(C1) 0.3570
#> P(C2) 0.3308
#> P(C3) 0.3122
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.9237 0.8621 0.1187
#> P(Y2|C) 0.9083 0.9219 0.1178
#> P(Y3|C) 0.9148 0.9571 0.0731
#> P(Y4|C) 0.8843 0.1481 0.0875
#> P(Y5|C) 0.8817 0.1340 0.1118
#> P(Y6|C) 0.9174 0.0889 0.1252
```

With low separation the measurement model can struggle to find the
global maximum. Use `iter.measurement` to trigger the number of random
restarts whenever entropy $`R^2`$ falls below `R2.threshold`.

``` r

d.low.measurement <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  iter.measurement = 10,
  R2.threshold = 0.9
)
summary(d.low.measurement)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -2019.2458
#> AIC            : 4078.4916
#> BIC            : 4162.7837
#> Entropy R²     : 0.3327
#> 
#> Class prevalences:
#>             
#> P(C1) 0.2753
#> P(C2) 0.4551
#> P(C3) 0.2696
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.7328 0.6418 0.3345
#> P(Y2|C) 0.5723 0.7549 0.3223
#> P(Y3|C) 0.6937 0.7101 0.3036
#> P(Y4|C) 0.6846 0.4588 0.2105
#> P(Y5|C) 0.6947 0.4058 0.3787
#> P(Y6|C) 0.8651 0.3551 0.2456
```

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) S3 method
delegates to `multilevLCA`’s item-profile plot.

``` r

plot(d.measurement)
```

![](tseLCA-workflow_files/figure-html/plot-measurement-1.png)

------------------------------------------------------------------------

## Two-step estimates

[`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md)
fixes the measurement parameters at their Step-1 values and estimates
multinomial logit coefficients $`\gamma`$ via EM. These two-step
estimates serve as starting values for Step 3 and are generally close to
the final three-step estimates.

``` r

d.fitZ <- fitZ_from_fit0(
  fit0 = d.measurement$measurement_model$fit0,
  data = d,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp"
)
# True slopes: -1 (C2) and +1 (C3) relative to C1
d.fitZ$mGamma
#>                   C2         C3
#> Intercept  2.1934130 -3.4524271
#> Zp        -0.9411383  0.8971774
```

Starting values from the high-separation fit can be passed to the
low-separation fit to help it converge.

``` r

d.low.fitZ <- fitZ_from_fit0(
  fit0 = d.low.measurement$measurement_model$fit0,
  data = d.low,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp",
  starting_val = d.fitZ$mGamma
)
d.low.fitZ$mGamma
#>                   C2         C3
#> Intercept  3.0446368 -3.6948005
#> Zp        -0.9832597  0.9487391
```

------------------------------------------------------------------------

## Three-step estimation

### ML estimator (default)

A single
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
call handles all three steps. By default it uses the ML correction of
Vermunt (2010) and modal class assignment.

``` r

d.three_step <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp"
)
summary(d.three_step)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0650
#> AIC            : 2758.1299
#> BIC            : 2926.7143
#> Entropy R²     : 0.8693  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.0411    0.3237  6.3050 < 0.001 ***
#> Zp:C2         -0.8821    0.1406 -6.2730 < 0.001 ***
#> Intercept:C3  -3.4836    0.5913 -5.8913 < 0.001 ***
#> Zp:C3          0.8985    0.1435  6.2606 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The standard [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) S3 methods work on any
`tseLCA` object.

``` r

coef(d.three_step)
#>                   C2         C3
#> Intercept  2.0410764 -3.4835616
#> Zp        -0.8820801  0.8984978
vcov(d.three_step)
#>              Intercept:C2        Zp:C2 Intercept:C3        Zp:C3
#> Intercept:C2  0.104798063 -0.041343485  -0.01048872  0.001510799
#> Zp:C2        -0.041343485  0.019772531   0.01121397 -0.001932282
#> Intercept:C3 -0.010488717  0.011213968   0.34964328 -0.082842095
#> Zp:C3         0.001510799 -0.001932282  -0.08284210  0.020597096
```

### Proportional assignment

With modal assignment (`use.modal.assignment = TRUE`, the default), the
Jacobian in the measurement-uncertainty correction is not mathematically
defined. Setting `use.modal.assignment = FALSE` uses soft posterior
weights throughout, giving an analytic Jacobian and is recommended when
separation is moderate or low. When `use.modal.assignment = TRUE`, the
Jacobian $`\frac{\partial\theta_2}{\partial\theta_1}`$ computed using
the full posterior weights (e.g., behaving as if
`use.modal.assignment = FALSE`) to maintain well-defined derivatives,
though three-step estimates would still be computed with modal
assignment as specified. The different is negligible when separation is
high.

``` r

d.three_step.prop <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
summary(d.three_step.prop)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0617
#> AIC            : 2758.1234
#> BIC            : 2926.7078
#> Entropy R²     : 0.8680  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.1096    0.3287  6.4183 < 0.001 ***
#> Zp:C2         -0.9121    0.1460 -6.2456 < 0.001 ***
#> Intercept:C3  -3.6508    0.6430 -5.6777 < 0.001 ***
#> Zp:C3          0.9367    0.1547  6.0563 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Simple (robust) standard errors

Setting `use.simple.cov = TRUE` skips the measurement-uncertainty
correction and returns the robust sandwich SEs from Step 3 only. When
separation is high the correction is negligible, so this is a useful
computational shortcut for large samples.

``` r

d.three_step.simple <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE
)
summary(d.three_step.simple)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0650
#> AIC            : 2758.1299
#> BIC            : 2926.7143
#> Entropy R²     : 0.8693  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.0411    0.3214  6.3504 < 0.001 ***
#> Zp:C2         -0.8821    0.1394 -6.3257 < 0.001 ***
#> Intercept:C3  -3.4836    0.5876 -5.9281 < 0.001 ***
#> Zp:C3          0.8985    0.1426  6.3008 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### BCH estimator

The BCH correction of Bolck, Croon & Hagenaars (2004) is available via
`use.bch = TRUE`. It works well with high separation but can produce an
ill-conditioned Hessian when separation is low (resulting in a
covariance matrix that is not positive semi-definite), in which case the
ML estimator is preferred.

``` r

d.three_step.bch <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.bch = TRUE
)
summary(d.three_step.bch)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : BCH
#> Log-likelihood : -1339.2863
#> AIC            : 2758.5726
#> BIC            : 2927.1569
#> Entropy R²     : 0.8700  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.9554    0.3111  6.2844 < 0.001 ***
#> Zp:C2         -0.8424    0.1304 -6.4613 < 0.001 ***
#> Intercept:C3  -3.4634    0.5697 -6.0790 < 0.001 ***
#> Zp:C3          0.8923    0.1385  6.4412 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

BCH with low-separation data can fail to produce a positive
semi-definite Hessian. The ML estimator with proportional assignment is
more reliable in this setting.

``` r

# Not run in vignette build (slow and and produces warnings)
bch.fail <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.bch = TRUE,
  maxIter.measurement = 2000,
  iter.measurement = 10
)
```

``` r

# Preferred approach for low separation
d.low.three_step.prop <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
summary(d.low.three_step.prop)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1979.3372
#> AIC            : 4038.6744
#> BIC            : 4207.2588
#> Entropy R²     : 0.3518  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  3.0267 -3.6919
#> Zp        -0.9767  0.9482
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   3.2034    2.2929  1.3971 0.1624     
#> Zp:C2         -1.0761    1.9088 -0.5638 0.5729     
#> Intercept:C3  -3.8431    2.9955 -1.2830 0.1995     
#> Zp:C3          0.9554    0.6034  1.5832 0.1134     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Choosing the reference class

By default, class 1 (`"C1"`) is the reference category for the
multinomial logit parameterization. The `rebase` argument changes this.
Estimates are reparameterized consistently: log-likelihoods are
invariant, and the coefficients satisfy the transitivity relation
$`\log(\pi_t / \pi_j) = \log(\pi_t / \pi_1) - \log(\pi_j / \pi_1)`$.

``` r

# Default: C1 as reference
summary(d.three_step.simple)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0650
#> AIC            : 2758.1299
#> BIC            : 2926.7143
#> Entropy R²     : 0.8693  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.0411    0.3214  6.3504 < 0.001 ***
#> Zp:C2         -0.8821    0.1394 -6.3257 < 0.001 ***
#> Intercept:C3  -3.4836    0.5876 -5.9281 < 0.001 ***
#> Zp:C3          0.8985    0.1426  6.3008 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r

d.three_step.simpleC2 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C2"
)
summary(d.three_step.simpleC2)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0650
#> AIC            : 2758.1299
#> BIC            : 2926.7143
#> Entropy R²     : 0.8693  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C1      C3
#> Intercept -2.1941 -5.6433
#> Zp         0.9413  1.8377
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C1  -2.0411    0.3214 -6.3504 < 0.001 ***
#> Zp:C1          0.8821    0.1394  6.3257 < 0.001 ***
#> Intercept:C3  -5.5246    0.6823 -8.0976 < 0.001 ***
#> Zp:C3          1.7806    0.2078  8.5696 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r

d.three_step.simpleC3 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C3"
)
summary(d.three_step.simpleC3)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0650
#> AIC            : 2758.1299
#> BIC            : 2926.7143
#> Entropy R²     : 0.8693  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C1      C2
#> Intercept  3.4492  5.6433
#> Zp        -0.8964 -1.8377
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C1   3.4836    0.5876  5.9281 < 0.001 ***
#> Zp:C1         -0.8985    0.1426 -6.3008 < 0.001 ***
#> Intercept:C2   5.5246    0.6823  8.0976 < 0.001 ***
#> Zp:C2         -1.7806    0.2078 -8.5696 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Passing a pre-fitted measurement model

The `step1` argument accepts any previously fitted `tseLCA` object or
the raw output of
[`lca_step1()`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md).
This is useful when you want to:

- Reuse an expensive measurement model across multiple structural
  models.
- Estimate the measurement model on a large reference sample and apply
  it to a smaller analysis sample.
- Inject custom two-step starting values computed via
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md).

``` r

# Reuse the measurement model estimated above
d.three_step.prop2 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.measurement$measurement_model
)
summary(d.three_step.prop2)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1339.0617
#> AIC            : 2758.1234
#> BIC            : 2926.7078
#> Entropy R²     : 0.8680  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.1934 -3.4524
#> Zp        -0.9411  0.8972
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.1096    0.3287  6.4183 < 0.001 ***
#> Zp:C2         -0.9121    0.1460 -6.2456 < 0.001 ***
#> Intercept:C3  -3.6508    0.6430 -5.6777 < 0.001 ***
#> Zp:C3          0.9367    0.1547  6.0563 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r

# Measurement model from a larger low-separation sample
d.low2000 <- generate_data(
  n = 2000,
  separation = "low",
  scenario = "covariate",
  seed = 2
)
d.low.measurement2000 <- three_step(
  data = d.low2000,
  Y.names = paste0("Y", 1:6),
  n_classes = 3
)

# Apply to the smaller sample; get.twostep.vcov returns multilevLCA's
# bias-corrected vcov for the two-step estimates
d.low.three_step.prop2 <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.low.measurement2000$measurement_model,
  get.twostep.vcov = TRUE
)
summary(d.low.three_step.prop2)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1983.8159
#> AIC            : 4047.6319
#> BIC            : 4216.2162
#> Entropy R²     : 0.3770  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.5856 -4.2916
#> Zp        -1.3548  1.0808
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.7061    1.1851  2.2835 0.0224  *  
#> Zp:C2         -1.3765    0.9785 -1.4068 0.1595     
#> Intercept:C3  -3.9495    2.0920 -1.8879 0.0590  .  
#> Zp:C3          1.0305    0.4740  2.1741 0.0297  *  
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

You can also compute two-step starting values separately and inject them
before calling
[`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

``` r

d.low.fitZ2 <- fitZ_from_fit0(
  fit0 = d.low.measurement2000$measurement_model$fit0,
  data = d.low,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp"
)
d.low.measurement2000$measurement_model$fitZ <- d.low.fitZ2

d.low.three_step.prop3 <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.low.measurement2000$measurement_model
)
summary(d.low.three_step.prop3)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1983.8159
#> AIC            : 4047.6319
#> BIC            : 4216.2162
#> Entropy R²     : 0.3770  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.5856 -4.2916
#> Zp        -1.3548  1.0808
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.7061    1.1851  2.2835 0.0224  *  
#> Zp:C2         -1.3765    0.9785 -1.4068 0.1595     
#> Intercept:C3  -3.9495    2.0920 -1.8879 0.0590  .  
#> Zp:C3          1.0305    0.4740  2.1741 0.0297  *  
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Missing data

`tseLCA` uses a two-pass row-filtering strategy that matches
`multilevLCA`’s approach for the measurement model while allowing more
observations into Steps 1 and 2 than Step 3.

``` r

set.seed(42)
d.new <- generate_data(500, separation = "high", seed = 3)
sparsity <- 0.1
missing <- 1 -
  matrix(
    rbinom(prod(dim(d.new)), size = 1, prob = sparsity),
    nrow = nrow(d.new),
    ncol = ncol(d.new)
  )
missing[missing == 0] <- NA_real_
d.sparse <- d.new * missing
head(d.sparse)
#>   Y1 Y2 Y3 Y4 Y5 Y6  X Zp
#> 1  0  0 NA  0  1  0  3  5
#> 2  1  1 NA  0  0  0  2  2
#> 3  1  1  1  1  1  0  1  4
#> 4  0  0  0  0  0  0  3  4
#> 5  1  1  1  0  0  0 NA  2
#> 6  1  1 NA  0  0  0  2  3
```

With `incomplete = FALSE` (the default), any row with a missing
indicator is dropped before the measurement model is estimated.

``` r

d.sparse.measurement <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = FALSE,
  verbose = TRUE
)
#> 242 row(s) dropped from measurement/classification steps (missing Y).
# Rows dropped = number of rows with at least one missing Y
sum(apply(d.sparse[, paste0("Y", 1:6)], 1, \(x) any(is.na(x))))
#> [1] 242
summary(d.sparse.measurement)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -742.0656
#> AIC            : 1524.1311
#> BIC            : 1595.1903
#> Entropy R²     : 0.9027
#> 
#> Class prevalences:
#>             
#> P(C1) 0.2995
#> P(C2) 0.3967
#> P(C3) 0.3037
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8241 0.8869 0.0834
#> P(Y2|C) 0.8600 0.9104 0.0811
#> P(Y3|C) 0.9110 0.9344 0.0633
#> P(Y4|C) 0.9035 0.0486 0.1302
#> P(Y5|C) 0.9349 0.1517 0.0795
#> P(Y6|C) 0.9331 0.1555 0.0762
```

With `incomplete = TRUE`, only fully-missing rows are dropped; partially
observed rows contribute to the measurement model via FIML.

``` r

d.sparse.measurement2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = TRUE,
  verbose = TRUE
)
summary(d.sparse.measurement2)
#> -- tseLCA Measurement Model --------------------------------
#> Latent classes : 3
#> Log-likelihood : -1342.8026
#> AIC            : 2725.6052
#> BIC            : 2809.8974
#> Entropy R²     : 0.8425
#> 
#> Class prevalences:
#>             
#> P(C1) 0.3049
#> P(C2) 0.3652
#> P(C3) 0.3299
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8797 0.8916 0.0925
#> P(Y2|C) 0.8888 0.8858 0.0677
#> P(Y3|C) 0.9337 0.8859 0.1453
#> P(Y4|C) 0.9079 0.0819 0.1339
#> P(Y5|C) 0.9536 0.1359 0.1056
#> P(Y6|C) 0.9583 0.1590 0.1176
```

Regardless of `incomplete`, Step 3 drops any row with a missing
covariate. The rows used in Step 3 are a subset of those used in Steps 1
and 2.

``` r

d.sparse.three_step <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  verbose = TRUE
)
#> 43 row(s) excluded from covariate step (missing Z).
#> fitZ EM converged in 9 iterations.
#> 43 row(s) excluded from covariate step (missing Z).
#> EM converged in 8 iterations.
# Additional rows dropped from Step 3 due to missing Zp
sum(is.na(d.sparse$Zp))
#> [1] 43
summary(d.sparse.three_step)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1088.3344
#> AIC            : 2256.6688
#> BIC            : 2421.6562
#> Entropy R²     : 0.8672  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.5023 -4.6555
#> Zp        -1.0494  1.2929
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.4803    0.3807  6.5149 < 0.001 ***
#> Zp:C2         -1.0136    0.1667 -6.0789 < 0.001 ***
#> Intercept:C3  -4.7754    0.7267 -6.5716 < 0.001 ***
#> Zp:C3          1.3164    0.1796  7.3290 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

A FIML measurement model can be passed in and then reused for the
covariate step on the same sparse data.

``` r

d.sparse.three_step2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  step1 = d.sparse.measurement2$measurement_model,
  verbose = TRUE
)
#> 43 row(s) excluded from covariate step (missing Z).
#> fitZ EM converged in 9 iterations.
#> 43 row(s) excluded from covariate step (missing Z).
#> EM converged in 8 iterations.
summary(d.sparse.three_step2)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1088.3344
#> AIC            : 2256.6688
#> BIC            : 2421.6562
#> Entropy R²     : 0.8672  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.5023 -4.6555
#> Zp        -1.0494  1.2929
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.4803    0.3807  6.5149 < 0.001 ***
#> Zp:C2         -1.0136    0.1667 -6.0789 < 0.001 ***
#> Intercept:C3  -4.7754    0.7267 -6.5716 < 0.001 ***
#> Zp:C3          1.3164    0.1796  7.3290 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Polytomous items

`tseLCA` supports polytomous indicators, following `multilevLCA`’s
convention that item categories are coded as consecutive integers
starting at 0.

Here we reproduce the example from the `poLCA` package.

``` r

data(election, package = "poLCA")
elec <- election
elec.items <- colnames(election)[1:12]

# Recode to 0-based integers as required by multilevLCA
elec[, elec.items] <- lapply(elec[, elec.items], \(x) as.integer(x) - 1L)
```

``` r

elec.measurement <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  #The poLCA example drops any row with a missing cell
  incomplete = FALSE
)

elec.three_step <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  Zp.names = c("PARTY"),
  step1 = elec.measurement$measurement_model,
  incomplete = FALSE,
  #With the neutral group as the base-category
  rebase = "C3"
)
#> Warning: lca_indiv_varmat: Infomat is singular even after removing boundary
#> parameters; returning NA matrix. Check for near-empty classes.
summary(elec.three_step)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -16278.0242
#> AIC            : 32852.0485
#> BIC            : 33617.2262
#> Entropy R²     : 0.7956  (covariate-adjusted)
#> 
#> Two-step (starting) estimates:
#>                C1      C2
#> Intercept -2.5781  1.8687
#> PARTY      0.4289 -0.6983
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value p.value
#> Intercept:C1  -2.4701        NA      NA      NA
#> PARTY:C1       0.4077        NA      NA      NA
#> Intercept:C2   1.7324        NA      NA      NA
#> PARTY:C2      -0.6727        NA      NA      NA
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r


party.x <- seq(from = 1, to = 7, length.out = 101)
pidmat <- cbind(1, party.x)
exb <- exp(pidmat %*% coef(elec.three_step))

matplot(
  party.x,
  (cbind(1, exb)) / (1 + rowSums(exb)),
  ylim = c(0, 1),
  type = "l",
  lwd = 3,
  col = 1,
  xlab = "Party ID: strong Democratic (1) to strong Republican (7)",
  ylab = "Probability of latent class membership",
  main = "Party ID as a predictor of candidate affinity class",
)
text(3.9, 0.60, "Other")
text(6.2, 0.6, "Bush affinity")
text(2.0, 0.65, "Gore affinity")
```

![](tseLCA-workflow_files/figure-html/elec-example-1.png)

------------------------------------------------------------------------

## Distal outcomes

For distal outcomes ($`Z_o \leftarrow X \rightarrow Y`$), supply
`Zo.name` and a `family` argument. The available families are
`"gaussian"` (default), `"poisson"`, and `"binomial"`. Both ML and BCH
estimators are available.

``` r

d.distal <- generate_data(
  n = 500,
  separation = "high",
  scenario = "distal",
  seed = 4
)
# True class means: mu = (0, 1, -1) for C1, C2, C3
```

``` r

d.distal.measurement <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3
)

# ML estimator
d.distal.three_step.ml <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  use.modal.assignment = FALSE,
  family = "gaussian"
)

# BCH estimator: closed-form M-step for distal outcomes
d.distal.three_step.bch <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  use.modal.assignment = FALSE,
  use.bch = TRUE,
  family = "gaussian"
)

summary(d.distal.three_step.ml)
#> -- tseLCA Three-Step Distal Outcome Model -------------------
#> Latent classes : 3
#> Estimator      : ML
#> Family         : gaussian
#> Log-likelihood : -2169.0110
#> AIC            : 4384.0220
#> BIC            : 4480.9580
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -1.0821    0.0817 -13.2495 < 0.001 ***
#> mu_C2 (mean)   1.0172    0.0819  12.4151 < 0.001 ***
#> mu_C3 (mean)   0.0254    0.0878   0.2887 0.7728     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
summary(d.distal.three_step.bch)
#> -- tseLCA Three-Step Distal Outcome Model -------------------
#> Latent classes : 3
#> Estimator      : BCH
#> Family         : gaussian
#> Log-likelihood : -2168.8685
#> AIC            : 4383.7370
#> BIC            : 4480.6730
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -1.0941    0.0867 -12.6134 < 0.001 ***
#> mu_C2 (mean)   0.9751    0.0823  11.8534 < 0.001 ***
#> mu_C3 (mean)   0.0578    0.0859   0.6733 0.5008     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Three-step estimation with both covariates (Zp) and distal outcomes (Zo)

Consistent with how most research in the social sciences construct the
relationships between $`Z_p`$ and $`X`$, and $`X`$ and $`Z_o`$, the
relationship between $`Z_p`$ and $`X`$ is estimated **first**, followed
by estimation between $`X`$ and $`Z_o`$, adjusting for the
covariate-adjusted posteriors in the estimation procedures for the
distal outcome model in step 3.

``` r

d.covariate <- generate_data(
  n = 500,
  separation = "high",
  scenario = "covariate",
  seed = 4
)
d.covariate$Zo <- draw_Zo(d.covariate$X, bk2018_params$distal_params)
head(d.covariate)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp         Zo
#> 1  1  1  1  1  0  0 2  3 -0.1624650
#> 2  1  1  1  1  1  1 1  3 -1.1591833
#> 3  1  1  1  1  1  1 1  3 -1.2055132
#> 4  0  0  0  0  0  0 3  4  1.8752276
#> 5  1  1  1  1  1  1 1  3 -2.5582369
#> 6  1  1  1  1  1  1 1  5 -0.4723262

d.covariate.three_step <- three_step(
  data = d.covariate,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  Zo.name = "Zo",
  use.modal.assignment = FALSE
)
summary(d.covariate.three_step)
#> -- tseLCA Three-Step Model: Covariate + Distal Outcome -----
#> Latent classes : 3
#> Estimator      : ML
#> Family         : gaussian
#> Log-likelihood : -1315.6596
#> AIC            : 2711.3193
#> BIC            : 2879.9036
#> 
#> Covariate -- two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.4973 -4.2196
#> Zp        -1.0177  1.1159
#> 
#> Covariate -- three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.6602    0.3998  6.6538 < 0.001 ***
#> Zp:C2         -1.0790    0.1632 -6.6134 < 0.001 ***
#> Intercept:C3  -4.6917    0.7070 -6.6365 < 0.001 ***
#> Zp:C3          1.2278    0.1735  7.0773 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Distal outcome -- three-step estimates:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -0.9851    0.0828 -11.8995 < 0.001 ***
#> mu_C2 (mean)   0.9298    0.0851  10.9217 < 0.001 ***
#> mu_C3 (mean)   0.1188    0.0722   1.6458 0.0998  .  
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Note that with covariates in a model with high separation, the standard
errors above should, on average, by systematically smaller for distal
outcome estimation than if there were no covariates in the model (see
below).

``` r

three_step(
  data = d.covariate,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  use.modal.assignment = FALSE
) |>
  vcov() |>
  diag() |>
  sqrt()
#>      mu_C1      mu_C2      mu_C3 
#> 0.08302639 0.08922331 0.07595896
```

------------------------------------------------------------------------

## References

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311. <https://doi.org/10.1177/0081175012470644>

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
<https://doi.org/10.1007/s11336-017-9592-7>

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, 12(1), 3–27.
<https://doi.org/10.1093/pan/mph001>

Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025).
Multilevel latent class analysis: State-of-the-art methodologies and
their implementation in the R package multilevLCA. *Multivariate
Behavioral Research*, 60(4), 731–747.
<https://doi.org/10.1080/00273171.2025.2473935>

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
<https://doi.org/10.1093/pan/mpq025>

------------------------------------------------------------------------

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tseLCA_1.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10        generics_0.1.4     tidyr_1.3.2        pracma_2.4.6      
#>  [5] hms_1.1.4          digest_0.6.39      magrittr_2.0.5     RColorBrewer_1.1-3
#>  [9] evaluate_1.0.5     iterators_1.0.14   fastmap_1.2.0      foreach_1.5.2     
#> [13] jsonlite_2.0.0     combinat_0.0-8     promises_1.5.0     purrr_1.2.2       
#> [17] codetools_0.2-20   textshaping_1.0.5  jquerylib_0.1.4    cli_3.6.6         
#> [21] shiny_1.14.0       labelled_2.16.0    rlang_1.2.0        cachem_1.1.0      
#> [25] yaml_2.3.12        otel_0.2.0         klaR_1.7-4         parallel_4.6.1    
#> [29] tools_4.6.1        dplyr_1.2.1        httpuv_1.6.17      forcats_1.0.1     
#> [33] vctrs_0.7.3        R6_2.6.1           mime_0.13          lifecycle_1.0.5   
#> [37] multilevLCA_2.1.4  tictoc_1.2.1       fs_2.1.0           MASS_7.3-65       
#> [41] miniUI_0.1.2       cluster_2.1.8.2    ragg_1.5.2         pkgconfig_2.0.3   
#> [45] desc_1.4.3         pkgdown_2.2.0      bslib_0.11.0       pillar_1.11.1     
#> [49] later_1.4.8        glue_1.8.1         Rcpp_1.1.1-1.1     systemfonts_1.3.2 
#> [53] haven_2.5.5        xfun_0.59          tibble_3.3.1       tidyselect_1.2.1  
#> [57] highr_0.12         rstudioapi_0.19.0  knitr_1.51         xtable_1.8-8      
#> [61] htmltools_0.5.9    rmarkdown_2.31     clustMixType_0.5-1 compiler_4.6.1    
#> [65] questionr_0.8.2
```
