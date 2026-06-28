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
#> 3  1  1  1  1  1  1 1  1
#> 4  1  1  1  1  1  0 2  2
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
#> 1  1  0  0  0  0  0 2  1
#> 2  0  0  0  1  1  0 3  4
#> 3  1  1  1  1  1  1 1  1
#> 4  1  1  1  1  1  1 2  2
#> 5  0  0  0  0  0  0 3  5
#> 6  1  1  1  0  1  1 1  3
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
#> Log-likelihood : -1467.9716
#> AIC            : 2975.9433
#> BIC            : 3060.2354
#> Entropy R²     : 0.8745
#> 
#> Class prevalences:
#>             
#> P(C1) 0.3468
#> P(C2) 0.2950
#> P(C3) 0.3582
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8920 0.8682 0.1187
#> P(Y2|C) 0.9340 0.8954 0.1283
#> P(Y3|C) 0.8893 0.8915 0.1077
#> P(Y4|C) 0.8473 0.0680 0.1007
#> P(Y5|C) 0.9131 0.0773 0.0908
#> P(Y6|C) 0.9018 0.0882 0.1040
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
#> Log-likelihood : -2006.5242
#> AIC            : 4053.0483
#> BIC            : 4137.3405
#> Entropy R²     : 0.4447
#> 
#> Class prevalences:
#>             
#> P(C1) 0.4551
#> P(C2) 0.1723
#> P(C3) 0.3726
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.6411 0.7435 0.3010
#> P(Y2|C) 0.8222 0.7303 0.2255
#> P(Y3|C) 0.6894 0.4495 0.3511
#> P(Y4|C) 0.6122 0.2426 0.2674
#> P(Y5|C) 0.6651 0.0661 0.3165
#> P(Y6|C) 0.6540 0.1105 0.3632
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
#>                   C2        C3
#> Intercept  1.7790263 -4.170998
#> Zp        -0.8401332  1.148593
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
#>                   C2        C3
#> Intercept  0.1647684 -3.396233
#> Zp        -0.5668117  0.937924
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
#> Log-likelihood : -1338.2884
#> AIC            : 2756.5768
#> BIC            : 2925.1611
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.8739    0.3159  5.9313 < 0.001 ***
#> Zp:C2         -0.8833    0.1331 -6.6364 < 0.001 ***
#> Intercept:C3  -4.3703    0.6918 -6.3169 < 0.001 ***
#> Zp:C3          1.1945    0.1735  6.8850 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The standard [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) S3 methods work on any
`tseLCA` object.

``` r

coef(d.three_step)
#>                   C2        C3
#> Intercept  1.8738903 -4.370297
#> Zp        -0.8833431  1.194500
vcov(d.three_step)
#>              Intercept:C2        Zp:C2 Intercept:C3        Zp:C3
#> Intercept:C2  0.099811915 -0.037977145  -0.02047535  0.003563972
#> Zp:C2        -0.037977145  0.017717066   0.01161711 -0.001880094
#> Intercept:C3 -0.020475347  0.011617106   0.47864889 -0.117651101
#> Zp:C3         0.003563972 -0.001880094  -0.11765110  0.030100056
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
#> Log-likelihood : -1338.3943
#> AIC            : 2756.7886
#> BIC            : 2925.3730
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.8671    0.3189  5.8552 < 0.001 ***
#> Zp:C2         -0.8756    0.1352 -6.4764 < 0.001 ***
#> Intercept:C3  -4.5000    0.6967 -6.4591 < 0.001 ***
#> Zp:C3          1.2258    0.1752  6.9950 < 0.001 ***
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
#> Log-likelihood : -1338.2884
#> AIC            : 2756.5768
#> BIC            : 2925.1611
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.8739    0.3145  5.9578 < 0.001 ***
#> Zp:C2         -0.8833    0.1323 -6.6760 < 0.001 ***
#> Intercept:C3  -4.3703    0.6863 -6.3681 < 0.001 ***
#> Zp:C3          1.1945    0.1722  6.9378 < 0.001 ***
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
#> Log-likelihood : -1338.4071
#> AIC            : 2756.8142
#> BIC            : 2925.3985
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.9274    0.3346  5.7596 < 0.001 ***
#> Zp:C2         -0.9084    0.1468 -6.1880 < 0.001 ***
#> Intercept:C3  -3.9789    0.6458 -6.1609 < 0.001 ***
#> Zp:C3          1.1066    0.1636  6.7662 < 0.001 ***
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
#> Log-likelihood : -1974.5489
#> AIC            : 4029.0977
#> BIC            : 4197.6821
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  0.1594 -3.3877
#> Zp        -0.5637  0.9358
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   0.9140    1.7455  0.5237 0.6005     
#> Zp:C2         -1.0336    1.6784 -0.6158 0.5380     
#> Intercept:C3  -4.2276    1.7326 -2.4400 0.0147  *  
#> Zp:C3          1.1703    0.3835  3.0515 0.0023  ** 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Choosing the reference class

By default, class 1 (`"C1"`) is the reference category for the
multinomial logit parameterisation. The `rebase` argument changes this.
Estimates are reparameterised consistently: log-likelihoods are
invariant, and the coefficients satisfy the transitivity relation
$`\log(\pi_t / \pi_j) = \log(\pi_t / \pi_1) - \log(\pi_j / \pi_1)`$.

``` r

# Default: C1 as reference
summary(d.three_step.simple)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1338.2884
#> AIC            : 2756.5768
#> BIC            : 2925.1611
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.8739    0.3145  5.9578 < 0.001 ***
#> Zp:C2         -0.8833    0.1323 -6.6760 < 0.001 ***
#> Intercept:C3  -4.3703    0.6863 -6.3681 < 0.001 ***
#> Zp:C3          1.1945    0.1722  6.9378 < 0.001 ***
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
#> Log-likelihood : -1338.2884
#> AIC            : 2756.5768
#> BIC            : 2925.1611
#> 
#> Two-step (starting) estimates:
#>                C1      C3
#> Intercept -1.7786 -5.9494
#> Zp         0.8400  1.9884
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C1  -1.8739    0.3145 -5.9578 < 0.001 ***
#> Zp:C1          0.8833    0.1323  6.6760 < 0.001 ***
#> Intercept:C3  -6.2442    0.7785 -8.0213 < 0.001 ***
#> Zp:C3          2.0778    0.2245  9.2571 < 0.001 ***
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
#> Log-likelihood : -1338.2884
#> AIC            : 2756.5768
#> BIC            : 2925.1611
#> 
#> Two-step (starting) estimates:
#>                C1      C2
#> Intercept  4.1707  5.9494
#> Zp        -1.1484 -1.9884
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C1   4.3703    0.6863  6.3681 < 0.001 ***
#> Zp:C1         -1.1945    0.1722 -6.9378 < 0.001 ***
#> Intercept:C2   6.2442    0.7785  8.0213 < 0.001 ***
#> Zp:C2         -2.0778    0.2245 -9.2571 < 0.001 ***
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
#> Log-likelihood : -1338.3943
#> AIC            : 2756.7886
#> BIC            : 2925.3730
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.7790 -4.1710
#> Zp        -0.8401  1.1486
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.8671    0.3189  5.8552 < 0.001 ***
#> Zp:C2         -0.8756    0.1352 -6.4764 < 0.001 ***
#> Intercept:C3  -4.5000    0.6967 -6.4591 < 0.001 ***
#> Zp:C3          1.2258    0.1752  6.9950 < 0.001 ***
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
#> Log-likelihood : -1986.2340
#> AIC            : 4052.4680
#> BIC            : 4221.0523
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.3687 -4.0525
#> Zp        -0.3242  1.0992
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.3449    1.1548  1.1646 0.2442     
#> Zp:C2         -0.2669    0.5540 -0.4817 0.6300     
#> Intercept:C3  -4.9820    1.8208 -2.7361 0.0062  ** 
#> Zp:C3          1.2694    0.3157  4.0206 < 0.001 ***
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
#> Log-likelihood : -1986.2340
#> AIC            : 4052.4680
#> BIC            : 4221.0523
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  1.3687 -4.0525
#> Zp        -0.3242  1.0992
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   1.3449    1.1548  1.1646 0.2442     
#> Zp:C2         -0.2669    0.5540 -0.4817 0.6300     
#> Intercept:C3  -4.9820    1.8208 -2.7361 0.0062  ** 
#> Zp:C3          1.2694    0.3157  4.0206 < 0.001 ***
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
#> 1  0  0 NA  0  0  0  3  5
#> 2  1  1 NA  1  1  1  1  2
#> 3  0  0  1  0  0  1  3  4
#> 4  1  1  1  1  1  1  1  4
#> 5  0  1  1  1  1  1 NA  2
#> 6  1  0 NA  0  0  0  3  3
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
#> Log-likelihood : -763.0984
#> AIC            : 1566.1968
#> BIC            : 1637.2559
#> Entropy R²     : 0.8835
#> 
#> Class prevalences:
#>             
#> P(C1) 0.2914
#> P(C2) 0.3924
#> P(C3) 0.3162
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.9038 0.9015 0.1078
#> P(Y2|C) 0.8920 0.9116 0.0938
#> P(Y3|C) 0.8558 0.8980 0.1073
#> P(Y4|C) 0.9546 0.0841 0.0945
#> P(Y5|C) 0.8999 0.1121 0.1224
#> P(Y6|C) 0.9116 0.1262 0.1677
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
#> Log-likelihood : -1372.3013
#> AIC            : 2784.6025
#> BIC            : 2868.8947
#> Entropy R²     : 0.8376
#> 
#> Class prevalences:
#>             
#> P(C1) 0.3064
#> P(C2) 0.3764
#> P(C3) 0.3171
#> attr(,"names")
#> [1] "C1" "C2" "C3"
#> 
#> Item-response probabilities (P(Y=1|class)):
#>             C1     C2     C3
#> P(Y1|C) 0.8697 0.8820 0.0864
#> P(Y2|C) 0.8965 0.8784 0.1405
#> P(Y3|C) 0.8527 0.8662 0.0796
#> P(Y4|C) 0.9424 0.0817 0.0902
#> P(Y5|C) 0.9159 0.1019 0.1156
#> P(Y6|C) 0.9106 0.1184 0.1592
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
#> EM converged in 6 iterations.
# Additional rows dropped from Step 3 due to missing Zp
sum(is.na(d.sparse$Zp))
#> [1] 43
summary(d.sparse.three_step)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1126.5227
#> AIC            : 2333.0454
#> BIC            : 2498.0327
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.6367 -4.1123
#> Zp        -1.0765  1.0957
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.6060    0.3945  6.6066 < 0.001 ***
#> Zp:C2         -1.0627    0.1670 -6.3627 < 0.001 ***
#> Intercept:C3  -4.2255    0.5606 -7.5373 < 0.001 ***
#> Zp:C3          1.1247    0.1436  7.8307 < 0.001 ***
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
#> EM converged in 6 iterations.
summary(d.sparse.three_step2)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -1126.5227
#> AIC            : 2333.0454
#> BIC            : 2498.0327
#> 
#> Two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.6367 -4.1123
#> Zp        -1.0765  1.0957
#> 
#> Three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.6060    0.3945  6.6066 < 0.001 ***
#> Zp:C2         -1.0627    0.1670 -6.3627 < 0.001 ***
#> Intercept:C3  -4.2255    0.5606 -7.5373 < 0.001 ***
#> Zp:C3          1.1247    0.1436  7.8307 < 0.001 ***
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
summary(elec.three_step)
#> -- tseLCA Three-Step Covariate Model -----------------------
#> Latent classes : 3
#> Estimator      : ML
#> Log-likelihood : -16278.0242
#> AIC            : 32852.0485
#> BIC            : 33617.2262
#> 
#> Two-step (starting) estimates:
#>                C1      C2
#> Intercept -2.5781  1.8687
#> PARTY      0.4289 -0.6983
#> 
#> Three-step estimates:
#>              Estimate Std.Error  z.value     p.value
#> Intercept:C1  -2.4701    0.3098  -7.9728 < 0.001 ***
#> PARTY:C1       0.4077    0.0564   7.2248 < 0.001 ***
#> Intercept:C2   1.7324    0.1671  10.3703 < 0.001 ***
#> PARTY:C2      -0.6727    0.0551 -12.2042 < 0.001 ***
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
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -1.0899    0.0795 -13.7040 < 0.001 ***
#> mu_C2 (mean)   1.0012    0.0818  12.2438 < 0.001 ***
#> mu_C3 (mean)   0.0450    0.0856   0.5257 0.5991     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
summary(d.distal.three_step.bch)
#> -- tseLCA Three-Step Distal Outcome Model -------------------
#> Latent classes : 3
#> Estimator      : BCH
#> Family         : gaussian
#> 
#> Distal outcome estimates by class:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -1.0906    0.0862 -12.6567 < 0.001 ***
#> mu_C2 (mean)   0.9291    0.0805  11.5410 < 0.001 ***
#> mu_C3 (mean)   0.0592    0.0837   0.7066 0.4798     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

## Three-step estimation with both covariates ($`Z_p`$) and distal outcomes ($`Z_o`$)

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
#> 1  0  0  0  0  0  0 3  3 -1.1624650
#> 2  1  1  1  1  0  0 1  3 -1.1591833
#> 3  1  0  1  1  1  1 1  3 -1.2055132
#> 4  1  1  1  1  1  1 1  4  0.8752276
#> 5  0  1  1  0  1  0 1  3 -2.5582369
#> 6  0  0  0  0  1  0 3  5  0.5276738

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
#> Log-likelihood : -1300.1948
#> AIC            : 2680.3896
#> BIC            : 2848.9739
#> 
#> Covariate -- two-step (starting) estimates:
#>                C2      C3
#> Intercept  2.4310 -4.4953
#> Zp        -1.1027  1.2103
#> 
#> Covariate -- three-step estimates:
#>              Estimate Std.Error z.value     p.value
#> Intercept:C2   2.4602    0.3855  6.3819 < 0.001 ***
#> Zp:C2         -1.1178    0.1771 -6.3119 < 0.001 ***
#> Intercept:C3  -4.6274    0.6030 -7.6738 < 0.001 ***
#> Zp:C3          1.2413    0.1499  8.2791 < 0.001 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Distal outcome -- three-step estimates:
#>              Estimate Std.Error  z.value     p.value
#> mu_C1 (mean)  -0.9346    0.0774 -12.0711 < 0.001 ***
#> mu_C2 (mean)   1.0685    0.0894  11.9577 < 0.001 ***
#> mu_C3 (mean)  -0.0007    0.0745  -0.0097 0.9923     
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Note that with covariates in the model, the standard errors above are
systematically smaller for distal outcome estimation than if there were
no covariates in the model (see below).

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
#> 0.07688217 0.09796084 0.07614852
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
#> [1] tseLCA_0.1.4
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
