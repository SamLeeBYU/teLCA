# Draw latent classes conditional on the covariate (scenario "covariate")

Draw latent classes conditional on the covariate (scenario "covariate")

## Usage

``` r
draw_classes_given_Zp(Zp, params)
```

## Arguments

- Zp:

  Numeric vector of length n. Covariate values.

- params:

  Multinomial logistic parameter list (see
  `bk2018_params$covariate_params`).

## Value

Integer vector of length n with class labels in `1:T`.

## Examples

``` r
Zp <- draw_Zp(1000)
X  <- draw_classes_given_Zp(Zp, bk2018_params$covariate_params)
table(X) # Should be roughly uniform
#> X
#>   1   2   3 
#> 294 341 365 
```
