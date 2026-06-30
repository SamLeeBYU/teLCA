# Draw a continuous distal outcome given true class memberships (scenario "distal")

Draw a continuous distal outcome given true class memberships (scenario
"distal")

## Usage

``` r
draw_Zo(X, params)
```

## Arguments

- X:

  Integer vector of length n. Latent class (1-indexed).

- params:

  List with `$mu` (length-T class means) and `$sigma` (SD). See
  `bk2018_params$distal_params`.

## Value

Numeric vector of length n.

## Examples

``` r
X  <- draw_classes(100, c(1/3, 1/3, 1/3))
Zo <- draw_Zo(X, bk2018_params$distal_params)
tapply(Zo, X, mean)   # should be close to true mu
#>         1         2         3 
#> -1.054481  1.211078 -0.040895 
```
