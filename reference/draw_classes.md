# Draw latent class memberships from their marginal distribution

Draw latent class memberships from their marginal distribution

## Usage

``` r
draw_classes(n, pi_)
```

## Arguments

- n:

  Integer. Sample size.

- pi\_:

  Numeric vector of length T. Class proportions (must sum to 1).

## Value

Integer vector of length n with values in `1:T`.

## Examples

``` r
# Draw 100 class labels from equal prevalences
draw_classes(100, c(1/3, 1/3, 1/3))
#>   [1] 2 1 3 2 3 1 3 1 2 3 3 2 3 1 1 3 3 3 1 3 3 3 3 1 1 1 2 2 1 1 2 3 3 1 1 2 2
#>  [38] 1 2 3 3 3 1 3 3 3 3 2 2 3 2 3 3 1 2 3 1 3 3 3 2 3 3 2 3 2 2 2 1 2 2 1 2 1
#>  [75] 1 2 1 3 3 1 3 2 2 3 1 3 3 3 3 3 1 1 1 1 2 1 1 2 2 1
```
