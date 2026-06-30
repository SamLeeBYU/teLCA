# Draw the covariate Zp ~ Uniform{1, 2, 3, 4, 5}

Draw the covariate Zp ~ Uniform{1, 2, 3, 4, 5}

## Usage

``` r
draw_Zp(n)
```

## Arguments

- n:

  Integer. Sample size.

## Value

Integer vector of length n.

## Examples

``` r
Zp <- draw_Zp(100)
table(Zp)
#> Zp
#>  1  2  3  4  5 
#> 21 19 17 23 20 
```
