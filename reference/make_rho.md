# Build the item-response probability matrix for the simulation

Returns a T x K matrix `rho` where `rho[t, k] = P(Y_k = 1 | X = t)`.

## Usage

``` r
make_rho(pi_)
```

## Arguments

- pi\_:

  Numeric scalar in (0.5, 1). Probability of the "likely" response. Use
  `bk2018_params$separation_levels` for the three simulation levels.

## Value

A 3 x 6 numeric matrix.

## Details

The three-class structure is:

- Class 1: `pi_` on all 6 items (high responders).

- Class 2: `pi_` on items 1-3, `1 - pi_` on items 4-6 (mixed).

- Class 3: `1 - pi_` on all 6 items (low responders).

## Examples

``` r
# High separation: P(Y=1|class) = 0.9 for the "high" class
make_rho(0.9)
#>        [,1] [,2] [,3] [,4] [,5] [,6]
#> class1  0.9  0.9  0.9  0.9  0.9  0.9
#> class2  0.9  0.9  0.9  0.1  0.1  0.1
#> class3  0.1  0.1  0.1  0.1  0.1  0.1

# Low separation
make_rho(0.7)
#>        [,1] [,2] [,3] [,4] [,5] [,6]
#> class1  0.7  0.7  0.7  0.7  0.7  0.7
#> class2  0.7  0.7  0.7  0.3  0.3  0.3
#> class3  0.3  0.3  0.3  0.3  0.3  0.3
```
