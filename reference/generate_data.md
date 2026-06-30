# Generate one dataset following the Bakk & Kuha (2018) simulation design

Generate one dataset following the Bakk & Kuha (2018) simulation design

## Usage

``` r
generate_data(
  n,
  separation = c("low", "mid", "high"),
  scenario = c("covariate", "distal"),
  params = bk2018_params,
  seed = NULL
)
```

## Arguments

- n:

  Integer. Sample size (paper uses 500, 1000, or 2000).

- separation:

  Character. One of `"low"`, `"mid"`, `"high"`. Maps to pi = 0.70, 0.80,
  0.90 respectively.

- scenario:

  Character. One of:

  `"covariate"`

  :   Zp (discrete, 1-5) predicts latent X via multinomial logit.

  `"distal"`

  :   Latent X predicts continuous Zo via linear regression.

- params:

  List of population parameters. Defaults to
  [bk2018_params](https://samleebyu.github.io/tseLCA/reference/bk2018_params.md).

- seed:

  Integer or `NULL`. Optional random seed for reproducibility.

## Value

A `data.frame` with columns:

- `Y1` .. `Y6`:

  Binary indicators (always present).

- `X`:

  True latent class, integer 1-3 (not observed in practice).

- `Zp`:

  Integer covariate 1-5 (scenario `"covariate"` only).

- `Zo`:

  Continuous distal outcome (scenario `"distal"` only).

## Examples

``` r
# Covariate scenario with high separation
d <- generate_data(n = 200, separation = "high", scenario = "covariate",
                   seed = 1)
head(d)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp
#> 1  1  1  1  1  1  1 1  1
#> 2  1  1  1  1  1  1 1  4
#> 3  1  1  1  0  0  0 2  1
#> 4  1  1  1  0  0  0 2  2
#> 5  0  0  0  0  0  0 3  5
#> 6  1  1  0  1  1  1 1  3
colMeans(d)
#>    Y1    Y2    Y3    Y4    Y5    Y6     X    Zp 
#> 0.580 0.610 0.585 0.350 0.350 0.360 2.035 2.925 

# Distal outcome scenario
d2 <- generate_data(n = 200, separation = "high", scenario = "distal",
                    seed = 2)
head(d2)
#>   Y1 Y2 Y3 Y4 Y5 Y6 X         Zo
#> 1  1  1  1  0  0  0 2  2.2387443
#> 2  1  1  1  1  1  1 1 -0.7681038
#> 3  0  0  0  0  0  0 3 -0.3144379
#> 4  0  1  1  0  0  0 2  2.4997037
#> 5  1  0  1  1  1  1 1 -0.9304256
#> 6  1  1  1  1  1  1 1  0.3340337
```
