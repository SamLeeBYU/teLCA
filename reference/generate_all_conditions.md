# Generate datasets for all 18 conditions in the simulation design

Iterates over the 2 scenarios x 3 separation levels x 3 sample sizes,
generating `n_rep` independent replications per condition. Seeds are
derived deterministically from `base_seed` so the entire experiment is
reproducible from a single integer.

## Usage

``` r
generate_all_conditions(
  n_rep = 500L,
  base_seed = 5262026L,
  params = bk2018_params,
  scenarios = c("covariate", "distal"),
  sep_levels = c("low", "mid", "high"),
  sample_sizes = c(500L, 1000L, 2000L),
  verbose = TRUE
)
```

## Arguments

- n_rep:

  Integer. Replications per condition (paper uses 500).

- base_seed:

  Integer. Base seed for reproducibility.

- params:

  Population parameters list. Defaults to
  [bk2018_params](https://samleebyu.github.io/tseLCA/reference/bk2018_params.md).

- scenarios:

  Character. Lists the scenario(s) ("covariate" and/or "distal") wanting
  to be simulated. Passed into
  [`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md).

- sep_levels:

  Character. Lists the separation level(s) ("low", "mid", "high")
  wanting to be simulated. Passed into
  [`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md).

- sample_sizes:

  Integer. Lists the sample size(s) wanting to be generated for each
  replication condition. Passed into
  [`generate_data()`](https://samleebyu.github.io/tseLCA/reference/generate_data.md).

- verbose:

  Logical. If `TRUE` (default), display a live CLI progress bar with
  per-rep status and ETA.

## Value

Nested list indexed as
`datasets[[scenario]][[separation]][[as.character(n)]]`, each element a
list of `n_rep` data frames.

## Examples

``` r
# \donttest{
# Generate 5 replicates for mid and high separation only
datasets <- generate_all_conditions(n_rep = 5L, base_seed = 1L,
                                    sep_levels = c("mid", "high"))
# Access a single replicate
head(datasets[["covariate"]][["high"]][["500"]][[1]])
#>   Y1 Y2 Y3 Y4 Y5 Y6 X Zp
#> 1  1  1  1  1  0  1 1  2
#> 2  0  1  1  1  0  0 1  3
#> 3  1  1  1  1  1  1 1  3
#> 4  1  1  1  1  1  1 1  3
#> 5  1  1  1  0  0  0 2  1
#> 6  0  1  1  0  0  0 2  2
# }
```
