# Prepare and validate data for tseLCA estimation

Prepare and validate data for tseLCA estimation

## Usage

``` r
clean_data(
  data,
  Y.names,
  Zp.names = NULL,
  Zo.name = NULL,
  incomplete = FALSE,
  include.intercept = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- Zp.names:

  Character vector of covariate column names, or `NULL`.

- Zo.name:

  Single distal outcome column name, or `NULL`.

- incomplete:

  Logical. If `TRUE`, use FIML for partially-observed Y.

- include.intercept:

  Logical. Prepend intercept column to Z.

- verbose:

  Logical. Print row-drop messages.

## Value

A named list with:

- Y.obs:

  N_Y x K expanded one-hot indicator matrix for Steps 1 & 2.

- mDesign:

  N_Y x K design/mask matrix (all 1s when incomplete = FALSE).

- ivItemcat:

  Integer vector of category counts per item.

- keep_Y:

  Integer indices of rows kept for Steps 1 & 2 (into original N).

- Z_mat:

  N_Z x Q covariate design matrix, or NULL.

- keep_step3_Z_in_Y:

  Positions of Z-complete rows within keep_Y.

- Zo_mat:

  N_Zo x 1 distal outcome matrix, or NULL.

- keep_step3_Zo_in_Y:

  Positions of Zo-complete rows within keep_Y.
