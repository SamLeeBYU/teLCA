# Normalise row/column names of a fitZ\$mGamma matrix

A plain `multiLCA` object uses `rownames` like `"gamma(Intercept|C)"`
and `"gamma(Zp|C)"`. This function strips the `gamma(...)` wrapper so
names match the clean format used throughout tseLCA (`"Intercept"`,
`"Zp"`, etc.) and ensures column names are `"C2"`, `"C3"`, etc.

## Usage

``` r
normalise_fitZ_names(fitZ, Zp.names = NULL, n_classes = NULL)
```

## Arguments

- fitZ:

  A fitZ-like list with at least `$mGamma`.

- Zp.names:

  Character vector of covariate column names (used to set clean rownames
  when the raw names can't be parsed). If `NULL`, rownames are stripped
  from the `gamma(X|C)` pattern only.

- n_classes:

  Integer. Total number of classes (used to derive clean column names if
  they are non-standard).

## Value

`fitZ` with normalised `$mGamma` row/col names.
