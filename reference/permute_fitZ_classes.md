# Permute class columns of a fitZ object to match a new reference class

Rebases a `fitZ` object (output of `fitZ_from_fit0` or
`fitZ_from_multiLCA`) so that `ref_idx` becomes the reference class.
This involves:

1.  Rebasing `$mGamma` – reconstructing the full T-column log-ratio
    matrix, subtracting the new reference column, and dropping it.

2.  Propagating through `$Varmat_cor` via the delta method – the
    rebasing transformation is linear (`gamma_new = A * gamma_old`) so
    the vcov transforms exactly as `A %*% V %*% t(A)`.

3.  Updating all column names.

## Usage

``` r
permute_fitZ_classes(fitZ, ref_idx)
```

## Arguments

- fitZ:

  Output of
  [`fitZ_from_fit0()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md)
  or
  [`fitZ_from_multiLCA()`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md).

- ref_idx:

  Integer. New reference class (1-based index into the T classes as
  currently ordered in `fitZ`).

## Value

`fitZ` with `$mGamma`, `$Varmat_cor`, and names updated.
