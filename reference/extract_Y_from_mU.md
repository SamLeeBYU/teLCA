# Extract Y.exp, mDesign, posteriors from a multilevLCA mU matrix

`fit0$mU` from multilevLCA stores data already in one-hot expanded form:
each item h occupies K_h consecutive columns (one per category),
followed by T columns of posterior class probabilities.

## Usage

``` r
extract_Y_from_mU(fit0, ivItemcat = NULL)
```

## Arguments

- fit0:

  Raw multilevLCA fit object with `$mU`, `$mPhi`, `$vPi`.

- ivItemcat:

  Integer vector of category counts per item (length H). If `NULL`,
  inferred from `fit0$mPhi` dimensions.

## Value

A list with:

- Y.exp:

  N x K_total expanded one-hot matrix (NAs replaced with 0).

- mDesign:

  N x K_total design/mask matrix. `NULL` if no missing.

- ivItemcat:

  Integer vector of category counts per item.

- u_post:

  N x T posterior class probability matrix from `mU`.

## Details

For dichotomous items (K_h=2) the two columns are stored. For polytomous
items (K_h\>2) all K_h columns are stored. This function first
compresses the expanded Y back to integer codes through `compress_Y`,
then re-expands consistently via `expand_Y` so downstream functions
receive the correct N x K_total matrix.
