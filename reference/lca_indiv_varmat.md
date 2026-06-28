# Individual-level BHHH Varmat for binary and polytomous LCA

Computes the outer-product (BHHH) information matrix and
variance-covariance matrix for LCA measurement model parameters in the
unconstrained (logit/log-ratio) space, matching multilevLCA's `$Varmat`.

## Usage

``` r
lca_indiv_varmat(Y.exp, mDesign.exp, fit0, ivItemcat, use.freq = TRUE)
```

## Arguments

- Y.exp:

  Expanded indicator matrix (N x sum(K_h)).

- mDesign.exp:

  Expanded design matrix (same dimensions as `Y.exp`), or `NULL` for
  complete data.

- fit0:

  Step-1 measurement model.

- ivItemcat:

  Number of categories for each item.

- use.freq:

  Logical. Collapse duplicate score vectors before computing the BHHH
  information matrix.

## Value

A list containing `Infomat`, `Varmat`, `SEs`, and the individual score
matrix `mScore`.

## Details

The score in unconstrained space is \\s\_{it} = u\_{it}(y_i - d_i \circ
p\_{it})\\, where \\d_i\\ is the missing-data design indicator matrix.
