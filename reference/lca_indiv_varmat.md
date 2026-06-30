# Individual-level BHHH variance matrix for binary and polytomous LCA

Computes the outer-product (BHHH) information matrix and
variance-covariance matrix for LCA measurement model parameters in the
unconstrained (logit/log-ratio) space, matching multilevLCA's `$Varmat`.

## Usage

``` r
lca_indiv_varmat(
  Y.exp,
  mDesign.exp,
  fit0,
  ivItemcat,
  boundary.tol = 0.01,
  use.freq = TRUE,
  u_post = NULL
)
```

## Arguments

- Y.exp:

  N x sum(K_h) expanded one-hot indicator matrix.

- mDesign.exp:

  Expanded design matrix (same dimensions as `Y.exp`), or `NULL` for
  complete data.

- fit0:

  Step-1 fit object with `$vPi` and `$mPhi`.

- ivItemcat:

  Integer vector of category counts per item.

- boundary.tol:

  Scalar tolerance for boundary detection. Default `1e-2`.

- use.freq:

  Logical. Collapse duplicate score rows before computing the
  cross-product, weighting by frequency. Default `TRUE`.

- u_post:

  Optional N x T matrix of posterior class probabilities. When supplied
  (e.g. extracted from `fit0$mU` via `extract_Y_from_mU`),
  `compute_posteriors` is skipped. Default `NULL`.

## Value

A list with the following elements:

- `Infomat`:

  Square BHHH information matrix of dimension p x p, where p = (T-1) +
  sum(ivItemcat - 1) \* T is the total number of free parameters.
  Boundary parameters have zero rows and columns.

- `Varmat`:

  Inverse of `Infomat` divided by N, giving the asymptotic
  variance-covariance matrix on the same scale as multilevLCA's
  `$Varmat`. Boundary parameters have zero rows and columns.

- `SEs`:

  Numeric vector of length p. Square root of the diagonal of `Varmat`;
  zero for boundary parameters.

- `mScore`:

  N x p matrix of individual score contributions in the unconstrained
  parameterization, used for sandwich variance propagation in `lca_vcov`
  and `lca_vcov_distal`.

## Details

The score in unconstrained space is \\s\_{it} = u\_{it}(y_i - d_i \circ
p\_{it})\\, where \\d_i\\ is the missing-data design indicator matrix.

Assumes `fit0$mPhi` follows the multilevLCA storage convention:

- Dichotomous item h (`ivItemcat[h] == 2`): 1 row = \\P(Y=1\|C)\\; the
  base level \\P(Y=0\|C)\\ is excluded.

- Polytomous item h (`ivItemcat[h] > 2`): `K_h` rows = \\P(Y=0\|C),
  \ldots, P(Y=K_h-1\|C)\\; the base level is included.

`expand_Y` produces one-hot columns in the same order so that
`expand_Phi(fit0$mPhi, ivItemcat)` aligns column-wise with
`expand_Y(mY, ivItemcat)`. Free (estimable) parameters per item are the
single \\P(Y=1\|C)\\ row for dichotomous items, and rows 2 through
\\K_h\\ for polytomous items (row 1, \\P(Y=0\|C)\\, is the reference).
Boundary parameters (within `boundary.tol` of 0 or 1) are treated as
fixed: their score columns are zeroed and they do not contribute to the
information matrix.
