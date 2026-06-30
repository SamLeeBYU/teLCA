# Permute class columns of a fit0 object so that class ref_idx is first

Reorders columns of mPhi and vPi so that the desired reference class
becomes column 1 before estimation. This ensures the multinomial logit
is parameterized with the correct baseline from the start.

## Usage

``` r
permute_fit0_classes(fit0, ref_idx)
```

## Arguments

- fit0:

  Raw multilevLCA fit object (has \$mPhi and \$vPi).

- ref_idx:

  Integer. Class index to move to position 1.

## Value

fit0 with columns permuted.
