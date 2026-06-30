# tseLCA: Three-Step Estimation for Latent Class Analysis

tseLCA implements bias-adjusted three-step estimators for structural
latent class models with covariates and distal outcomes. Building on the
efficient measurement-model estimation in multilevLCA (Lyrvall et al.,
2025), tseLCA adds modern three-step estimators, classification-error
corrections, and analytic sandwich variance estimation that propagates
measurement uncertainty from the latent class measurement stage through
to the final structural parameter estimates.

In contrast to one-step approaches such as poLCA, where including
covariates may alter the underlying latent class definitions, three-step
methods fix the measurement model before estimating structural
relationships and adjust for classification error at the final stage.
tseLCA also allows the measurement and structural models to be estimated
on different datasets, enabling researchers to calibrate a measurement
model on a large reference sample and apply it to a separate analysis
sample.

## The three-step approach

1.  **Measurement model**: estimate class-conditional item-response
    probabilities \\\phi\\ and class prevalences \\\pi\\ using
    multilevLCA (Lyrvall et al., 2025 ).

2.  **Classification-error matrix**: assign posterior class
    probabilities and compute the T x T misclassification matrix \\P(W =
    s \mid X = t)\\, with standard errors corrected for
    classification-error propagation (Bakk, Oberski & Vermunt, 2014).

3.  **Structural model**: estimate covariate effects using two-step
    starting values (Bakk & Kuha, 2018) and/or distal outcome (Bakk,
    Tekle & Vermunt) means with either the ML correction (Vermunt, 2010)
    or the BCH correction (Bolck, Croon & Hagenaars, 2004).

## Main functions

- [`three_step`](https://samleebyu.github.io/tseLCA/reference/three_step.md):

  Full three-step estimation pipeline. Accepts covariates (`Zp.names`),
  distal outcomes (`Zo.name`), or both. Handles Steps 1–3 in a single
  call, with optional pre-fitted Step-1 input through `step1`.

- [`lca_step1`](https://samleebyu.github.io/tseLCA/reference/lca_step1.md):

  Standalone Step-1 measurement model estimation via multilevLCA.
  Returns a reusable fit object that can be passed to
  [`three_step()`](https://samleebyu.github.io/tseLCA/reference/three_step.md)
  to avoid re-estimating the measurement model across multiple
  structural specifications.

- [`fitZ_from_fit0`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_fit0.md):

  Two-step covariate estimation by fixing measurement parameters at
  their Step-1 values and estimating multinomial logit coefficients
  \\\gamma\\ with an EM algorithm. Returns starting values for Step 3.
  Custom starting values can be supplied via `starting_val`.

- [`fitZ_from_multiLCA`](https://samleebyu.github.io/tseLCA/reference/fitZ_from_multiLCA.md):

  Two-step covariate estimation via `multiLCA(fixedpars = 1)`, returning
  multilevLCA's bias-corrected standard errors. Called automatically
  when `get.twostep.vcov = TRUE` in
  [`three_step`](https://samleebyu.github.io/tseLCA/reference/three_step.md).

- [`generate_data`](https://samleebyu.github.io/tseLCA/reference/generate_data.md):

  Simulate data replicating the Bakk & Kuha (2018) three-class design:
  six binary indicators across three separation levels (`"low"`,
  `"medium"`, `"high"`) and two scenarios (`"covariate"` or `"distal"`).

## Estimators

- ML (default, `use.bch = FALSE`):

  The Vermunt (2010) ML correction uses a weighted pseudo-likelihood
  with the misclassification matrix as a bias adjustment. Preferred when
  class separation is low or moderate.

- BCH (`use.bch = TRUE`):

  The Bolck, Croon & Hagenaars (2004) correction inverts the
  misclassification matrix to obtain direct class weights. Works well
  under high separation but may produce an ill-conditioned Hessian
  (non-positive semi-definite covariance matrix) when separation is low;
  use the ML estimator in that case.

## Variance estimation

- Full correction (`use.simple.cov = FALSE`, default):

  Analytic propagation of Step-1 measurement uncertainty through the
  classification-error correction, following Bakk, Oberski & Vermunt
  (2014). Uses soft (proportional) posteriors for the Jacobian
  \\\partial\theta_2/\partial\theta_1\\ regardless of
  `use.modal.assignment`. Recommended when separation is moderate or
  low.

- Simple/robust (`use.simple.cov = TRUE`):

  Sandwich SEs from Step 3 only, ignoring measurement uncertainty. A
  useful computational shortcut when separation is high and the
  correction is negligible.

## Class assignment

- Modal (`use.modal.assignment = TRUE`, default):

  Each observation is assigned to its most probable class (hard
  assignment). The Jacobian for variance correction is still computed
  from soft posteriors.

- Proportional (`use.modal.assignment = FALSE`):

  Soft posterior weights are used throughout Steps 2 and 3. Recommended
  when separation is moderate or low, and required for a mathematically
  well-defined analytic Jacobian.

## Supported features

- Binary and polytomous indicators, following multilevLCA coding
  conventions.

- Gaussian, Poisson, and binomial distal outcome families.

- Full-information maximum likelihood (FIML) for partially observed
  indicator patterns (`incomplete = TRUE`). Step 3 always performs
  listwise deletion on missing covariates or distal outcomes.

- Flexible measurement and structural samples: fit the measurement model
  on a reference sample and apply it to a different analysis sample via
  the `step1` argument.

- Arbitrary reference class for the multinomial logit parameterization
  via the `rebase` argument. Log-likelihoods are invariant to this
  choice.

- Joint covariate and distal outcome estimation (`Zp.names` and
  `Zo.name` supplied together). The covariate model is estimated first;
  covariate-adjusted posteriors are then used as priors in the distal
  outcome step.

- S3 methods (`print`, `summary`, `coef`, `vcov`, `plot`) for all four
  return subclasses: `tseLCA_measurement`, `tseLCA_covariate`,
  `tseLCA_distal`, `tseLCA_both`.

## Getting started

    # Introductory vignette
    vignette("tseLCA-workflow", package = "tseLCA")

## References

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, 22(4), 520–540.
<https://www.jstor.org/stable/24573086>

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, 12(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025).
Multilevel latent class analysis: State-of-the-art methodologies and
their implementation in the R package multilevLCA. *Multivariate
Behavioral Research*, 60(4), 731–747.
[doi:10.1080/00273171.2025.2473935](https://doi.org/10.1080/00273171.2025.2473935)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

## See also

Useful links:

- <https://samleebyu.github.io/tseLCA/>

- <https://github.com/SamLeeBYU/tseLCA>

- Report bugs at <https://github.com/SamLeeBYU/tseLCA/issues>

## Author

Sam Lee <samlee@arizona.edu>
