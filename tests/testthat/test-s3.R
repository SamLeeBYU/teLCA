# tests/testthat/test-s3.R
# S3 dispatch tests using minimal hand-crafted tseLCA objects.
# No multilevLCA calls needed — tests dispatch, not estimation.

# ── Minimal mock objects ───────────────────────────────────────────────────────

.mock_vcov <- function(n) {
  V <- diag(runif(n, 0.01, 0.1))
  rownames(V) <- colnames(V) <- paste0("p", seq_len(n))
  V
}

.mock_covariate <- function() {
  coefs <- matrix(c(1.0, -0.5, 0.3, 0.8), nrow = 2L, ncol = 2L)
  rownames(coefs) <- c("Intercept", "Zp")
  colnames(coefs) <- c("C2", "C3")
  vcov <- .mock_vcov(4L)
  rownames(vcov) <- colnames(vcov) <-
    c("Intercept:C2", "Zp:C2", "Intercept:C3", "Zp:C3")
  structure(
    list(
      three_step = coefs,
      three_step_vcov = vcov,
      two_step = coefs + 0.1,
      two_step_vcov = NULL,
      llik = -300.0,
      AIC = 620.0,
      BIC = 640.0,
      n_classes = 3L,
      measurement_model = list(
        fit0 = list(vPi = c(0.4, 0.3, 0.3), mPhi = matrix(0.7, 6L, 3L))
      )
    ),
    class = c("tseLCA_covariate", "tseLCA")
  )
}

.mock_distal <- function() {
  mu <- c(mu_C1 = -1.0, mu_C2 = 0.0, mu_C3 = 1.0)
  vcov <- .mock_vcov(3L)
  rownames(vcov) <- colnames(vcov) <- names(mu)
  structure(
    list(
      three_step = mu,
      three_step_vcov = vcov,
      family = "gaussian",
      n_classes = 3L,
      measurement_model = list(
        fit0 = list(vPi = c(0.4, 0.3, 0.3), mPhi = matrix(0.7, 6L, 3L))
      )
    ),
    class = c("tseLCA_distal", "tseLCA")
  )
}

.mock_measurement <- function() {
  structure(
    list(
      measurement_model = list(
        fit0 = list(
          vPi = c(0.5, 0.5),
          mPhi = matrix(c(0.8, 0.2), nrow = 2L, ncol = 2L),
          AIC = 400.0,
          BIC = 420.0,
          LLKSeries = matrix(-200, 1L, 1L),
          R2entr = 0.85
        )
      ),
      llik = -200.0,
      AIC = 400.0,
      BIC = 420.0,
      R2entr = 0.85,
      n_classes = 2L
    ),
    class = c("tseLCA_measurement", "tseLCA")
  )
}

.mock_both <- function() {
  cov <- .mock_covariate()
  dis <- .mock_distal()
  structure(
    list(
      covariate = cov,
      distal = list(
        three_step = dis$three_step,
        three_step_vcov = dis$three_step_vcov
      ),
      family = "gaussian",
      n_classes = 3L
    ),
    class = c("tseLCA_both", "tseLCA")
  )
}

# ── class tests ───────────────────────────────────────────────────────────────

test_that("mock objects have correct classes", {
  expect_s3_class(.mock_measurement(), "tseLCA_measurement")
  expect_s3_class(.mock_measurement(), "tseLCA")
  expect_s3_class(.mock_covariate(), "tseLCA_covariate")
  expect_s3_class(.mock_distal(), "tseLCA_distal")
  expect_s3_class(.mock_both(), "tseLCA_both")
})

# ── coef() ────────────────────────────────────────────────────────────────────

test_that("coef.tseLCA_covariate returns three_step by default", {
  obj <- .mock_covariate()
  expect_identical(coef(obj), obj$three_step)
})

test_that("coef.tseLCA_covariate returns two_step when requested", {
  obj <- .mock_covariate()
  expect_identical(coef(obj, which = "two_step"), obj$two_step)
})

test_that("coef.tseLCA_distal returns the named mu vector", {
  obj <- .mock_distal()
  co <- coef(obj)
  expect_length(co, 3L)
  expect_named(co, c("mu_C1", "mu_C2", "mu_C3"))
})

test_that("coef.tseLCA_both dispatches to covariate and distal correctly", {
  obj <- .mock_both()
  expect_identical(coef(obj, which = "covariate"), obj$covariate$three_step)
  expect_identical(coef(obj, which = "distal"), obj$distal$three_step)
  both <- coef(obj, which = "both")
  expect_named(both, c("covariate", "distal"))
})

test_that("coef.tseLCA_measurement returns prevalences and item probs", {
  obj <- .mock_measurement()
  co <- coef(obj)
  expect_named(co, c("prevalences", "item_probs"))
  expect_length(co$prevalences, 2L)
})

# ── vcov() ────────────────────────────────────────────────────────────────────

test_that("vcov.tseLCA_covariate returns three_step_vcov by default", {
  obj <- .mock_covariate()
  expect_identical(vcov(obj), obj$three_step_vcov)
})

test_that("vcov.tseLCA_covariate errors informatively for missing two_step vcov", {
  obj <- .mock_covariate() # two_step_vcov is NULL in mock
  expect_error(vcov(obj, which = "two_step"), regexp = "get.twostep.vcov")
})

test_that("vcov.tseLCA_distal returns three_step_vcov", {
  obj <- .mock_distal()
  expect_identical(vcov(obj), obj$three_step_vcov)
})

test_that("vcov.tseLCA_both dispatches correctly", {
  obj <- .mock_both()
  expect_identical(
    vcov(obj, which = "covariate"),
    obj$covariate$three_step_vcov
  )
  expect_identical(vcov(obj, which = "distal"), obj$distal$three_step_vcov)
  both <- vcov(obj, which = "both")
  expect_named(both, c("covariate", "distal"))
})

# ── print() and summary() smoke tests ─────────────────────────────────────────

test_that("print.tseLCA_measurement produces output", {
  expect_output(print(.mock_measurement()), regexp = "tseLCA")
})

test_that("print.tseLCA_covariate produces output with coefficient table", {
  expect_output(print(.mock_covariate()), regexp = "Estimate")
})

test_that("print.tseLCA_distal produces output", {
  expect_output(print(.mock_distal()), regexp = "Estimate")
})

test_that("print.tseLCA_both produces output for both branches", {
  out <- capture_output(print(.mock_both()))
  expect_match(out, "Covariate")
  expect_match(out, "Distal")
})

test_that("summary.tseLCA_covariate shows two_step estimates", {
  out <- capture_output(summary(.mock_covariate()))
  expect_match(out, "Two-step")
  expect_match(out, "Intercept")
})

test_that("p-value significance stars appear in coefficient table", {
  obj <- .mock_covariate()
  # Force a very small p-value by giving a large z-value via tiny SE
  obj$three_step_vcov <- diag(rep(1e-6, 4L))
  rownames(obj$three_step_vcov) <- colnames(obj$three_step_vcov) <-
    c("Intercept:C2", "Zp:C2", "Intercept:C3", "Zp:C3")
  out <- capture_output(print(obj))
  expect_match(out, "\\*")
})
