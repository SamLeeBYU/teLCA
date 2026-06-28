# tseLCA/R/lca_measurement.R
#
# Step-1 measurement model via multilevLCA::multiLCA.
#
# Exports:
#   lca_step1()        - measurement model fit + optional two-step fitZ
#   fitZ_from_fit0()   - pure-R EM for gamma with mPhi fixed (default fitZ path)
#   fitZ_from_multiLCA() - two-step via multiLCA(fixedpars=1, Z=...) (used when
#                          get.twostep.vcov = TRUE in three_step())

# -- lca_step1 -----------------------------------------------------------------

#' Individual-level BHHH Varmat for binary and polytomous LCA
#'
#' Computes the outer-product (BHHH) information matrix and variance-covariance
#' matrix for LCA measurement model parameters in the unconstrained
#' (logit/log-ratio) space, matching \pkg{multilevLCA}'s \code{$Varmat}.
#'
#' The score in unconstrained space is
#' \eqn{s_{it} = u_{it}(y_i - d_i \circ p_{it})},
#' where \eqn{d_i} is the missing-data design indicator matrix.
#'
#' @param Y.exp Expanded indicator matrix (N x sum(K_h)).
#' @param mDesign.exp Expanded design matrix (same dimensions as \code{Y.exp}),
#'   or \code{NULL} for complete data.
#' @param fit0 Step-1 measurement model.
#' @param ivItemcat Number of categories for each item.
#' @param use.freq Logical. Collapse duplicate score vectors before computing
#'   the BHHH information matrix.
#'
#' @return A list containing \code{Infomat}, \code{Varmat},
#'   \code{SEs}, and the individual score matrix \code{mScore}.
#' Individual-level BHHH Varmat for binary and polytomous LCA
#'
#' Computes the outer-product (BHHH) information matrix and variance-covariance
#' matrix for LCA measurement model parameters in the unconstrained
#' (logit/log-ratio) space.
#'
#' Assumes \code{fit0$mPhi} has the following structure (from multilevLCA):
#' \itemize{
#'   \item Dichotomous item h (\code{ivItemcat[h] == 2}): 1 row =
#'     \eqn{P(Y=1|C)}. The base level \eqn{P(Y=0|C)} is excluded.
#'   \item Polytomous item h (\code{ivItemcat[h] > 2}): \code{K_h} rows =
#'     \eqn{P(Y=0|C), \ldots, P(Y=K_h-1|C)}. The base level IS included.
#' }
#' \code{expand_Y} produces one-hot columns in the same order so that
#' \code{expand_Phi(fit0$mPhi, ivItemcat)} aligns column-wise with
#' \code{expand_Y(mY, ivItemcat)}.
#'
#' Free (estimable) parameters per item:
#' \itemize{
#'   \item Dichotomous: the single row of \code{mPhi} (\eqn{P(Y=1|C)}).
#'   \item Polytomous: rows 2..K_h of \code{mPhi} (\eqn{P(Y=1|C), \ldots}).
#'     Row 1 (\eqn{P(Y=0|C)}) is the reference and is not a free parameter.
#' }
#'
#' Boundary parameters (within \code{boundary.tol} of 0 or 1) are treated as
#' fixed: their score columns are zeroed so they do not contribute to the
#' information matrix.
#'
#' @param Y.exp       Expanded one-hot matrix (N x sum(K_h)).
#' @param mDesign.exp Expanded design matrix (same dims), or \code{NULL}.
#' @param fit0        Step-1 fit with \code{$vPi} and \code{$mPhi}.
#' @param ivItemcat   Integer vector of category counts per item.
#' @param boundary.tol Scalar tolerance for boundary detection. Default
#'   \code{1e-2}.
#' @param use.freq    Collapse duplicate score rows before cross-product.
#'   Default \code{TRUE}.
#'
#' @return List with \code{$Infomat}, \code{$Varmat}, \code{$SEs},
#'   \code{$mScore}.
lca_indiv_varmat <- function(
  Y.exp,
  mDesign.exp,
  fit0,
  ivItemcat,
  boundary.tol = 1e-2,
  use.freq = TRUE
) {
  pi_ <- fit0$vPi
  phi <- fit0$mPhi
  T <- length(pi_)
  N <- nrow(Y.exp)

  if (is.null(mDesign.exp)) {
    mDesign.exp <- matrix(1L, N, ncol(Y.exp))
  }

  # ---- Boundary flags --------------------------------------------------------
  # pi_bdry: length T. Reference class (t=1) is not a free parameter.
  pi_bdry <- pi_ <= boundary.tol | pi_ >= (1 - boundary.tol)

  # phi_bdry: n_mPhi_rows x T, where n_mPhi_rows = sum(ifelse(K==2, 1, K)).
  # Rows align with mPhi row order.
  phi_bdry <- phi <= boundary.tol | phi >= (1 - boundary.tol)

  pi_[pi_bdry] <- pmax(pmin(pi_[pi_bdry], 1 - 1e-6), 1e-6)
  phi[phi_bdry] <- pmax(pmin(phi[phi_bdry], 1 - 1e-6), 1e-6)

  # ---- Build theta1 and compute posteriors via compute_posteriors ------------
  # theta1 = c(pi_2..pi_T, phi_free) -- the same flattened free-parameter
  # vector used throughout lca_step2. We extract phi_free from mPhi using
  # the same free_idx logic as lca_step2.
  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )
  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))
  phi_free <- phi[free_idx, , drop = FALSE]
  theta1 <- c(pi_[-1L], phi_free)

  u_post <- compute_posteriors(Y.exp, mDesign.exp, theta1, ivItemcat, T)

  # ---- Expand phi for residual computation -----------------------------------
  phi_exp <- expand_Phi(phi, ivItemcat) # K_total x T

  # ---- Build free_cols (columns of phi_exp for free parameters) -------------
  # free_cols: indices into phi_exp rows for the free categories.
  # phi_exp has K_total rows; free categories are the non-reference rows
  # within each item block. The mapping from free_idx (mPhi rows) to
  # phi_exp rows differs between dichotomous and polytomous items:
  #
  #   Dichotomous (K=2): phi_exp block = [P(Y=0), P(Y=1)], 2 rows.
  #     free_idx points to the single mPhi row = P(Y=1) = phi_exp row 2
  #     within the block (col_start + 1).
  #
  #   Polytomous (K>2): phi_exp block = [P(Y=0)..P(Y=K-1)], K rows.
  #     free_idx points to mPhi rows 2..K within the block (P(Y=1)..P(Y=K-1))
  #     = phi_exp rows col_start+1 .. col_start+K-1.
  #
  # In both cases the free phi_exp rows are exactly col_start+1..col_start+K-1
  # (dropping col_start = reference P(Y=0) for binary and polytomous alike).

  free_cols <- integer(0L)
  col_start <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    free_cols <- c(free_cols, (col_start + 1L):(col_start + K_h - 1L))
    col_start <- col_start + K_h
  }
  n_free_phi <- length(free_cols) # = sum(ivItemcat - 1) = nrow(phi_free) for all items

  # phi_bdry_free: n_free_phi x T, boundary flags for free phi parameters.
  # free_idx already selects the free mPhi rows in item order, so:
  phi_bdry_free <- phi_bdry[free_idx, , drop = FALSE]

  # ---- Pi scores (T-1 columns, one per non-reference class t=2..T) -----------
  s_u_pi <- sweep(u_post[, -1L, drop = FALSE], 2L, pi_[-1L], "-")
  pi_free_bdry <- pi_bdry[-1L] # drop reference class t=1
  if (any(pi_free_bdry)) {
    s_u_pi[, pi_free_bdry] <- 0
  }

  # ---- Phi scores (n_free_phi * T columns, class-major) ----------------------
  Y_free <- Y.exp[, free_cols, drop = FALSE]
  D_free <- mDesign.exp[, free_cols, drop = FALSE]

  s_u_phi <- matrix(0, N, n_free_phi * T)

  for (t in seq_len(T)) {
    idx <- ((t - 1L) * n_free_phi + 1L):(t * n_free_phi)
    resid <- Y_free -
      D_free *
        matrix(
          phi_exp[free_cols, t],
          nrow = N,
          ncol = n_free_phi,
          byrow = TRUE
        )
    s_col <- u_post[, t] * resid

    # Zero boundary free parameters for this class
    bdry_t <- phi_bdry_free[, t]
    if (any(bdry_t)) {
      s_col[, bdry_t] <- 0
    }

    s_u_phi[, idx] <- s_col
  }

  S <- cbind(s_u_pi, s_u_phi) # N x p

  # ---- Identify active (non-boundary) columns --------------------------------
  # Boundary parameters have all-zero score columns -- removing them before
  # qr.solve avoids rank deficiency, then we restore zero rows/cols after.
  active <- which(colSums(S != 0) > 0L)
  p_full <- ncol(S)

  # ---- BHHH information matrix (on active columns only) ----------------------
  S_active <- S[, active, drop = FALSE]

  if (use.freq) {
    S_char <- apply(S_active, 1L, paste, collapse = "\r")
    uniq <- !duplicated(S_char)
    freq <- tabulate(match(S_char, S_char[uniq]))
    Infomat_active <- crossprod(S_active[uniq, , drop = FALSE] * sqrt(freq)) / N
  } else {
    Infomat_active <- crossprod(S_active) / N
  }

  Varmat_active <- tryCatch(
    qr.solve(Infomat_active) / N,
    error = function(e) {
      warning(
        "lca_indiv_varmat: Infomat is singular even after removing boundary ",
        "parameters; returning NA matrix. Check for near-empty classes.",
        call. = FALSE
      )
      matrix(NA_real_, length(active), length(active))
    }
  )

  # ---- Restore full-size Infomat and Varmat ----------------------------------
  # Boundary parameters get zero rows/cols in Infomat (no information)
  # and zero rows/cols in Varmat (variance treated as zero / fixed).
  Infomat <- matrix(0, p_full, p_full)
  Infomat[active, active] <- Infomat_active

  Varmat <- matrix(0, p_full, p_full)
  Varmat[active, active] <- Varmat_active

  list(
    Infomat = Infomat,
    Varmat = Varmat,
    SEs = sqrt(diag(Varmat)),
    mScore = S
  )
}

#' Fit the LCA measurement model (Step 1)
#'
#' Estimates the latent class measurement model via \pkg{multilevLCA} and,
#' optionally, fixes `mPhi` and estimates covariate effects (two-step
#' initialization) via `fitZ_from_fit0()`.
#'
#' @param data A data.frame containing at minimum the indicator columns.
#' @param Y.names Character vector of item column names.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or `NULL`.
#' @param maxIter.measurement Maximum EM iterations. Default `5000L`.
#' @param measurement.tol Convergence tolerance. Default `1e-8`.
#' @param covariate.tol Convergence tolerance for the `fitZ` BFGS M-step.
#' @param iter.measurement Number of random restarts when entropy R\eqn{^2} is low.
#' @param R2.threshold Entropy R\eqn{^2} below which restarts are triggered.
#' @param use.two.step Logical. If `TRUE`, also estimate `fitZ` via
#'   `fitZ_from_fit0()`.
#' @param estimate.one.step Logical. If `FALSE`, skip the unconditional EM and
#'   only compute `fitZ`.
#' @param incomplete Logical. FIML for partially missing indicators.
#' @param maxIter.fitZ Maximum BFGS-EM iterations for `fitZ_from_fit0()`.
#' @param include.intercept Logical. Prepend intercept to covariate design matrix.
#' @param rebase Character or integer specifying the reference latent class.
#'   Use `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
#'   measurement model is permuted so this class becomes column 1, making it
#'   the reference for all downstream multinomial logit parameterizations.
#' @param verbose Logical. Print progress messages. Default `FALSE`.
#'
#' @return A list with `$fit0` (multilevLCA measurement model) and `$fitZ`
#'   (two-step covariate model from `fitZ_from_fit0`, or `NULL`).
#' @examples
#' \donttest{
#' d <- generate_data(200, "high", "covariate", seed = 1)
#'
#' # Measurement model only
#' s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' s1$fit0$vPi    # estimated class prevalences
#' s1$fit0$mPhi   # item-response probabilities
#'
#' # With two-step covariate initialization
#' s1z <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                  Zp.names = "Zp", use.two.step = TRUE, verbose = TRUE)
#' s1z$fitZ$mGamma   # two-step gamma estimates
#' }
#' @export
lca_step1 <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.two.step = TRUE,
  estimate.one.step = TRUE,
  incomplete = FALSE,
  maxIter.fitZ = 200L,
  include.intercept = TRUE,
  rebase = "C1",
  verbose = FALSE
) {
  run_measurement_fit <- function(extra_args = list()) {
    args <- c(
      list(
        data,
        Y.names,
        n_classes,
        extout = TRUE,
        incomplete = incomplete,
        maxIter = maxIter.measurement,
        tol = measurement.tol,
        verbose = FALSE
      ),
      extra_args
    )
    fit <- do.call(multilevLCA::multiLCA, args)
    if (nrow(fit$LLKSeries) == maxIter.measurement) {
      args$maxIter <- 2L * maxIter.measurement
      fit <- do.call(multilevLCA::multiLCA, args)
      if (verbose) {
        warning(sprintf(
          "Measurement model hit %d iterations; retried with %d. Low separation is likely the cause.",
          maxIter.measurement,
          2L * maxIter.measurement
        ))
      }
      if (nrow(fit$LLKSeries) == 2L * maxIter.measurement) {
        warning(
          "Measurement model still failed to converge even after running more iterations. Consider increasing maxIter.measurement and or measurement.tol"
        )
      }
    }
    fit
  }

  best_fit <- function(initial, run_fn) {
    ll0 <- initial$LLKSeries[nrow(initial$LLKSeries), 1L]
    if (is.null(initial$R2entr) || initial$R2entr >= R2.threshold) {
      return(initial)
    }
    if (verbose) {
      warning(sprintf(
        "Measurement model has low entropy R\u00b2 (%.3f < %.3f). Running %d additional random restarts.",
        initial$R2entr,
        R2.threshold,
        iter.measurement
      ))
    }
    if (iter.measurement > 0L) {
      cands <- lapply(seq_len(iter.measurement), function(r) run_fn())
      cand_lls <- vapply(
        cands,
        function(f) f$LLKSeries[nrow(f$LLKSeries), 1L],
        numeric(1L)
      )
      best_r <- which.max(cand_lls)
      if (cand_lls[best_r] > ll0) {
        if (verbose) {
          message(sprintf(
            "Restart %d improved log-likelihood to %.4f.",
            best_r,
            cand_lls[best_r]
          ))
        }
        cands[[best_r]]
      } else {
        if (verbose) {
          message("No restart improved on the initial measurement model.")
        }
        initial
      }
    } else {
      initial
    }
  }

  fit0 <- if (estimate.one.step) {
    best_fit(initial = run_measurement_fit(), run_fn = run_measurement_fit)
  } else {
    NULL
  }

  # Permute classes so the desired reference is column 1
  if (!is.null(fit0)) {
    ref_idx <- parse_rebase(rebase, n_classes)
    fit0 <- permute_fit0_classes(fit0, ref_idx)
  }

  fitZ <- if (use.two.step && !is.null(Zp.names) && !is.null(fit0)) {
    fitZ_from_fit0(
      fit0 = fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = covariate.tol,
      maxIter = maxIter.fitZ,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  } else {
    NULL
  }

  list(fit0 = fit0, fitZ = fitZ)
}


# -- fitZ_from_fit0 ------------------------------------------------------------

#' Estimate covariate effects with measurement parameters fixed (two-step EM)
#'
#' Fixes `mPhi` at `fit0$mPhi` and estimates multinomial logit coefficients
#' `mGamma` (Q x (T-1)) via an EM algorithm with a BFGS M-step.
#'
#' @param fit0 Output of `lca_step1()$fit0`.
#' @param data A data.frame.
#' @param Y.names Character vector of item column names.
#' @param Zp.names Character vector of covariate column names.
#' @param tol Convergence tolerance.
#' @param maxIter Maximum EM iterations.
#' @param incomplete Logical.
#' @param include.intercept Logical.
#' @param rebase Character or integer. Reference class for the multinomial logit
#'   parameterization (e.g. `"C1"`, `"C2"`, or an integer). Default `"C1"`.
#'   Must match the `rebase` used in `lca_step1()` so class column ordering
#'   is consistent.
#' @param starting_val Optional Q x (T-1) starting value matrix for `mGamma`.
#' @param verbose Logical. Print convergence messages. Default `FALSE`.
#'
#' @return A list with `$mGamma` (Q x (T-1)), `$mPhi`, `$vOmega`,
#'   `$LLKSeries`, `$converged`, `$n_obs`.
#' @examples
#' \donttest{
#' d  <- generate_data(200, "high", "covariate", seed = 1)
#' s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#'
#' # Estimate two-step gamma with mPhi fixed at Step-1 values
#' fZ <- fitZ_from_fit0(
#'   fit0     = s1$fit0,
#'   data     = d,
#'   Y.names  = paste0("Y", 1:6),
#'   Zp.names = "Zp",
#'   verbose  = TRUE
#' )
#' fZ$mGamma   # Q x (T-1) coefficient matrix
#' fZ$converged
#' }
#' @export
fitZ_from_fit0 <- function(
  fit0,
  data,
  Y.names,
  Zp.names,
  tol = 1e-6,
  maxIter = 200L,
  incomplete = FALSE,
  include.intercept = TRUE,
  rebase = "C1",
  starting_val = NULL,
  verbose = FALSE
) {
  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    incomplete = incomplete,
    include.intercept = include.intercept,
    verbose = verbose
  )
  mY <- cd$Y.obs # expanded N_Y x K
  mDesign <- cd$mDesign
  ivItemcat <- cd$ivItemcat
  # For fitZ we need the Z rows that overlap with the Y-kept rows
  mZ <- cd$Z_mat # N_Z x Q, already complete-case
  # Restrict Y to the Z-complete rows
  mY <- mY[cd$keep_step3_Z_in_Y, , drop = FALSE]
  if (!is.null(mDesign)) {
    mDesign <- mDesign[cd$keep_step3_Z_in_Y, , drop = FALSE]
  }

  mPhi <- expand_Phi(fit0$mPhi, ivItemcat)
  iT <- ncol(mPhi)
  iN <- nrow(mY)
  iP <- ncol(mZ)

  phi_clamped <- pmax(pmin(mPhi, 1 - 1e-10), 1e-10)
  log_p_it <- if (is.null(mDesign)) {
    mY %*% log(phi_clamped)
  } else {
    (mDesign * mY) %*% log(phi_clamped)
  }

  softmax_rows <- function(mat) {
    mat <- mat - apply(mat, 1L, max)
    ex <- exp(mat)
    ex / rowSums(ex)
  }

  gamma <- matrix(0, nrow = iP, ncol = iT - 1L)
  if (!is.null(starting_val)) {
    if (!isTRUE(all.equal(dim(gamma), dim(starting_val)))) {
      warning(sprintf("starting_val dimensions must be %d x %d.", iP, iT - 1L))
    } else {
      gamma <- starting_val
    }
  }

  ll_prev <- -Inf
  LLKSeries <- numeric(0L)

  for (iter in seq_len(maxIter)) {
    eta_full <- cbind(0, mZ %*% gamma)
    pi_mat <- softmax_rows(eta_full)
    log_joint <- log_p_it + log(pi_mat)
    log_marg <- apply(log_joint, 1L, function(row) {
      mx <- max(row)
      mx + log(sum(exp(row - mx)))
    })
    ll_curr <- sum(log_marg)
    LLKSeries <- c(LLKSeries, ll_curr)
    w_mat <- exp(log_joint - log_marg)

    gamma_new <- tryCatch(
      {
        obj <- function(g_vec) {
          g_mat <- matrix(g_vec, nrow = iP, ncol = iT - 1L)
          pi_ <- softmax_rows(cbind(0, mZ %*% g_mat))
          -sum(w_mat * log(pi_))
        }
        gr <- function(g_vec) {
          g_mat <- matrix(g_vec, nrow = iP, ncol = iT - 1L)
          pi_ <- softmax_rows(cbind(0, mZ %*% g_mat))
          resid <- w_mat[, -1L, drop = FALSE] - pi_[, -1L, drop = FALSE]
          -as.vector(t(mZ) %*% resid)
        }
        res <- optim(par = as.vector(gamma), fn = obj, gr = gr, method = "BFGS")
        matrix(res$par, nrow = iP, ncol = iT - 1L)
      },
      error = function(e) {
        warning(
          "fitZ_from_fit0: optim failed at iter ",
          iter,
          ": ",
          conditionMessage(e),
          ". Keeping previous gamma."
        )
        gamma
      }
    )

    if (iter > 1L && abs(ll_curr - ll_prev) < tol) {
      gamma <- gamma_new
      if (verbose) {
        message(sprintf("fitZ EM converged in %d iterations.", iter))
      }
      break
    }
    gamma <- gamma_new
    ll_prev <- ll_curr
  }

  converged <- (length(LLKSeries) < maxIter) ||
    (abs(LLKSeries[length(LLKSeries)] - LLKSeries[length(LLKSeries) - 1L]) <
      tol)
  if (!converged) {
    warning(
      "fitZ_from_fit0: gamma EM did not converge in ",
      maxIter,
      " iterations."
    )
  }

  pi_final <- softmax_rows(cbind(0, mZ %*% gamma))
  vOmega <- colMeans(pi_final)

  rownames(gamma) <- colnames(mZ)
  # Column names reflect the non-reference classes
  # (all classes except the reference, in ascending order)
  ref_idx <- parse_rebase(rebase, iT)
  non_ref_classes <- seq_len(iT)[-ref_idx]
  colnames(gamma) <- paste0("C", non_ref_classes)

  list(
    mGamma = gamma,
    mPhi = mPhi,
    vOmega = vOmega,
    LLKSeries = matrix(LLKSeries, ncol = 1L),
    converged = converged,
    n_obs = iN
  )
}


# -- fitZ_from_multiLCA --------------------------------------------------------

#' Estimate two-step covariate model via multilevLCA (optional reference path)
#'
#' Calls `multilevLCA::multiLCA` with `fixedpars = 1` and `Z = Zp.names` to
#' fit the two-step covariate model.  This is the original multilevLCA approach
#' and is used when `get.twostep.vcov = TRUE` in [tseLCA::three_step()] to obtain
#' multilevLCA's corrected standard errors for the two-step gamma estimates.
#'
#' @param data A data.frame.
#' @param Y.names Character vector of item column names.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names.
#' @param maxIter.measurement Maximum EM iterations.
#' @param measurement.tol Convergence tolerance.
#' @param covariate.tol NR tolerance for the covariate model.
#' @param iter.measurement Number of random restarts.
#' @param R2.threshold Entropy R\eqn{^2} restart threshold.
#' @param incomplete Logical.
#' @param rebase Character or integer. Reference class for column naming of
#'   `$mGamma`. Must match the `rebase` used in [tseLCA::three_step()] so
#'   coefficient labels are consistent. Default `"C1"`.
#' @param verbose Logical.
#'
#' @return A list with `$mGamma`, `$mPhi`, `$vOmega`, `$LLKSeries`, and
#'   `$raw_fit` (the full multilevLCA output, including `$Varmat_cor` and
#'   `$SEs_cor_gamma` if available).
#' @examples
#' \donttest{
#' d <- generate_data(200, "high", "covariate", seed = 1)
#'
#' # Two-step estimation via multiLCA (fixedpars = 1)
#' fZ_ml <- fitZ_from_multiLCA(
#'   data                = d,
#'   Y.names             = paste0("Y", 1:6),
#'   n_classes           = 3,
#'   Zp.names            = "Zp",
#'   maxIter.measurement = 5000L,
#'   measurement.tol     = 1e-8,
#'   covariate.tol       = 1e-6,
#'   iter.measurement    = 10L,
#'   R2.threshold        = 0.70
#' )
#' fZ_ml$mGamma           # two-step estimates
#' fZ_ml$raw_fit$Varmat_cor   # multilevLCA corrected vcov
#' }
#' @export
fitZ_from_multiLCA <- function(
  data,
  Y.names,
  n_classes,
  Zp.names,
  maxIter.measurement,
  measurement.tol,
  covariate.tol,
  iter.measurement,
  R2.threshold,
  incomplete = FALSE,
  rebase = "C1",
  verbose = FALSE
) {
  run_fit <- function() {
    args <- list(
      data,
      Y.names,
      n_classes,
      Z = Zp.names,
      extout = TRUE,
      incomplete = incomplete,
      maxIter = maxIter.measurement,
      tol = measurement.tol,
      NRtol = covariate.tol,
      fixedpars = 1L,
      verbose = FALSE
    )
    fit <- do.call(multilevLCA::multiLCA, args)
    if (nrow(fit$LLKSeries) == maxIter.measurement) {
      args$maxIter <- 2L * maxIter.measurement
      fit <- do.call(multilevLCA::multiLCA, args)
      if (verbose) {
        warning(sprintf(
          "fitZ multiLCA hit %d iterations; retried with %d.",
          maxIter.measurement,
          2L * maxIter.measurement
        ))
      }
    }
    fit
  }

  initial <- run_fit()
  ll0 <- initial$LLKSeries[nrow(initial$LLKSeries), 1L]

  if (!is.null(initial$R2entr) && initial$R2entr < R2.threshold) {
    if (verbose) {
      warning(sprintf(
        "fitZ multiLCA has low entropy R\u00b2 (%.3f < %.3f). Running %d additional random restarts.",
        initial$R2entr,
        R2.threshold,
        iter.measurement
      ))
    }
    if (iter.measurement > 0L) {
      cands <- lapply(seq_len(iter.measurement), function(r) run_fit())
      cand_lls <- vapply(
        cands,
        function(f) f$LLKSeries[nrow(f$LLKSeries), 1L],
        numeric(1L)
      )
      best_r <- which.max(cand_lls)
      if (cand_lls[best_r] > ll0) {
        if (verbose) {
          message(sprintf(
            "fitZ restart %d improved log-likelihood to %.4f.",
            best_r,
            cand_lls[best_r]
          ))
        }
        initial <- cands[[best_r]]
      } else {
        if (verbose) {
          message(
            "No fitZ restart improved on the initial multiLCA covariate fit."
          )
        }
      }
    }
  }

  raw <- initial
  mGamma <- raw$mGamma
  rownames(mGamma) <- c("Intercept", Zp.names)
  ref_idx <- parse_rebase(rebase, n_classes)
  non_ref_classes <- seq_len(n_classes)[-ref_idx]
  colnames(mGamma) <- paste0("C", non_ref_classes)

  list(
    mGamma = mGamma,
    mPhi = raw$mPhi,
    vOmega = as.vector(raw$vPi_avg),
    LLKSeries = raw$LLKSeries,
    raw_fit = raw
  )
}
