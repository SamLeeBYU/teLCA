# tseLCA/R/clean.R
#
# Shared data-cleaning utilities used by three_step() and fitZ_from_fit0().
#
# The central function is clean_data() which takes the raw data.frame and
# returns consistently prepared matrices for Y (expanded one-hot), Z
# (covariate design), Z0 (distal outcome), and mDesign (FIML mask).
#
# Row-filtering philosophy
# ------------------------
# Steps 1 & 2 (measurement model, posteriors):
#   Use as much Y data as possible.
#   - Drop rows where ALL Y items are missing (completely uninformative).
#   - If incomplete = FALSE, also drop rows where ANY Y is missing.
#   - Missingness in Z or Z0 is irrelevant at this stage.
#
# Step 3 covariate:
#   Restrict to rows that (a) passed the Y filter AND (b) have complete Z.
#
# Step 3 distal:
#   Restrict to rows that (a) passed the Y filter AND (b) have complete Z0.

#' Prepare and validate data for tseLCA estimation
#'
#' @keywords internal
#'
#' @param data       A data.frame.
#' @param Y.names    Character vector of item column names.
#' @param Zp.names   Character vector of covariate column names, or `NULL`.
#' @param Zo.name    Single distal outcome column name, or `NULL`.
#' @param incomplete Logical. If `TRUE`, use FIML for partially-observed Y.
#' @param include.intercept Logical. Prepend intercept column to Z.
#' @param verbose    Logical. Print row-drop messages.
#'
#' @return A named list with:
#' \describe{
#'   \item{Y.obs}{N_Y x K expanded one-hot indicator matrix for Steps 1 & 2.}
#'   \item{mDesign}{N_Y x K design/mask matrix (all 1s when incomplete = FALSE).}
#'   \item{ivItemcat}{Integer vector of category counts per item.}
#'   \item{keep_Y}{Integer indices of rows kept for Steps 1 & 2 (into original N).}
#'   \item{Z_mat}{N_Z x Q covariate design matrix, or NULL.}
#'   \item{keep_step3_Z_in_Y}{Positions of Z-complete rows within keep_Y.}
#'   \item{Z0_mat}{N_Z0 x 1 distal outcome matrix, or NULL.}
#'   \item{keep_step3_Z0_in_Y}{Positions of Z0-complete rows within keep_Y.}
#' }
clean_data <- function(
  data,
  Y.names,
  Zp.names = NULL,
  Zo.name = NULL,
  incomplete = FALSE,
  include.intercept = TRUE,
  verbose = FALSE
) {
  Y_raw <- as.matrix(data[, Y.names, drop = FALSE])
  N_full <- nrow(Y_raw)

  ivItemcat <- apply(Y_raw, 2L, \(x) length(na.omit(unique(x))))

  # ---- Y row filter -----------------------------------------------------------
  all_Y_missing <- rowSums(!is.na(Y_raw)) == 0L
  any_Y_missing <- rowSums(is.na(Y_raw)) > 0L

  drop_Y <- all_Y_missing
  if (!incomplete) {
    drop_Y <- drop_Y | any_Y_missing
  }

  keep_Y <- which(!drop_Y)

  if (any(drop_Y) && verbose) {
    message(sprintf(
      "%d row(s) dropped from measurement/classification steps (missing Y).",
      sum(drop_Y)
    ))
  }

  Y.obs <- Y_raw[keep_Y, , drop = FALSE]

  # ---- Expand Y and build mDesign ---------------------------------------------
  if (incomplete) {
    Y.obs_exp <- expand_Y(Y.obs, ivItemcat)
    mDesign <- (!is.na(Y.obs_exp)) * 1L
    Y.obs_exp[is.na(Y.obs_exp)] <- 0L
  } else {
    Y.obs[is.na(Y.obs)] <- 0L
    Y.obs_exp <- expand_Y(Y.obs, ivItemcat)
    mDesign <- NULL
  }

  # ---- Covariate design matrix ------------------------------------------------
  if (!is.null(Zp.names)) {
    Z_mat_full <- if (include.intercept) {
      m <- cbind(1, as.matrix(data[, Zp.names, drop = FALSE]))
      colnames(m) <- c("Intercept", Zp.names)
      m
    } else {
      as.matrix(data[, Zp.names, drop = FALSE])
    }
    any_Z_missing_full <- !complete.cases(Z_mat_full)

    drop_step3_Z <- drop_Y | any_Z_missing_full
    keep_step3_Z <- which(!drop_step3_Z)
    Z_mat <- Z_mat_full[keep_step3_Z, , drop = FALSE]
    keep_step3_Z_in_Y <- match(keep_step3_Z, keep_Y)

    if (any(is.na(keep_step3_Z_in_Y))) {
      stop("Internal error: Z rows not a subset of Y rows.", call. = FALSE)
    }

    if (sum(any_Z_missing_full[keep_Y]) > 0L && verbose) {
      message(sprintf(
        "%d row(s) excluded from covariate step (missing Z).",
        sum(any_Z_missing_full[keep_Y])
      ))
    }

    # Check for linear dependence
    Z_rank <- qr(Z_mat)$rank
    if (Z_rank < ncol(Z_mat)) {
      stop(
        sprintf(
          paste0(
            "Covariate design matrix is rank-deficient (rank %d, %d columns). ",
            "Check for perfectly collinear predictors or a redundant intercept."
          ),
          Z_rank,
          ncol(Z_mat)
        ),
        call. = FALSE
      )
    }
  } else {
    Z_mat <- NULL
    keep_step3_Z_in_Y <- seq_along(keep_Y)
  }

  # ---- Distal outcome matrix --------------------------------------------------
  if (!is.null(Zo.name)) {
    Z0_mat_full <- as.matrix(data[, Zo.name, drop = FALSE])
    any_Z0_missing_full <- !complete.cases(Z0_mat_full)

    drop_step3_Z0 <- drop_Y | any_Z0_missing_full
    keep_step3_Z0 <- which(!drop_step3_Z0)
    Z0_mat <- Z0_mat_full[keep_step3_Z0, , drop = FALSE]
    keep_step3_Z0_in_Y <- match(keep_step3_Z0, keep_Y)

    if (any(is.na(keep_step3_Z0_in_Y))) {
      stop("Internal error: Z0 rows not a subset of Y rows.", call. = FALSE)
    }

    if (sum(any_Z0_missing_full[keep_Y]) > 0L && verbose) {
      message(sprintf(
        "%d row(s) excluded from distal step (missing Z0).",
        sum(any_Z0_missing_full[keep_Y])
      ))
    }
  } else {
    Z0_mat <- NULL
    keep_step3_Z0_in_Y <- seq_along(keep_Y)
  }

  list(
    Y.obs = Y.obs_exp,
    mDesign = mDesign,
    ivItemcat = ivItemcat,
    keep_Y = keep_Y,
    Z_mat = Z_mat,
    keep_step3_Z_in_Y = keep_step3_Z_in_Y,
    Z0_mat = Z0_mat,
    keep_step3_Z0_in_Y = keep_step3_Z0_in_Y,
    keep_step3_Z0 = if (!is.null(Zo.name)) keep_step3_Z0 else integer(0L)
  )
}


#' Parse and validate the rebase argument
#'
#' @param rebase Character like "C2" or integer class index.
#' @param T      Total number of classes.
#' @return Integer class index (1-based) to use as reference.
#' @keywords internal
parse_rebase <- function(rebase, T) {
  if (is.character(rebase)) {
    if (!grepl("^C[0-9]+$", rebase)) {
      stop(
        sprintf(
          '`rebase` must be "C1", "C2", ... "C%d" or an integer. Got: "%s".',
          T,
          rebase
        ),
        call. = FALSE
      )
    }
    idx <- as.integer(sub("^C", "", rebase))
  } else {
    idx <- as.integer(rebase)
  }
  if (idx < 1L || idx > T) {
    stop(
      sprintf(
        "`rebase` must be between 1 and %d. Got: %d.",
        T,
        idx
      ),
      call. = FALSE
    )
  }
  idx
}

#' Normalise row/column names of a fitZ$mGamma matrix
#'
#' A plain `multiLCA` object uses `rownames` like `"gamma(Intercept|C)"` and
#' `"gamma(Zp|C)"`. This function strips the `gamma(...)` wrapper so names
#' match the clean format used throughout tseLCA (`"Intercept"`, `"Zp"`, etc.)
#' and ensures column names are `"C2"`, `"C3"`, etc.
#'
#' @param fitZ  A fitZ-like list with at least `$mGamma`.
#' @param Zp.names Character vector of covariate column names (used to set
#'   clean rownames when the raw names can't be parsed). If `NULL`, rownames
#'   are stripped from the `gamma(X|C)` pattern only.
#' @param n_classes Integer. Total number of classes (used to derive clean
#'   column names if they are non-standard).
#' @return `fitZ` with normalised `$mGamma` row/col names.
#' @keywords internal
normalise_fitZ_names <- function(fitZ, Zp.names = NULL, n_classes = NULL) {
  if (is.null(fitZ) || is.null(fitZ$mGamma)) {
    return(fitZ)
  }

  mG <- fitZ$mGamma

  # ---- Normalise rownames ----------------------------------------------------
  rn <- rownames(mG)
  if (!is.null(rn)) {
    # Strip "gamma(X|C)" -> "X"
    rn_clean <- sub("^gamma\\((.+)\\|C\\)$", "\\1", rn)
    # Also handle "gamma(X)" without the |C suffix
    rn_clean <- sub("^gamma\\((.+)\\)$", "\\1", rn_clean)
    rownames(mG) <- rn_clean
  } else if (!is.null(Zp.names)) {
    rownames(mG) <- c("Intercept", Zp.names)
  }

  # ---- Normalise colnames ----------------------------------------------------
  cn <- colnames(mG)
  if (!is.null(cn)) {
    # Already clean ("C2", "C3", ...) -- leave alone
    if (!all(grepl("^C[0-9]+$", cn))) {
      # Non-standard: derive from n_classes if available
      if (!is.null(n_classes)) {
        colnames(mG) <- paste0("C", seq_len(ncol(mG)) + 1L)
      }
    }
  } else if (!is.null(n_classes)) {
    colnames(mG) <- paste0("C", seq_len(ncol(mG)) + 1L)
  }

  fitZ$mGamma <- mG
  fitZ
}

#' Permute class columns of a fitZ object to match a new reference class
#'
#' Rebases a `fitZ` object (output of `fitZ_from_fit0` or
#' `fitZ_from_multiLCA`) so that `ref_idx` becomes the reference class.
#' This involves:
#' \enumerate{
#'   \item Rebasing `$mGamma` -- reconstructing the full T-column log-ratio
#'     matrix, subtracting the new reference column, and dropping it.
#'   \item Propagating through `$Varmat_cor` via the delta method -- the
#'     rebasing transformation is linear (`gamma_new = A * gamma_old`)
#'     so the vcov transforms exactly as `A %*% V %*% t(A)`.
#'   \item Updating all column names.
#' }
#'
#' @param fitZ    Output of `fitZ_from_fit0()` or `fitZ_from_multiLCA()`.
#' @param ref_idx Integer. New reference class (1-based index into the T
#'   classes as currently ordered in `fitZ`).
#' @return `fitZ` with `$mGamma`, `$Varmat_cor`, and names updated.
#' @keywords internal
permute_fitZ_classes <- function(fitZ, ref_idx) {
  if (is.null(fitZ)) {
    return(NULL)
  }

  # Normalise names first so all downstream logic sees clean "Intercept"/"Zp"
  # rownames and "C2"/"C3" colnames regardless of whether fitZ came from
  # fitZ_from_fit0, fitZ_from_multiLCA, or a raw multiLCA call.
  fitZ <- normalise_fitZ_names(fitZ, n_classes = ncol(fitZ$mGamma) + 1L)

  mGamma <- fitZ$mGamma # Q x (T-1): cols = non-ref classes (C2..CT)
  Q <- nrow(mGamma)
  T <- ncol(mGamma) + 1L # total number of classes

  if (ref_idx == 1L) {
    return(fitZ)
  } # already the right reference

  # ---- Step 1: rebase mGamma --------------------------------------------------
  # Reconstruct Q x T full log-ratio matrix (column 1 = 0, reference)
  gamma_full <- cbind(0, mGamma) # Q x T

  # Subtract the new reference column
  new_ref_col <- gamma_full[, ref_idx, drop = FALSE] # Q x 1
  gamma_rebased <- gamma_full - as.vector(new_ref_col) # Q x T, col ref_idx = 0

  # Drop the new reference column; keep remaining classes in ascending order
  keep_cols <- seq_len(T)[-ref_idx] # T-1 indices
  gamma_new <- gamma_rebased[, keep_cols, drop = FALSE]
  colnames(gamma_new) <- paste0("C", keep_cols)
  rownames(gamma_new) <- rownames(mGamma)
  fitZ$mGamma <- gamma_new

  # ---- Step 2: propagate Varmat_cor via delta method --------------------------
  # The rebasing is a linear map on the vec(mGamma) parameter vector.
  # Stacking columns: theta = vec(mGamma), length = Q*(T-1).
  # After rebasing: theta_new = (A_kron_I_Q) * theta_old
  # where A is the (T-1) x (T-1) contrast matrix acting on class columns.
  #
  # A[i, j] = I(keep_cols[i] == old_col_j) - I(ref_idx == old_col_j + 1)
  # Simpler: A = gamma_rebased[, keep_cols] expressed as a linear function
  # of gamma_full[, -1] (the original non-ref columns).
  #
  # gamma_full[, keep_cols] = gamma_full[, keep_cols]
  #                         - gamma_full[, ref_idx] * 1'
  # i.e. A_col[j] = e_{keep_cols[j]-1} - e_{ref_idx-1}  (in the T-1 basis)
  # where e_k is the k-th standard basis vector of R^{T-1},
  # with the convention that the "C1" column has index 0 (not in basis).
  #
  # For the original columns j = 2..T (indexed 1..T-1 in mGamma):
  #   A[i, j] = I(keep_cols[i] - 1 == j)       (direct contribution)
  #           - I(ref_idx - 1   == j)            (subtracted reference)

  if (!is.null(fitZ$Varmat_cor) || !is.null(fitZ$raw_fit$Varmat_cor)) {
    V <- if (!is.null(fitZ$Varmat_cor)) {
      fitZ$Varmat_cor
    } else {
      fitZ$raw_fit$Varmat_cor
    } # Q*(T-1) x Q*(T-1)

    # Build (T-1) x (T-1) column-space contrast matrix A
    # Original columns are indexed 1..(T-1) corresponding to classes C2..CT
    old_non_ref <- seq_len(T - 1L) # 1..(T-1) indexing into mGamma columns
    # keep_cols are class indices (1-based), need to map to mGamma col indices
    keep_mGamma_cols <- keep_cols - 1L # subtract 1 because C1 is not in mGamma
    ref_mGamma_col <- ref_idx - 1L # column of the new ref in old mGamma

    A <- matrix(0, T - 1L, T - 1L)
    for (i in seq_len(T - 1L)) {
      j_direct <- keep_mGamma_cols[i]
      if (j_direct >= 1L && j_direct <= T - 1L) {
        A[i, j_direct] <- 1
      }
      if (ref_mGamma_col >= 1L && ref_mGamma_col <= T - 1L) {
        A[i, ref_mGamma_col] <- A[i, ref_mGamma_col] - 1
      }
    }

    # Kronecker: A acts on class columns, I_Q acts on predictor rows
    # Full transformation: (A kron I_Q)
    A_kron <- kronecker(A, diag(Q)) # Q*(T-1) x Q*(T-1)
    V_new <- A_kron %*% V %*% t(A_kron)

    # Update names
    param_names <- as.vector(outer(
      rownames(mGamma),
      paste0("C", keep_cols),
      paste,
      sep = ":"
    ))
    rownames(V_new) <- param_names
    colnames(V_new) <- param_names

    if (!is.null(fitZ$Varmat_cor)) {
      fitZ$Varmat_cor <- V_new
    }
    if (!is.null(fitZ$raw_fit$Varmat_cor)) {
      fitZ$raw_fit$Varmat_cor <- V_new
    }
  }

  fitZ
}

#' Permute class columns of a fit0 object so that class ref_idx is first
#'
#' Reorders columns of mPhi and vPi so that the desired reference class
#' becomes column 1 before estimation. This ensures the multinomial logit
#' is parameterised with the correct baseline from the start.
#'
#' @param fit0    Raw multilevLCA fit object (has $mPhi and $vPi).
#' @param ref_idx Integer. Class index to move to position 1.
#' @return fit0 with columns permuted.
#' @keywords internal
permute_fit0_classes <- function(fit0, ref_idx) {
  if (ref_idx == 1L) {
    return(fit0)
  }
  T <- ncol(fit0$mPhi)
  ord <- c(ref_idx, seq_len(T)[-ref_idx])
  fit0$mPhi <- fit0$mPhi[, ord, drop = FALSE]
  fit0$vPi <- fit0$vPi[ord]
  fit0
}
