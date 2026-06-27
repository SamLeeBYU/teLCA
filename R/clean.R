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
