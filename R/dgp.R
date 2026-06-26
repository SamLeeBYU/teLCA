# teLCA/R/dgp.R
#
# Data-generating process for the simulation study in:
#   Bakk & Kuha (2018). "Two-Step Estimation of Models between Latent Classes
#   and External Variables." Psychometrika, 83(4), 871-892.
#
# Design (Section 3):
#   - 3-class model (T = 3), 6 dichotomous indicators (K = 6)
#   - 3 separation levels: pi = .70 / .80 / .90
#     => entropy R^2 ~ .36 / .65 / .90 (low / medium / high)
#   - 3 sample sizes: n = 500 / 1,000 / 2,000
#   - 2 structural scenarios:
#       "covariate" : Zp (observed) -> X (latent)   [multinomial logistic]
#       "distal"    : X  (latent)   -> Zo (outcome)  [linear regression]
#   => 18 conditions

# -- Default population parameters ---------------------------------------------

#' Default population parameters for the Bakk & Kuha (2018) simulation
#'
#' A list of pre-specified parameters used by [teLCA::generate_data()] and related
#' functions.  All elements correspond to the values stated in the paper
#' (Section 3, p. 879).
#'
#' \describe{
#'   \item{`class_props`}{Length-3 vector of equal class proportions (1/3 each).}
#'   \item{`separation_levels`}{Named vector mapping `"low"`, `"mid"`, `"high"`
#'     to the probability of a "likely" response (0.70, 0.80, 0.90).}
#'   \item{`covariate_params`}{List with `$b0` (intercepts) and `$b` (slopes)
#'     for the multinomial logit P(X=t | Zp).  Intercepts `b02` and `b03` are
#'     set so that marginal class sizes average to 1/3 when Zp ~ Uniform\{1..5\}.}
#'   \item{`distal_params`}{List with `$mu` (class means, c(-1, 1, 0)) and
#'     `$sigma` (residual SD = 1) for the distal outcome model.}
#' }
#'
#' @export
bk2018_params <- list(
  class_props = c(1 / 3, 1 / 3, 1 / 3),

  separation_levels = c(low = 0.70, mid = 0.80, high = 0.90),

  # Intercepts b02 ~ 2.3446, b03 ~ -3.6554 are the unique solution to
  # colMeans(P(X|Zp)) = (1/3, 1/3, 1/3) for Zp ~ Uniform\{1..5\} and slopes
  # b = (0, -1, 1).  Derived once by numerical optimisation; hardcoded here
  # so the package does not run an optimiser at load time.
  covariate_params = list(
    b0 = c(0, 2.3445911086, -3.6553700529), # intercepts (ref class = 0)
    b = c(0, -1, 1) # slopes     (ref class = 0)
  ),

  distal_params = list(
    mu = c(-1, 1, 0), # class-specific means
    sigma = 1 # common residual SD
  )
)


# -- Measurement model ---------------------------------------------------------

#' Build the item-response probability matrix for the simulation
#'
#' Returns a T x K matrix `rho` where `rho[t, k] = P(Y_k = 1 | X = t)`.
#'
#' The three-class structure is:
#' \itemize{
#'   \item Class 1: `pi` on all 6 items (high responders).
#'   \item Class 2: `pi` on items 1-3, `1 - pi` on items 4-6 (mixed).
#'   \item Class 3: `1 - pi` on all 6 items (low responders).
#' }
#'
#' @param pi Numeric scalar in (0.5, 1). Probability of the "likely" response.
#'   Use `bk2018_params$separation_levels` for the three simulation levels.
#'
#' @return A 3 x 6 numeric matrix.
#' @export
make_rho <- function(pi) {
  if (!is.numeric(pi) || length(pi) != 1L || pi <= 0.5 || pi >= 1) {
    stop("`pi` must be a single numeric value in (0.5, 1).", call. = FALSE)
  }
  p_low <- 1 - pi
  rbind(
    class1 = rep(pi, 6),
    class2 = c(rep(pi, 3), rep(p_low, 3)),
    class3 = rep(p_low, 6)
  )
}


# -- Structural model helpers --------------------------------------------------

#' Compute multinomial logistic class probabilities given a covariate
#'
#' Evaluates P(X = t | Zp) for each observation using a multinomial logit
#' with one covariate and class-specific intercepts and slopes.
#'
#' @param Zp     Numeric vector of length n. Covariate values.
#' @param params List with elements `$b0` (length-T intercepts, reference = 0)
#'   and `$b` (length-T slopes, reference = 0).  See `bk2018_params$covariate_params`.
#'
#' @return An n x T matrix of class probabilities (rows sum to 1).
#' @export
mnl_probs <- function(Zp, params) {
  n <- length(Zp)
  T_ <- length(params$b0)
  eta <- outer(Zp, params$b) +
    matrix(params$b0, nrow = n, ncol = T_, byrow = TRUE)
  eta <- eta - apply(eta, 1L, max) # numerical stability
  exp_eta <- exp(eta)
  exp_eta / rowSums(exp_eta) # n x T
}


# -- Drawing functions ---------------------------------------------------------

#' Draw latent class memberships from their marginal distribution
#'
#' @param n  Integer. Sample size.
#' @param pi Numeric vector of length T. Class proportions (must sum to 1).
#'
#' @return Integer vector of length n with values in `1:T`.
#' @export
draw_classes <- function(n, pi) {
  sample(seq_along(pi), size = n, replace = TRUE, prob = pi)
}


#' Draw binary indicators given true class memberships
#'
#' @param X   Integer vector of length n. True latent class (1-indexed).
#' @param rho T x K matrix. `rho[t, k] = P(Y_k = 1 | X = t)`.
#'
#' @return An n x K integer matrix of 0/1 values.
#' @export
draw_indicators <- function(X, rho) {
  n <- length(X)
  K <- ncol(rho)
  Y <- matrix(NA_integer_, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    Y[i, ] <- rbinom(K, size = 1L, prob = rho[X[i], ])
  }
  Y
}


#' Draw the covariate Zp ~ Uniform\{1, 2, 3, 4, 5\}
#'
#' @param n Integer. Sample size.
#'
#' @return Integer vector of length n.
#' @export
draw_Zp <- function(n) {
  sample(1:5, size = n, replace = TRUE)
}


#' Draw latent classes conditional on the covariate (scenario "covariate")
#'
#' @param Zp     Numeric vector of length n. Covariate values.
#' @param params Multinomial logistic parameter list (see `bk2018_params$covariate_params`).
#'
#' @return Integer vector of length n with class labels in `1:T`.
#' @export
draw_classes_given_Zp <- function(Zp, params) {
  probs <- mnl_probs(Zp, params) # n x T
  vapply(
    seq_len(nrow(probs)),
    function(i) sample(ncol(probs), 1L, prob = probs[i, ]),
    integer(1L)
  )
}


#' Draw a continuous distal outcome given true class memberships (scenario "distal")
#'
#' @param X      Integer vector of length n. Latent class (1-indexed).
#' @param params List with `$mu` (length-T class means) and `$sigma` (SD).
#'   See `bk2018_params$distal_params`.
#'
#' @return Numeric vector of length n.
#' @export
draw_Zo <- function(X, params) {
  rnorm(length(X), mean = params$mu[X], sd = params$sigma)
}


# -- Main data-generating function ---------------------------------------------

#' Generate one dataset following the Bakk & Kuha (2018) simulation design
#'
#' @param n          Integer. Sample size (paper uses 500, 1000, or 2000).
#' @param separation Character. One of `"low"`, `"mid"`, `"high"`.
#'   Maps to pi = 0.70, 0.80, 0.90 respectively.
#' @param scenario   Character. One of:
#'   \describe{
#'     \item{`"covariate"`}{Zp (discrete, 1-5) predicts latent X via multinomial logit.}
#'     \item{`"distal"`}{Latent X predicts continuous Zo via linear regression.}
#'   }
#' @param params     List of population parameters.  Defaults to [teLCA::bk2018_params].
#' @param seed       Integer or `NULL`. Optional random seed for reproducibility.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`Y1` .. `Y6`}{Binary indicators (always present).}
#'   \item{`X`}{True latent class, integer 1-3 (not observed in practice).}
#'   \item{`Zp`}{Integer covariate 1-5 (scenario `"covariate"` only).}
#'   \item{`Zo`}{Continuous distal outcome (scenario `"distal"` only).}
#' }
#' @export
generate_data <- function(
  n,
  separation = c("low", "mid", "high"),
  scenario = c("covariate", "distal"),
  params = bk2018_params,
  seed = NULL
) {
  separation <- match.arg(separation)
  scenario <- match.arg(scenario)

  if (!is.numeric(n) || length(n) != 1L || n < 1L || n != round(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  pi_val <- params$separation_levels[[separation]]
  rho <- make_rho(pi_val)

  if (scenario == "covariate") {
    Zp <- draw_Zp(n)
    X <- draw_classes_given_Zp(Zp, params$covariate_params)
    Y <- draw_indicators(X, rho)
    df <- as.data.frame(Y)
    names(df) <- paste0("Y", seq_len(ncol(Y)))
    df$X <- X
    df$Zp <- Zp
  } else {
    X <- draw_classes(n, params$class_props)
    Y <- draw_indicators(X, rho)
    Zo <- draw_Zo(X, params$distal_params)
    df <- as.data.frame(Y)
    names(df) <- paste0("Y", seq_len(ncol(Y)))
    df$X <- X
    df$Zo <- Zo
  }

  df
}


# -- Batch generator -----------------------------------------------------------

#' Generate datasets for all 18 conditions in the simulation design
#'
#' Iterates over the 2 scenarios x 3 separation levels x 3 sample sizes,
#' generating `n_rep` independent replications per condition.  Seeds are
#' derived deterministically from `base_seed` so the entire experiment is
#' reproducible from a single integer.
#'
#' @param n_rep        Integer. Replications per condition (paper uses 500).
#' @param base_seed    Integer. Base seed for reproducibility.
#' @param params       Population parameters list.  Defaults to [teLCA::bk2018_params].
#' @param scenarios    Character. Lists the scenario(s) ("covariate" and/or "distal")
#'   wanting to be simulated. Passed into [generate_data()].
#' @param sep_levels   Character. Lists the separation level(s) ("low", "mid", "high")
#'   wanting to be simulated. Passed into [generate_data()].
#' @param sample_sizes Integer. Lists the sample size(s) wanting to be generated
#'   for each replication condition. Passed into [generate_data()].
#' @param verbose      Logical. If `TRUE` (default), display a live CLI progress
#'   bar with per-rep status and ETA.
#'
#' @return Nested list indexed as
#'   `datasets[[scenario]][[separation]][[as.character(n)]]`,
#'   each element a list of `n_rep` data frames.
#'
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done cli_alert_success
#' @export
generate_all_conditions <- function(
  n_rep = 500L,
  base_seed = 05262026L,
  params = bk2018_params,
  scenarios = c("covariate", "distal"),
  sep_levels = c("low", "mid", "high"),
  sample_sizes = c(500L, 1000L, 2000L),
  verbose = TRUE
) {
  datasets <- list()
  condition_index <- 0L
  n_conditions <- length(scenarios) * length(sep_levels) * length(sample_sizes)
  total_reps <- n_conditions * n_rep

  if (verbose) {
    cli::cli_progress_bar(
      name = "Generating simulation datasets",
      total = total_reps,
      format = paste0(
        "{cli::pb_spin} {cli::pb_name} | ",
        "{cli::pb_bar} {cli::pb_percent} | ",
        "Rep {cli::pb_current}/{cli::pb_total} | ",
        "Elapsed: {cli::pb_elapsed} | ETA: {cli::pb_eta}"
      ),
      format_done = paste0(
        "{cli::pb_name} | ",
        "{cli::pb_total} reps across {n_conditions} conditions | ",
        "Total time: {cli::pb_elapsed}"
      )
    )
  }

  for (sc in scenarios) {
    datasets[[sc]] <- list()
    for (sep in sep_levels) {
      datasets[[sc]][[sep]] <- list()
      for (nn in sample_sizes) {
        condition_index <- condition_index + 1L
        key <- as.character(nn)
        reps <- vector("list", n_rep)

        for (r in seq_len(n_rep)) {
          if (verbose) {
            cli::cli_progress_update(
              status = sprintf(
                "condition %2d/%2d | scenario=%-10s sep=%-4s n=%-5d rep=%d/%d",
                condition_index,
                n_conditions,
                sc,
                sep,
                nn,
                r,
                n_rep
              )
            )
          }

          seed_r <- base_seed * 1e6 + condition_index * 1e3 + r
          reps[[r]] <- generate_data(
            n = nn,
            separation = sep,
            scenario = sc,
            params = params,
            seed = as.integer(seed_r %% .Machine$integer.max)
          )
        }

        datasets[[sc]][[sep]][[key]] <- reps
      }
    }
  }

  if (verbose) {
    cli::cli_progress_done()
  }

  datasets
}
