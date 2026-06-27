#####################################################
### Simulation Script for the tseLCA package:
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee
### E-Mail: samlee@arizona.edu
#####################################################

###################################################
### preliminaries
###################################################

rm(list = ls())
gc()
r_opts <- options(
  prompt = "R> ",
  continue = "+  ",
  width = 77,
  digits = 4,
  useFancyQuotes = FALSE,
  warn = 1
)

# Loading libraries and installing if unavailable

# Install the development version from GitHub
if (!require("tseLCA")) {
  if (!require("pak")) {
    install.packages("pak")
  }
  pak::pak("SamLeeBYU/tseLCA")
}

library(tseLCA)

###################################################
### generate all simulation conditions
###################################################

output.dir <- "tseLCA_output/simulation"
dataset_path <- file.path(output.dir, "sim_datasets.rds")

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

if (!file.exists(dataset_path)) {
  message("Generating simulation data (this will take a few minutes)...")

  datasets <- generate_all_conditions(
    n_rep = 500L,
    base_seed = 06262026L,
    sep_levels = c("low", "mid", "high")
  )
  saveRDS(datasets, file = dataset_path)
  message("Completed data generation and saved datasets to: ", dataset_path)
} else {
  message("Loading existing simulation data from: ", dataset_path)
  datasets <- readRDS(dataset_path)
}

###################################################
### obtain measurement models for all conditions
###################################################

# Pre-computing measurement models for each replicate saves time when
# comparing estimators (BCH vs ML, proportional/modal assignment) because
# lca_step1() is the computational bottleneck. Every dataset is unique
# (different seed per replicate x condition), so no measurement models
# can be shared.
#
# We call multiLCA to obtain two-step vcov estimates.
# While multiLCA cannot accomodate fixed parameter values as input,
# we create covariate-adjusted predictions for the latent class to obtain
# consistent results for multiLCA in the presence of low-separation.
#
# The alternative would be to call multiLCA a bunch of times and take the model
# with the best log-likelihood (which is what we already do for the measurement model).
# This saves some computation time.
#
# Structure mirrors datasets: measurement_models[[scenario]][[sep]][[n]][[rep]]

measurement_path <- file.path(output.dir, "measurement_models.rds")

if (!file.exists(measurement_path)) {
  message(
    "Fitting measurement models... This will take several hours due to low-separation cases."
  )

  scenarios <- c("covariate", "distal")
  sep_levels <- c("low", "mid", "high")
  sample_sizes <- c("500", "1000", "2000")

  total_conditions <- length(scenarios) *
    length(sep_levels) *
    length(sample_sizes)
  n_rep <- length(datasets[[scenarios[1]]][[sep_levels[1]]][[sample_sizes[1]]])
  total_reps <- total_conditions * n_rep

  measurement_models <- list()

  cli::cli_progress_bar(
    name = "Fitting measurement models",
    total = total_reps,
    format = paste0(
      "{cli::pb_name} | {cli::pb_bar} {cli::pb_percent} | ",
      "Rep {cli::pb_current}/{cli::pb_total} | ",
      "Elapsed: {cli::pb_elapsed} | ETA: {cli::pb_eta}"
    )
  )

  mGamma.init <- do.call(rbind, bk2018_params$covariate_params)[, -1L]
  p.xz.sim <- function(Z_mat, params) {
    eta_full <- cbind(0, Z_mat %*% params)
    row_max <- apply(eta_full, 1, max)
    exp_eta <- exp(eta_full - row_max)
    exp_eta / rowSums(exp_eta)
  }
  for (sc in scenarios) {
    measurement_models[[sc]] <- list()
    for (sep in sep_levels) {
      measurement_models[[sc]][[sep]] <- list()
      for (nn in sample_sizes) {
        reps_data <- datasets[[sc]][[sep]][[nn]]
        n_rep <- length(reps_data)
        reps_fit <- vector("list", n_rep)

        for (r in seq_len(n_rep)) {
          cli::cli_progress_update(
            status = sprintf(
              "scenario=%-10s sep=%-4s n=%-5s rep=%d/%d",
              sc,
              sep,
              nn,
              r,
              n_rep
            )
          )

          reps_fit[[r]] <- tryCatch(
            {
              m.r <- three_step(
                data = reps_data[[r]],
                Y.names = paste0("Y", 1:6),
                n_classes = 3L,
                maxIter.measurement = 5000,
                iter.measurement = 20L,
                R2.threshold = 0.6,
                verbose = FALSE
              )$measurement_model
              if (sc == "covariate") {
                c.fitZ <- fitZ_from_fit0(
                  fit0 = m.r$fit0,
                  data = reps_data[[r]],
                  Y.names = paste0("Y", 1:6),
                  Zp.names = "Zp",
                  maxIter = 500,
                  starting_val = mGamma.init
                )

                pi_adj <- p.xz.sim(cbind(1, reps_data[[r]]$Zp), c.fitZ$mGamma)

                Y_mat <- as.matrix(reps_data[[r]][, paste0("Y", 1:6)])
                mPhi.init <- m.r$fit0$mPhi
                log_lik_items <- Y_mat %*%
                  log(mPhi.init) +
                  (1 - Y_mat) %*% log(1 - mPhi.init)
                log_W <- log(pi_adj) + log_lik_items
                log_W <- log_W -
                  apply(log_W, 1, function(x) {
                    m <- max(x)
                    m + log(sum(exp(x - m)))
                  })
                W_init <- exp(log_W)

                #Covariate-adjusted classes
                reps_data[[r]]$startval <- apply(W_init, 1, which.max)
                has_all_classes <- length(unique(reps_data[[r]]$startval)) == 3L
                c.r <- multilevLCA::multiLCA(
                  data = reps_data[[r]],
                  Y = paste0("Y", 1:6),
                  iT = 3L,
                  Z = "Zp",
                  startval = if (has_all_classes) "startval" else NULL,
                  extout = TRUE,
                  verbose = FALSE
                )
                #c.r$LLKSeries[nrow(c.r$LLKSeries)]
                m.r$fitZ <- c.r

                m.r$fitZ_converged <- tail(c.r$LLKSeries, 2) |>
                  diff() |>
                  abs() <
                  1e-8
                m.r$fitZ_iters <- c.r$iter
              }
              m.r
            },
            error = function(e) {
              cli::cli_alert_warning(sprintf(
                "three_step failed: scenario=%s sep=%s n=%s rep=%d: %s",
                sc,
                sep,
                nn,
                r,
                conditionMessage(e)
              ))
              NULL
            }
          )
        }

        measurement_models[[sc]][[sep]][[nn]] <- reps_fit
      }
    }
  }

  cli::cli_progress_done()

  saveRDS(measurement_models, file = measurement_path)
  cli::cli_alert_success("Measurement models saved to: {measurement_path}")
} else {
  cli::cli_alert_info(
    "Loading existing measurement models from: {measurement_path}"
  )
  measurement_models <- readRDS(measurement_path)
}
