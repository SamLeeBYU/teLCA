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

if (!require("cli")) {
  install.packages("pak")
}

# Install the development version from GitHub
if (!require("pak")) {
  install.packages("pak")
}

if (!require("tseLCA")) {
  pak::pak("SamLeeBYU/tseLCA")
}

library(tseLCA)

###################################################
### generate all simulation conditions
###################################################

dataset_path <- file.path("data/simulation", "sim_datasets.rds")

if (!dir.exists("data/simulation")) {
  dir.create("data/simulation")
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
# Structure mirrors datasets: measurement_models[[scenario]][[sep]][[n]][[rep]]

measurement_path <- file.path("data/simulation", "measurement_models.rds")

if (!file.exists(measurement_path)) {
  message("Fitting measurement models...")

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
            three_step(
              data = reps_data[[r]],
              Y.names = paste0("Y", 1:6),
              n_classes = 3L,
              maxIter.measurement = 4000,
              iter.measurement = 15L,
              R2.threshold = 0.5,
              verbose = FALSE
            ),
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

        measurement_models[[sc]][[sep]][[nn]] <- reps_fit$measurement_model
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
