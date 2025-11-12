# =============================================================================
# Fit sampled data and analyze power
# =============================================================================
# This script fits models to sampled data and evaluates power to detect
# MPA recovery trends under different scenarios.

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

library(purrr)

# =============================================================================
# Configuration
# =============================================================================

cleaned_data_dir <- here::here("data-generated", "cleaned-species-data")
sample_dir <- here::here("data-generated", "sampled-data")
results_dir <- here::here("data-generated", "power-results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- TRUE
#N_WORKERS <- if (USE_PARALLEL) 20 else 1
N_WORKERS <- NULL

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 40
}

if (Sys.info()['user'] == "jilliandunic") {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 8
}

# Number of replicates to process (default 10 for testing, set to 50 for full run)
N_REPLICATES <- 10

# Setup parallel processing
if (USE_PARALLEL) {
  if (is.null(N_WORKERS)) N_WORKERS <- future::availableCores() / 2

  if (Sys.info()['user'] %in% c("dunic", "anderson")) {
    future::plan(future::multicore, workers = N_WORKERS)
    message("Using ", N_WORKERS, " parallel workers (multicore)")
  } else {
    future::plan(future::multisession, workers = N_WORKERS)
    message("Using ", N_WORKERS, " parallel workers (multisession)")
  }
  map_fn <- furrr::future_map_dfr
} else {
  future::plan(future::sequential)
  map_fn <- purrr::map_dfr
  message("Using sequential processing")
}

#' Fit sdmTMB model to sampled data
#'
#' @param dat Sampled data for one replicate
#' @param formula Model formula
#' @param spatial Spatial random field specification
#' @param spatiotemporal Spatiotemporal random field specification
#' @param family Distribution family
#' @param silent Suppress sdmTMB messages
#'
#' @return Fitted sdmTMB model or error object
fit_simulation <- function(dat,
                           formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
                           spatial = "on",
                           spatiotemporal = "iid",
                           family = betabinomial(link = "cloglog"),
                           cutoff = 10,
                           control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
                           silent = TRUE) {

  survey_type <- unique(dat$survey_abbrev)

  # Prepare data
  if (grepl("HBLL", survey_type)) { # future proofing to allow use of SYN surveys
    weights <- dat$hook_count
    offset <- NULL
  } else {
    weights <- NULL
    offset <- dat$offset
  }
  # Create mesh
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  # Fit model
  fit <- local(tryCatch({
    sdmTMB(
      formula = formula,
      data = dat,
      mesh = mesh,
      family = family,
      spatial = spatial,
      spatiotemporal = spatiotemporal,
      time = "year",
      weights = weights, # if HBLL use weights, otherwise use offset
      offset = offset,
      silent = silent,
      control = control
    )
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  }))

  return(fit)
}

#' Extract MPA trend estimate from fitted model
#'
#' @param fit Fitted sdmTMB model
#' @param trend_param Name of trend parameter in model
#'
#' @return List with estimate, se, and confidence interval
extract_trend_estimate <- function(fit, trend_param = "restricted:year_covariate") {
  if (!is.null(fit$error) && fit$error) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      converged = FALSE,
      sanity = NA_character_,
      error_msg = fit$message
    ))
  }

  # Extract coefficient
  coefs <- tidy(fit, conf.int = TRUE)
  trend_row <- coefs |> filter(term == trend_param)

  if (nrow(trend_row) == 0) {
    stop("Parameter ", trend_param, " not found in model")
  }

  return(list(
    estimate = trend_row$estimate,
    se = trend_row$std.error,
    ci_lower = trend_row$conf.low,
    ci_upper = trend_row$conf.high,
    converged = TRUE,
    sanity = summarise_sanity(fit),
    error_msg = NA_character_
  ))
}

#' Prepare combined historical and simulated data for model fitting
#'
#' @param species Species name
#' @param survey_abbrev Survey abbreviation (e.g., "HBLL-OUT-N")
#' @param sim_data Simulated data for one replicate
#' @param cleaned_data_dir Directory with cleaned historical data
#'
#' @return Combined data ready for model fitting
prepare_combined_data <- function(species, survey_abbrev, sim_data, cleaned_data_dir) {
  # Load MPA spatial data
  simple_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa.rds"))

  # Load historical data
  hdat0 <- readRDS(file.path(cleaned_data_dir, paste0(sp_to_hyphens(species), "-", survey_abbrev, ".rds")))

  # Join with MPA and prepare
  hdat <- st_join(XY_to_sf(hdat0, crs_to = st_crs(simple_mpa)), simple_mpa, join = st_within) |>
    mutate(
      restricted = ifelse(is.na(uid), 0, 1),
      last_sampled_year = max(year),
      year_covariate = 0,
      historical = TRUE
    ) |>
    st_drop_geometry()

  # Prepare simulated data
  sim_data_prep <- sim_data |>
    mutate(
      catch_count = observed,
      historical = FALSE
    )

  # Combine and create final structure
  combined_data <- bind_rows(hdat, sim_data_prep) |>
    select(survey_abbrev, X, Y, block_id, restricted, historical,
           year, year_covariate, last_sampled_year,
           catch_count, hook_count, offset) |>
    mutate(
      catch_prop = catch_count / hook_count,
      fyear_value = ifelse(historical, year, last_sampled_year),
      fyear = as.factor(fyear_value)
    )

  return(combined_data)
}


# # Testing
# species <- "yelloweye rockfish"
# f <- list.files(file.path(sample_dir, sp_to_hyphens(species)))
#
# scenario <- "HBLL-OUT-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo.rds"
#
# sim_dat0 <- readRDS(file.path(sample_dir, sp_to_hyphens(species), scenario)) |>
#   mutate(catch_count = observed, historical = FALSE)
#
# sim_dat <- sim_dat0 |> filter(replicate == 1)
#
# simple_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa.rds"))
# hdat0 <- readRDS(file.path(cleaned_data_dir, paste0(sp_to_hyphens(species), "-HBLL-OUT-N.rds")))
# hdat <- st_join(XY_to_sf(hdat0, crs_to = st_crs(simple_mpa)), simple_mpa, join = st_within) |>
#   mutate(restricted = ifelse(is.na(uid), 0, 1),
#          last_sampled_year = max(year),
#          year_covariate = 0,
#          historical = TRUE) |>
#   st_drop_geometry()
#
# d <- bind_rows(hdat, sim_dat) |>
#   select(survey_abbrev, X, Y, block_id, restricted, historical,
#          year, year_covariate, last_sampled_year, # year_covariate is time since implementation
#          catch_count, hook_count, offset) |>
#   mutate(catch_prop = catch_count / hook_count,
#          weights = hook_count / mean(hook_count),
#          fyear_value = ifelse(historical, year, last_sampled_year), # last sampled year should be the intercept for the future simulated years
#          fyear = as.factor(fyear_value))
#
# test <- fit_simulation(
#   dat = d,
#   formula = catch_prop ~ fyear + restricted + year_covariate + restricted:year_covariate,
#   spatial = "on",
#   spatiotemporal = "iid",
#   cutoff = 10,
#   control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
#   silent = FALSE
#   )
# meep()
# sanity(test)
# test
#
# extract_trend_estimate(test)

# I want to load in the models that fit these scenarios
# Then extract the trend estimates for each model
# I want to do this for each of the three surveys, so this is the part of the
# filenames I will probably need to capture: _mpa1.01509500351865_no_AR1_twenty_years_status-quo.rds
# "HBLL-INS-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo.rds"
# "HBLL-INS-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo-20-effort.rds"
# "HBLL-INS-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo-40-effort-every-5-years.rds"
# "HBLL-INS-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo-40-effort.rds"
# "HBLL-INS-N_mpa1.01509500351865_no_AR1_twenty_years_status-quo-no-sampling-in-mpas.rds"

# I want to fit the model to each of the replicates and extract the trend estimates
# I want to make this efficient to run on a 40 core server


# =============================================================================
# Main Execution
# =============================================================================

# Load recovery rates
recovery_rates <- readRDS(file.path("data-generated", "recovery-rates.rds"))

# Define scope
species <- "yelloweye rockfish"
surveys <- c("HBLL-OUT-N", "HBLL-OUT-S", "HBLL-INS-N")
plans <- c(
  "status-quo",
  "status-quo-20-effort",
  "status-quo-40-effort-every-5-years",
  "status-quo-40-effort",
  "status-quo-no-sampling-in-mpas"
)

# Get MPA rates for this species
species_rates <- recovery_rates |>
  filter(species == !!species) |>
  select(species, case, linear_mpa_rate)

message("\n=== Recovery rates for ", species, " ===")
print(species_rates)

# Build scenario file list by finding files that match patterns
species_dir <- file.path(sample_dir, sp_to_hyphens(species))
ar1_scenario <- "no_AR1"
time_scenario <- "twenty_years"

# Get all files in species directory
all_files <- list.files(species_dir, pattern = "\\.rds$", full.names = FALSE)

# Create expected scenarios and match to existing files
scenario_files <- tidyr::expand_grid(
  rate_case = species_rates$case,
  rate = species_rates$linear_mpa_rate,
  survey = surveys,
  plan = plans
) |>
  mutate(
    # Create pattern to match files (rate might have different precision in filename)
    file_pattern = paste0(
      "^", survey, "_mpa", rate, ".*_", ar1_scenario, "_", time_scenario, "_", plan, "\\.rds$"
    ),
    # Find matching file
    filename = purrr::map_chr(file_pattern, function(pattern) {
      matches <- grep(pattern, all_files, value = TRUE)
      if (length(matches) > 0) return(matches[1])
      return(NA_character_)
    }),
    filepath = file.path(species_dir, filename)
  ) |>
  filter(!is.na(filename)) |>
  select(rate_case, rate, survey, plan, filename, filepath) |>
  slice(1)

message("\n=== Starting Power Analysis Fitting ===")
message("Species: ", species)
message("Found ", nrow(scenario_files), " scenario files to process")
message("  Rates: ", length(unique(scenario_files$rate)))
message("  Surveys: ", length(unique(scenario_files$survey)))
message("  Plans: ", length(unique(scenario_files$plan)))
message("Processing ", N_REPLICATES, " replicates per scenario\n")

# Process each scenario
all_fitted_results <- purrr::map_dfr(1:nrow(scenario_files), function(i) {
  scenario <- scenario_files[i, ]

  message("========================================")
  message("Scenario ", i, "/", nrow(scenario_files))
  message("Rate: ", scenario$rate, " (", scenario$rate_case, ")")
  message("Survey: ", scenario$survey, " | Plan: ", scenario$plan)
  message("========================================")

  # Create cache filename
  plan_clean <- scenario$plan |>
    gsub("-", "_", x = _)
  rate_clean <- gsub("\\.", "", as.character(scenario$rate))

  cache_file <- file.path(
    results_dir,
    paste0(
      sp_to_hyphens(species), "-",
      scenario$survey, "-",
      "mpa", rate_clean, "-",
      ar1_scenario, "-",
      time_scenario, "-",
      plan_clean,
      "-trend-estimates.rds"
    )
  )

  # Check cache
  if (file.exists(cache_file)) {
    message("  Cache hit\n")
    return(readRDS(cache_file))
  }

  message("  Loading sampled data...")
  # Load sampled data (all replicates)
  sampled_data <- readRDS(scenario$filepath)

  # Subset replicates
  available_reps <- unique(sampled_data$replicate)
  selected_reps <- head(available_reps, N_REPLICATES)

  message("  Processing ", length(selected_reps), " of ", length(available_reps), " replicates")

  # Fit each replicate (in parallel if enabled)
  rep_results <- map_fn(selected_reps, function(rep) {
    # Filter to this replicate
    rep_sim_data <- sampled_data |> filter(replicate == rep)

    # Combine with historical data
    combined_data <- prepare_combined_data(
      species = species,
      survey_abbrev = scenario$survey,
      sim_data = rep_sim_data,
      cleaned_data_dir = cleaned_data_dir
    )

    # Fit model with collapse_spatial_variance = TRUE (handles collapsed fields automatically)
    fit <- fit_simulation(
      dat = combined_data,
      formula = catch_prop ~ fyear + restricted + year_covariate + restricted:year_covariate,
      spatial = "on",
      spatiotemporal = "iid",
      cutoff = 10,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
      silent = TRUE
    )

    # Extract trend
    trend <- extract_trend_estimate(fit)

    # Return results
    tibble(
      species = species,
      survey_abbrev = scenario$survey,
      sim_mpa_trend = scenario$rate,
      sim_mpa_case = scenario$rate_case,
      sim_ar1_scenario = ar1_scenario,
      sim_time_scenario = time_scenario,
      plan = scenario$plan,
      replicate = rep,
      estimate = trend$estimate,
      se = trend$se,
      ci_lower = trend$ci_lower,
      ci_upper = trend$ci_upper,
      converged = trend$converged,
      sanity = trend$sanity,
      error_msg = trend$error_msg,
      fit_family = clean_family_name(fit),
      fit_spatial = fit$spatial,
      fit_spatiotemporal = fit$spatiotemporal,
      formula = deparse(fit$formula)
    )
  }, .options = if (USE_PARALLEL) {
    furrr::furrr_options(
      seed = TRUE,
      packages = c("sdmTMB", "dplyr", "tibble", "broom.mixed", "sf"),
      globals = c("fit_simulation", "extract_trend_estimate",
                  "summarise_sanity", "clean_family_name", "prepare_combined_data",
                  "sp_to_hyphens", "XY_to_sf", "cleaned_data_dir", "species")
    )
  } else {
    list()
  })

  # Save cache
  saveRDS(rep_results, cache_file)
  message("  Saved: ", basename(cache_file), "\n")

  return(rep_results)
})

message("=== All fitting complete ===")
message("Fitted results cached in: ", results_dir)

# Save combined results
saveRDS(all_fitted_results, file.path(results_dir, "all-fitted-results.rds"))
message("Combined results saved: ", file.path(results_dir, "all-fitted-results.rds"))

# Reset to sequential processing
future::plan(future::sequential)
