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

USE_PARALLEL <- FALSE
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
} else {
  future::plan(future::sequential)
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
#' @param historical_data Pre-loaded and prepared historical data for this survey
#'
#' @return Combined data ready for model fitting
prepare_combined_data <- function(species, survey_abbrev, sim_data, historical_data) {
  # Prepare simulated data
  sim_data_prep <- sim_data |>
    mutate(
      catch_count = observed,
      historical = FALSE
    )

  # Combine and create final structure
  combined_data <- bind_rows(historical_data, sim_data_prep) |>
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
  select(rate_case, rate, survey, plan, filename, filepath)

# =============================================================================
# Pre-load historical data for all surveys (sf joins expensive on hake)
# =============================================================================
# Directory with pre-processed historical data
hist_data_dir <- here::here("data-generated", "historical-data-processed")

# Get unique surveys from scenario files
unique_surveys <- unique(scenario_files$survey)

# Load pre-processed historical data for each survey
historical_data_list <- purrr::map(unique_surveys, function(survey_abbrev) {
cache_file <- file.path(hist_data_dir, paste0(sp_to_hyphens(species), "-", survey_abbrev, "-processed.rds"))

  if (!file.exists(cache_file)) {
    stop("Pre-processed data not found: ", cache_file, "\n",
         "Run R/00-preprocess-historical-data.R on your local machine first!")
  }

  message("  Loading: ", survey_abbrev)
  hdat <- readRDS(cache_file)
  message("    ", nrow(hdat), " rows, ", sum(hdat$restricted), " in MPAs")

  return(hdat)
}) |>
  setNames(unique_surveys)

message("  Loaded ", length(historical_data_list), " historical datasets\n")

# =============================================================================
# Check existing caches for resume capability
# =============================================================================

message("=== Checking Existing Caches ===")

# Add cache filenames and check completion status
scenario_files <- scenario_files |>
  mutate(
    plan_clean = gsub("-", "_", plan),
    rate_clean = gsub("\\.", "", as.character(rate)),
    cache_file = file.path(
      results_dir,
      paste0(
        sp_to_hyphens(species), "-",
        survey, "-",
        "mpa", rate_clean, "-",
        ar1_scenario, "-",
        time_scenario, "-",
        plan_clean,
        "-trend-estimates.rds"
      )
    ),
    cache_exists = file.exists(cache_file)
  )

# For each scenario, determine which replicates are already complete
scenario_files <- scenario_files |>
  rowwise() |>
  mutate(
    completed_reps = list(
      if (cache_exists) {
        cached_results <- readRDS(cache_file)
        unique(cached_results$replicate)
      } else {
        integer(0)
      }
    ),
    n_complete = length(completed_reps),
    n_remaining = N_REPLICATES - n_complete
  ) |>
  ungroup()

# Summary of resume status
total_complete <- sum(scenario_files$n_complete)
total_remaining <- sum(scenario_files$n_remaining)
message("  Complete: ", total_complete, " replicate fits")
message("  Remaining: ", total_remaining, " replicate fits")
message("  Resume: ", sum(scenario_files$n_complete > 0), " scenarios have partial results\n")

message("=== Starting Power Analysis Fitting ===")
message("Species: ", species)
message("Found ", nrow(scenario_files), " scenario files to process")
message("  Rates: ", length(unique(scenario_files$rate)))
message("  Surveys: ", length(unique(scenario_files$survey)))
message("  Plans: ", length(unique(scenario_files$plan)))
message("Processing ", N_REPLICATES, " replicates per scenario")
message("Total jobs: ", total_remaining, " (", total_complete, " already complete)\n")

# =============================================================================
# Create flattened job grid (scenario × incomplete replicates)
# =============================================================================

message("=== Creating Job Grid ===")

# Expand scenario_files to create one row per incomplete replicate
job_grid <- scenario_files |>
  rowwise() |>
  mutate(
    # Get incomplete replicate numbers for this scenario
    incomplete_reps = list(setdiff(1:N_REPLICATES, completed_reps))
  ) |>
  ungroup() |>
  # Filter to scenarios with incomplete work
  filter(n_remaining > 0) |>
  # Expand to one row per incomplete replicate
  tidyr::unnest(incomplete_reps) |>
  rename(replicate = incomplete_reps) |>
  # Add job ID for progress tracking
  mutate(job_id = row_number())

message("  Created ", nrow(job_grid), " jobs\n")

# Exit early if no work to do
if (nrow(job_grid) == 0) {
  message("=== All scenarios complete! ===")
  # Load existing results
  all_fitted_results <- purrr::map_dfr(scenario_files$cache_file, readRDS)
} else {

  # =============================================================================
  # Pre-load sampled data for incomplete scenarios
  # =============================================================================

  message("=== Pre-loading Sampled Data ===")

  unique_filepaths <- unique(job_grid$filepath)
  message("  Loading ", length(unique_filepaths), " scenario files...")

  sampled_data_list <- purrr::map(unique_filepaths, function(fp) {
    readRDS(fp)
  }) |>
    setNames(unique_filepaths)

  message("  Loaded sampled data for ", length(sampled_data_list), " scenarios\n")

  # =============================================================================
  # Process all incomplete jobs in parallel
  # =============================================================================

  message("=== Processing Jobs ===")
  message("  Running ", nrow(job_grid), " jobs across ", if (USE_PARALLEL) N_WORKERS else 1, " workers\n")

  # Define the worker function (same for both parallel and sequential)
  process_job <- function(i) {
    job <- job_grid[i, ]

    # Get data for this job
    sim_data_all_reps <- sampled_data_list[[job$filepath]]
    rep_sim_data <- sim_data_all_reps |> filter(replicate == job$replicate)
    hist_data <- historical_data_list[[job$survey]]

    # Combine data
    combined_data <- prepare_combined_data(
      species = species,
      survey_abbrev = job$survey,
      sim_data = rep_sim_data,
      historical_data = hist_data
    )

    # Fit model
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
      survey_abbrev = job$survey,
      sim_mpa_trend = job$rate,
      sim_mpa_case = job$rate_case,
      sim_ar1_scenario = ar1_scenario,
      sim_time_scenario = time_scenario,
      plan = job$plan,
      replicate = job$replicate,
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
  }

  # Run jobs (parallel or sequential)
  if (USE_PARALLEL) {
    new_results <- furrr::future_map_dfr(
      1:nrow(job_grid),
      process_job,
      .options = furrr::furrr_options(
        seed = TRUE,
        packages = c("sdmTMB", "dplyr", "tibble", "broom.mixed", "sf"),
        globals = c("fit_simulation", "extract_trend_estimate",
                    "summarise_sanity", "clean_family_name", "prepare_combined_data",
                    "sp_to_hyphens", "XY_to_sf", "species", "ar1_scenario", "time_scenario",
                    "sampled_data_list", "historical_data_list", "job_grid", "process_job")
      ),
      .progress = TRUE
    )
  } else {
    new_results <- purrr::map_dfr(1:nrow(job_grid), process_job, .progress = TRUE)
  }

  message("\n=== Merging and Caching Results ===")

  # Group new results by scenario
  new_results_df <- bind_rows(new_results)

  # Get unique scenario identifiers from job_grid
  scenarios_to_update <- job_grid |>
    distinct(rate, rate_case, survey, plan, cache_file)

  # For each scenario, merge new results with existing cache (if any) and save
  all_fitted_results <- purrr::map_dfr(1:nrow(scenarios_to_update), function(i) {
    scenario <- scenarios_to_update[i, ]

    # Get new results for this scenario
    scenario_new_results <- new_results_df |>
      filter(
        sim_mpa_trend == scenario$rate,
        survey_abbrev == scenario$survey,
        plan == scenario$plan
      )

    # Merge with existing results if cache exists
    if (file.exists(scenario$cache_file)) {
      scenario_existing_results <- readRDS(scenario$cache_file)
      scenario_all_results <- bind_rows(scenario_existing_results, scenario_new_results)
      n_new <- nrow(scenario_new_results)
      n_existing <- nrow(scenario_existing_results)
      message("  ", basename(scenario$cache_file), ": ", n_new, " new + ", n_existing, " existing = ", nrow(scenario_all_results), " total")
    } else {
      scenario_all_results <- scenario_new_results
      message("  ", basename(scenario$cache_file), ": ", nrow(scenario_all_results), " new")
    }

    # Save updated cache
    saveRDS(scenario_all_results, scenario$cache_file)

    return(scenario_all_results)
  })

  # Also add any scenarios that were already 100% complete
  complete_scenarios <- scenario_files |>
    filter(n_remaining == 0)

  if (nrow(complete_scenarios) > 0) {
    complete_results <- purrr::map_dfr(complete_scenarios$cache_file, readRDS)
    all_fitted_results <- bind_rows(all_fitted_results, complete_results)
  }
}

message("=== All fitting complete ===")
message("Fitted results cached in: ", results_dir)

# Save combined results
saveRDS(all_fitted_results, file.path(results_dir, "all-fitted-results.rds"))
message("Combined results saved: ", file.path(results_dir, "all-fitted-results.rds"))

# Reset to sequential processing
future::plan(future::sequential)
