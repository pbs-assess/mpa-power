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

sample_dir <- here::here("data-generated", "sampled-data")
results_dir <- here::here("data-generated", "power-results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- TRUE
#N_WORKERS <- if (USE_PARALLEL) 20 else 1
N_WORKERS <- NULL

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

# =============================================================================
# Helper Functions
# =============================================================================

#' Fit sdmTMB model to sampled data
#'
#' @param sampled_dat Sampled data for one replicate
#' @param formula Model formula
#' @param spatial Spatial random field specification
#' @param spatiotemporal Spatiotemporal random field specification
#' @param family Distribution family
#' @param silent Suppress sdmTMB messages
#'
#' @return Fitted sdmTMB model or error object
fit_simulation <- function(sampled_dat,
                           formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
                           spatial = "on",
                           spatiotemporal = "iid",
                           family = betabinomial(link = "cloglog"),
                           cutoff = 10,
                           silent = TRUE) {

  survey_type <- unique(sampled_dat$survey_abbrev)

  # Prepare data
  if (grepl("HBLL", survey_type)) { # future proofing to allow use of SYN surveys
    dat <- sampled_dat |>
      mutate(
        obs_id = factor(row_number()),
        catch_prop = observed / hook_count,
        fyear = as.factor(year)
      )
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
      weights = dat$hook_count,
      silent = silent
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


# Get list of sampled data files (species-level)
sampled_files <- list.files(sample_dir, pattern = "-all-sampled\\.rds$", full.names = TRUE)

message("Found ", length(sampled_files), " sampled species files")

# Process each species file
all_fitted_results <- purrr::map_dfr(sampled_files, function(file) {
  species <- tools::file_path_sans_ext(basename(file)) |>
    gsub("-all-sampled", "", x = _) |>
    sp_from_hyphens()

  message("\n========================================")
  message("Fitting simulation for species: ", species)
  message("========================================")

  # Load sampled data (all params × plans for this species)
  sampled_data <- readRDS(file)

  # Get unique param × plan combinations
  param_plan_combos <- sampled_data |>
    distinct(survey_abbrev, sim_mpa_trend, sim_ar1_scenario, sim_time_scenario, plan)

  message("  ", nrow(param_plan_combos), " param × plan combinations")

  # Process each combination
  purrr::map_dfr(1:nrow(param_plan_combos), function(i) {
    combo <- param_plan_combos[i, ]

    message("  Combo ", i, "/", nrow(param_plan_combos),
            ": mpa=", combo$sim_mpa_trend,
            ", ar1=", combo$sim_ar1_scenario,
            ", time=", combo$sim_time_scenario,
            ", plan=", combo$plan)

    # Create cache filename (sanitize plan name for filesystem)
    plan_clean <- combo$plan |>
      gsub(" ", "_", x = _) |>
      gsub("\\+", "plus", x = _) |>
      gsub("%", "pct", x = _) |>
      gsub("[^a-zA-Z0-9_-]", "", x = _)

    cache_file <- file.path(
      results_dir,
      paste0(
        sp_to_hyphens(species), "-",
        "mpa", combo$sim_mpa_trend, "-",
        combo$sim_ar1_scenario, "-",
        combo$sim_time_scenario, "-",
        plan_clean,
        "-trend-estimates.rds"
      )
    )

    if (file.exists(cache_file)) {
      message("    Cache hit")
      return(readRDS(cache_file))
    }

    message("    Fitting...")
    # Filter to this combination
    combo_data <- sampled_data |>
      filter(
        survey_abbrev == combo$survey_abbrev,
        sim_mpa_trend == combo$sim_mpa_trend,
        sim_ar1_scenario == combo$sim_ar1_scenario,
        sim_time_scenario == combo$sim_time_scenario,
        plan == combo$plan
      )

    replicates <- unique(combo_data$replicate)

    # Fit each replicate (in parallel if enabled)
    rep_results <- map_fn(replicates, function(rep) {
      rep_data <- combo_data |> filter(replicate == rep)

      # Fit model for replicate
      fit <- fit_simulation(rep_data, silent = TRUE)

      # Turn off and refit collapsed random fields
      if (inherits(fit, "sdmTMB")) {
        collapse_check <- update_collapsed_rf(fit)
        if (collapse_check$needs_refit) {
          message("        Refitting with collapsed fields turned off: ", combo$survey_abbrev)
          fit <- fit_simulation(
            rep_data,
            spatial = collapse_check$spatial,
            spatiotemporal = collapse_check$spatiotemporal,
            silent = TRUE
          )
        }
      }

      # Extract trend
      trend <- extract_trend_estimate(fit)

      # Return results
      tibble(
        species = species,
        survey_abbrev = combo$survey_abbrev,
        sim_mpa_trend = combo$sim_mpa_trend,
        sim_ar1_scenario = combo$sim_ar1_scenario,
        sim_time_scenario = combo$sim_time_scenario,
        plan = combo$plan,
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
        packages = c("sdmTMB", "dplyr", "tibble", "broom.mixed"),
        globals = c("fit_simulation", "extract_trend_estimate", "update_collapsed_rf",
                    "summarise_sanity", "clean_family_name", "combo_data", "combo", "species")
      )
    } else {
      list()
    })

    # Save cache
    saveRDS(rep_results, cache_file)
    message("      Saved: ", basename(cache_file))

    return(rep_results)
  })
})

message("\n=== All fitting complete ===")
message("Fitted results cached in: ", results_dir)

# Save combined results
saveRDS(all_fitted_results, file.path(results_dir, "all-fitted-results.rds"))
message("\nCombined results saved: ", file.path(results_dir, "all-fitted-results.rds"))

# Reset to sequential processing
future::plan(future::sequential)
