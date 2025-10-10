# =============================================================================
# Fit sampled data and analyze power
# =============================================================================
# This script fits models to sampled data and evaluates power to detect
# MPA recovery trends under different scenarios.

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

library(purrr)

# stash historical data here for now
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets
sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) |>
  mutate(
    obs_id = factor(row_number()),
    catch_prop = catch_count / hook_count,
    log_depth = log(depth_m)
  )
# Prepare historical data for comparison and future modelling
historical <- sp_dat |>
  XY_to_sf(crs_to = 32609) |>
  st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
  mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
  st_join(hbll_grid_poly |> select(block_id, grouping_code), join = st_within) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, fishing_event_id, latitude, longitude, X, Y,
    block_id, fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y,
    depth_m, offset, hook_count,
    catch_count, restricted) |>
  mutate(after = 0) |>
  left_join(hbll_allocations,
    by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id")) |>
  mutate(observed = catch_count, replicate = 0, d = "historical") |>
  group_by(survey_abbrev) |>
  mutate(year_counter = year - min(year) + 1) |>
  ungroup()

# =============================================================================
# Configuration
# =============================================================================

sample_dir <- here::here("data-generated", "sampled-data")
results_dir <- here::here("data-generated", "power-results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- FALSE
#N_WORKERS <- if (USE_PARALLEL) 20 else 1
N_WORKERS <- NULL

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

    # Fit each replicate
    rep_results <- purrr::map_dfr(replicates, function(rep) {
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

# #' Process one replicate: fit model and extract results
# #'
# #' @param sampled_dat Sampled data for one replicate
# #' @param true_mpa_trend True MPA trend (multiplicative)
# #' @param replicate Replicate number
# #' @param silent Suppress messages
# #'
# #' @return Tibble with results
# process_replicate <- function(sampled_dat, true_mpa_trend, replicate, silent = TRUE) {

#   if (!silent) message("  Processing replicate ", replicate)

#   # Fit model
#   fit <- fit_sampled_model(sampled_dat, silent = silent)

#   # Extract trend
#   trend <- extract_trend_estimate(fit)

#   # Evaluate detection
#   detection <- evaluate_trend_detection(trend, true_mpa_trend)

#   # Compile results
#   tibble(
#     replicate = replicate,
#     estimate = trend$estimate,
#     se = trend$se,
#     ci_lower = trend$ci_lower,
#     ci_upper = trend$ci_upper,
#     converged = trend$converged,
#     error_msg = trend$error_msg,
#     detected = detection$detected,
#     ci_contains_true = detection$ci_contains_true,
#     within_threshold = detection$within_threshold,
#     bias = detection$bias,
#     relative_bias = detection$relative_bias
#   )
# }

# #' Analyze power for one species-parameter-design combination
# #'
# #' @param sampled_result Result object from 03-sample-simulated-data.R
# #' @param use_parallel Use parallel processing for replicates
# #' @param n_workers Number of parallel workers
# #'
# #' @return Tibble with power analysis results
# analyze_power <- function(sampled_result, use_parallel = FALSE, n_workers = NULL) {

#   species <- sampled_result$species
#   params <- sampled_result$params
#   sampled_data <- sampled_result$sampled_data

#   message("\n=== Analyzing power: ", species, " ===")
#   message("  True MPA trend: ", params$mpa_trend)
#   message("  AR1 scenario: ", params$ar1_scenario)
#   message("  Time scenario: ", params$time_scenario)

#   # Setup parallel processing if requested
#   if (use_parallel) {
#     if (is.null(n_workers)) n_workers <- floor(parallel::detectCores() / 2)
#     future::plan(future::multisession, workers = n_workers)
#     map_fn <- furrr::future_map_dfr
#     message("  Using ", n_workers, " parallel workers")
#   } else {
#     future::plan(future::sequential)
#     map_fn <- purrr::map_dfr
#   }

#   # Process each design separately
#   designs <- unique(sampled_data$design_id)

#   results_by_design <- purrr::map_dfr(designs, function(design) {
#     message("  Design: ", design)

#     design_data <- sampled_data |> filter(design_id == design)
#     replicates <- unique(design_data$replicate)

#     # Process all replicates
#     rep_results <- map_fn(replicates, function(rep) {
#       rep_data <- design_data |> filter(replicate == rep)
#       process_replicate(
#         sampled_dat = rep_data,
#         true_mpa_trend = params$mpa_trend,
#         replicate = rep,
#         silent = TRUE
#       )
#     })

#     # Add design info
#     rep_results |>
#       mutate(
#         design_id = design,
#         sampling_freq = unique(design_data$sampling_freq),
#         spatial_design = unique(design_data$spatial_design)
#       )
#   })

#   # Calculate power metrics by design
#   power_summary <- results_by_design |>
#     group_by(design_id, sampling_freq, spatial_design) |>
#     summarise(
#       n_replicates = n(),
#       n_converged = sum(converged, na.rm = TRUE),
#       power = mean(detected, na.rm = TRUE),
#       coverage = mean(ci_contains_true, na.rm = TRUE),
#       within_threshold = mean(within_threshold, na.rm = TRUE),
#       mean_bias = mean(bias, na.rm = TRUE),
#       median_bias = median(bias, na.rm = TRUE),
#       mean_relative_bias = mean(relative_bias, na.rm = TRUE),
#       mean_estimate = mean(estimate, na.rm = TRUE),
#       sd_estimate = sd(estimate, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     mutate(
#       species = species,
#       true_mpa_trend = params$mpa_trend,
#       ar1_scenario = params$ar1_scenario,
#       time_scenario = params$time_scenario
#     )

#   return(list(
#     summary = power_summary,
#     detailed = results_by_design
#   ))
# }

# # =============================================================================
# # Main execution
# # =============================================================================

# # Get list of sampled data files
# sampled_files <- list.files(sample_dir, pattern = "-sampled\\.rds$", full.names = TRUE)

# message("Found ", length(sampled_files), " sampled data files")

# # For testing, process first file
# # In production, iterate over all files
# test_file <- sampled_files[1]
# message("\nProcessing test file: ", basename(test_file))

# sampled_result <- readRDS(test_file)

# # Analyze power
# power_results <- analyze_power(
#   sampled_result = sampled_result,
#   use_parallel = USE_PARALLEL,
#   n_workers = N_WORKERS
# )

# # Display results
# message("\n=== Power Analysis Summary ===")
# print(power_results$summary)

# # Save results
# output_base <- tools::file_path_sans_ext(basename(test_file))
# output_summary <- file.path(results_dir, paste0(output_base, "-power-summary.rds"))
# output_detailed <- file.path(results_dir, paste0(output_base, "-power-detailed.rds"))

# saveRDS(power_results$summary, output_summary)
# saveRDS(power_results$detailed, output_detailed)

# message("\nResults saved:")
# message("  Summary: ", output_summary)
# message("  Detailed: ", output_detailed)

# =============================================================================
# Process all files (commented out for testing)
# =============================================================================

# all_power_results <- purrr::map(sampled_files, function(file) {
#   message("\nProcessing: ", basename(file))
#   sampled_result <- readRDS(file)
#
#   power_results <- analyze_power(
#     sampled_result = sampled_result,
#     use_parallel = USE_PARALLEL,
#     n_workers = N_WORKERS
#   )
#
#   # Save individual results
#   output_base <- tools::file_path_sans_ext(basename(file))
#   saveRDS(power_results$summary,
#           file.path(results_dir, paste0(output_base, "-power-summary.rds")))
#   saveRDS(power_results$detailed,
#           file.path(results_dir, paste0(output_base, "-power-detailed.rds")))
#
#   return(power_results$summary)
# })
#
# # Combine all summaries
# combined_summary <- bind_rows(all_power_results)
# saveRDS(combined_summary, file.path(results_dir, "all-power-summary.rds"))
#
# message("\n=== All analyses complete ===")
# message("Combined results saved to: ", file.path(results_dir, "all-power-summary.rds"))
