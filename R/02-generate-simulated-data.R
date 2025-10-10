# =============================================================================
# Generate simulated data for power analysis
# =============================================================================
# This script generates simulated data across multiple species and parameter
# combinations. The simulations are cached for reuse in sampling experiments.

source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

library(purrr)

# =============================================================================
# Configuration
# =============================================================================

USE_PARALLEL <- FALSE  # Set to TRUE for HPC
# N_WORKERS <- if (USE_PARALLEL) 40 else 1
N_WORKERS <- NULL

# Output directory
sim_dir <- here::here("data-generated", "sim-data")
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)


# Grid for data simulation
# ------------------------
restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  st_drop_geometry() |>
  mutate(log_depth = log(depth_m))

# Load allocations for status quo sampling
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL"))

#' Create parameter grid for simulations
#'
#' @param mpa_trend Vector of MPA trend values (multiplicative)
#' @param ar1_scenarios Tibble with columns: ar1_scenario, rho_V, sigma_V
#' @param time_scenarios Tibble with columns: time_scenario, year_covariate (list column)
#' @param formula_scenarios Tibble with columns: formula_scenario, formula (list column)
#' @param phi Vector of dispersion parameters (NA = use fitted value)
#' @param nreps Number of replicates per parameter combination
#'
#' @return Tibble with parameter grid
create_sim_param_grid <- function(mpa_trend,
                                  ar1_scenarios,
                                  time_scenarios,
                                  formula_scenarios,
                                  phi = c(NA),
                                  start_seed = 42,
                                  nreps = 20) {

  # Create grid
  param_grid <- expand_grid(
    mpa_trend = mpa_trend,
    ar1_scenario = ar1_scenarios$ar1_scenario,
    time_scenario = time_scenarios$time_scenario,
    formula_scenario = formula_scenarios$formula_scenario,
    phi = phi,
    replicate = 1:nreps
  ) |>
    left_join(ar1_scenarios, by = "ar1_scenario") |>
    left_join(time_scenarios, by = "time_scenario") |>
    left_join(formula_scenarios, by = "formula_scenario") |>
    mutate(seed = start_seed + row_number() - 1)

  message("Created parameter grid with ", nrow(param_grid), " rows")
  message("  - ", length(unique(param_grid$mpa_trend)), " MPA trends")
  message("  - ", length(unique(param_grid$ar1_scenario)), " AR1 scenarios")
  message("  - ", length(unique(param_grid$time_scenario)), " time scenarios")
  message("  - ", length(unique(param_grid$formula_scenario)), " formula scenarios")
  message("  - ", nreps, " replicates per combination")

  return(param_grid)
}

#' Generate file name for simulated data
#'
#' @param species Species name
#' @param param_row Single row from parameter grid
#' @param sim_hash Hash for cache validation
#'
#' @return Character string with file name
generate_sim_filename <- function(species, param_row, sim_hash) {
  sp <- sp_to_hyphens(species)

  # Create descriptive name parts
  mpa_str <- paste0("mpa", param_row$mpa_trend)
  ar1_str <- param_row$ar1_scenario
  time_str <- param_row$time_scenario
  formula_str <- if (param_row$formula_scenario != "standard") {
    param_row$formula_scenario
  } else {
    NULL
  }

  # Combine name parts
  name_parts <- c(sp, mpa_str, ar1_str, time_str, formula_str,
                  substr(sim_hash, 1, 8))
  fname <- paste(name_parts, collapse = "-")
  fname <- gsub("[^a-zA-Z0-9_.-]", "-", fname)

  return(paste0(fname, ".rds"))
}

#' Run simulations for a single species across parameter grid
#'
#' @param sp_name Species name
#' @param param_grid Parameter grid from create_sim_param_grid()
#' @param sim_dir Directory to save simulated data
#' @param check_cache Check for cached simulations
#'
#' @return List with simulation results and metadata
run_species_simulations <- function(sp_name,
                                   param_grid,
                                   sim_dir = "data-generated/sim-data",
                                   check_cache = TRUE) {

  message("\n========================================")
  message("Processing species: ", sp_name)
  message("========================================")

  # Load fits
  fits <- fit_species(sp_name, check_cache = TRUE, silent = TRUE,
                     refit_on_collapse = TRUE)

  # Separate passed vs failed
  fits_passed <- purrr::keep(fits, ~ isTRUE(.x$sanity_check$passed))
  fits_failed <- purrr::keep(fits, ~ isFALSE(.x$sanity_check$passed))

  if (length(fits_failed) > 0) {
    failed_surveys <- purrr::map_chr(fits_failed,
                                     ~ unique(.x$data$survey_abbrev))
    message("Warning: Skipping failed surveys for ", sp_name, ": ",
            paste(failed_surveys, collapse = ", "))
  }

  if (length(fits_passed) == 0) {
    warning("No valid fits for ", sp_name, ". Skipping species.")
    return(NULL)
  }

  message("Valid surveys: ", length(fits_passed))

  # Create survey list from passed fits
  survey_names <- names(fits_passed)
  surveys <- purrr::map(survey_names, ~ {
    list(
      fit = fits_passed[[.x]],
      abbrev = unique(fits_passed[[.x]]$data$survey_abbrev),
      tag_prefix = tolower(gsub("fit_", "", .x))
    )
  })

  # Group parameter grid by unique parameter combinations (excluding replicate)
  param_combos <- param_grid |>
    distinct(mpa_trend, ar1_scenario, time_scenario, formula_scenario,
             phi, rho_V, sigma_V, formula, year_covariate)

  message("Running ", nrow(param_combos), " parameter combinations with ",
          max(param_grid$replicate), " replicates each")

  # Process each parameter combination
  results <- purrr::map(1:nrow(param_combos), function(i) {
    param_combo <- param_combos[i, ]

    message("\n--- Parameter combination ", i, "/", nrow(param_combos), " ---")
    message("  MPA trend: ", param_combo$mpa_trend)
    message("  AR1: ", param_combo$ar1_scenario)
    message("  Time: ", param_combo$time_scenario)

    # Get replicates for this parameter combination
    combo_reps <- param_grid |>
      filter(
        mpa_trend == param_combo$mpa_trend,
        ar1_scenario == param_combo$ar1_scenario,
        time_scenario == param_combo$time_scenario,
        formula_scenario == param_combo$formula_scenario
      )

    # Create hash and check cache BEFORE running simulations
    hash_params <- list(
      species = sp_name,
      mpa_trend = param_combo$mpa_trend,
      ar1_scenario = param_combo$ar1_scenario,
      time_scenario = param_combo$time_scenario,
      formula = param_combo$formula[[1]],
      rho_V = param_combo$rho_V,
      sigma_V = param_combo$sigma_V,
      phi = param_combo$phi,
      nreps = nrow(combo_reps)
    )
    sim_hash <- digest::digest(hash_params)

    # Generate filename
    fname <- generate_sim_filename(sp_name, param_combo, sim_hash)
    fpath <- file.path(sim_dir, fname)

    # Check cache - return early if found
    if (check_cache && file.exists(fpath)) {
      message("  Cache hit: ", fname)
      return(list(
        species = sp_name,
        param_combo = param_combo,
        file = fpath,
        from_cache = TRUE
      ))
    }

    message("  Cache miss: running simulations")

    # Run all replicates for this parameter combination
    sim_dat_all_reps <- purrr::pmap_dfr(combo_reps, function(...) {
      row <- list(...)

      message("  - Running replicate ", row$replicate, " (seed: ", row$seed, ")")

      # Run simulation for each survey
      survey_results <- purrr::map_dfr(surveys, function(survey_config) {
        simulate_hbll(
          fit = survey_config$fit,
          restricted_df = restricted_df,
          sim_dir = "data-generated/sim-dat",
          check_cache = FALSE,  # Don't check individual caches - we cache combined file
          save_sim = FALSE,     # Don't save individual sim files
          formula = row$formula[[1]],
          seed = row$seed,
          year_covariate = row$year_covariate[[1]],
          mpa_trend = log(row$mpa_trend),  # Convert to log scale
          rho_V = if (is.na(row$rho_V)) NULL else row$rho_V,
          sigma_V = if (is.na(row$sigma_V)) NULL else row$sigma_V,
          fixed_spatial_re = TRUE,
          fixed_spatiotemporal_re = FALSE,
          phi = if (is.na(row$phi)) NULL else row$phi,
          tag = paste0(survey_config$tag_prefix, "-rep", row$replicate)
        )
      })

      # Add replicate column and remove fyear columns
      survey_results |>
        select(!contains("fyear")) |>
        mutate(replicate = row$replicate)
    })

    # Post-process combined data: add spatial joins and calendar years
    sim_dat_all_reps <- sim_dat_all_reps |>
      left_join(hbll_grid |> select(X, Y, block_id, grouping_code),
                by = c("X", "Y")) |>
      left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
      mutate(
        year_counter = year,  # Store original simulation year
        year = last_sampled_year + year,  # Convert to calendar year
        d = "simulated"
      ) |>
      left_join(hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
      mutate(spatial_grouping_id = ifelse(pfma %in% c("5A", "4B"), "5A4B", pfma))

    # Add parameter metadata as attributes
    attr(sim_dat_all_reps, "sim_params") <- list(
      species = sp_name,
      mpa_trend = param_combo$mpa_trend,
      ar1_scenario = param_combo$ar1_scenario,
      time_scenario = param_combo$time_scenario,
      formula_scenario = param_combo$formula_scenario,
      rho_V = param_combo$rho_V,
      sigma_V = param_combo$sigma_V,
      phi = param_combo$phi,
      year_covariate = param_combo$year_covariate[[1]],
      formula = param_combo$formula[[1]],
      nreps = nrow(combo_reps),
      created_date = Sys.time()
    )

    # Save combined file
    saveRDS(sim_dat_all_reps, fpath)
    message("  Saved: ", fname)

    return(list(
      species = sp_name,
      param_combo = param_combo,
      file = fpath,
      from_cache = FALSE
    ))
  })

  message("\n========================================")
  message("Completed: ", sp_name)
  message("========================================\n")

  return(results)
}

# =============================================================================
# Main execution
# =============================================================================

# Species list
species_list <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish"
)

# species_list <- "yelloweye rockfish"

# =============================================================================
# Define simulation parameter scenarios
# =============================================================================

# AR1 temporal variation scenarios
# - NA_real_ means no AR1 process (converted to NULL when passed to simulate_hbll)
# - rho_V: AR(1) correlation parameter
# - sigma_V: marginal standard deviation for AR(1) process
ar1_scenarios <- tribble(
  ~ar1_scenario, ~rho_V, ~sigma_V,
  "no_AR1", NA_real_, NA_real_#,           # No temporal AR1 variation
  # "moderate_AR1", 0.5, 0.2                # Moderate temporal autocorrelation
)

# Time scenarios
# - year_covariate must be a list column
# - Example with multiple: tribble(~time_scenario, ~year_covariate, "short", list(1:11), "long", list(1:21))
time_scenarios <- tribble(
  ~time_scenario, ~year_covariate,
  "twenty_years", list(1:21)                      # 21 years of simulated data
)

# Formula scenarios
# - formula must be a list column
# - Example with depth: tribble(~formula_scenario, ~formula, "with_depth", list(~ 1 + restricted * year_covariate + poly(log_depth, 2)))
formula_scenarios <- tribble(
  ~formula_scenario, ~formula,
  "standard", list(~ 1 + restricted * year_covariate)  # MPA × time interaction
)

# Create parameter grid
# This creates all combinations of:
#   - 3 MPA trends (no change, 2%, 5% annual increase)
#   - 2 AR1 scenarios (none, moderate)
#   - 1 time scenario (21 years, from default)
#   - 1 formula scenario (restricted:year_covariate interaction, from default)
#   - 5 replicates per combination
# Total: 3 × 2 × 1 × 1 × 5 = 30 simulations per species
param_grid <- create_sim_param_grid(
  mpa_trend = c(1.01, 1.02, 1.05),           # Multiplicative annual trend in MPAs
  ar1_scenarios = ar1_scenarios,          # Explicitly defined above
  time_scenarios = time_scenarios,        # Explicitly defined above
  formula_scenarios = formula_scenarios,  # Explicitly defined above
  nreps = 50                               # Use 5 for testing, increase to 20+ for production
)

message("\n=== Parameter Grid Summary ===")
message("This will generate ", nrow(param_grid), " simulations per species:")
message("  ", length(unique(param_grid$mpa_trend)), " MPA trends × ",
        length(unique(param_grid$ar1_scenario)), " AR1 scenarios × ",
        length(unique(param_grid$time_scenario)), " time scenarios × ",
        length(unique(param_grid$replicate)), " replicates")
message("  Species to process: ", length(species_list))
message("  Total simulations: ", nrow(param_grid) * length(species_list))

# Setup parallel processing
map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)

# Run simulations for all species
if (USE_PARALLEL) {
  all_results <- map_fn(
    species_list,
    run_species_simulations,
    param_grid = param_grid,
    sim_dir = sim_dir,
    check_cache = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )
} else {
  all_results <- map_fn(
    species_list,
    run_species_simulations,
    param_grid = param_grid,
    sim_dir = sim_dir,
    check_cache = TRUE
  )
}

# Create summary
sim_summary <- purrr::map_dfr(all_results, function(sp_results) {
  if (is.null(sp_results)) return(NULL)
  purrr::map_dfr(sp_results, function(x) {
    tibble(
      species = x$species,
      mpa_trend = x$param_combo$mpa_trend,
      ar1_scenario = x$param_combo$ar1_scenario,
      time_scenario = x$param_combo$time_scenario,
      file = basename(x$file),
      from_cache = x$from_cache
    )
  })
})

# Save summary
saveRDS(sim_summary, file.path(sim_dir, "simulation-summary.rds"))
message("\n=== Simulation Summary ===")
print(sim_summary)
message("\nSummary saved to: ", file.path(sim_dir, "simulation-summary.rds"))

# Reset to sequential processing
future::plan(future::sequential)
