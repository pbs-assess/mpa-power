# =============================================================================
# Generate simulated data for power analysis
# =============================================================================
# This script generates simulated data across multiple species and parameter
# combinations. The simulations are cached for reuse in sampling experiments.

source(here::here("R", "00-setup.R"))
# source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# =============================================================================
# Configuration
# =============================================================================
USE_PARALLEL <- FALSE
N_WORKERS <- NULL

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 40 #NULL
}

if (Sys.info()['user'] == "jilliandunic") {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 8
}

# Output directory
sim_dir <- here::here("data-generated", "sim-data")
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

# Load and validate recovery rates
# ---------------------------------
recovery_rates <- readRDS(here::here("data-generated", "recovery-rates-lambda.rds"))
message("Loaded recovery rates for ", length(unique(recovery_rates$species)), " species")

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
generate_sim_filename <- function(species, survey_abbrev, param_row, sim_hash) {
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
  name_parts <- c(sp, survey_abbrev, mpa_str, ar1_str, time_str, formula_str,
                  substr(sim_hash, 1, 8))
  fname <- paste(name_parts, collapse = "-")
  fname <- gsub("[^a-zA-Z0-9_.-]", "-", fname)

  return(paste0(fname, ".rds"))
}

#' Load and prepare survey fits for a species
#'
#' @param sp_name Species name
#'
#' @return List with survey_fits (list of survey configs) or NULL if no valid fits
prepare_species_fits <- function(sp_name) {
  message("Loading fits for species: ", sp_name)

  # Load fits
  fits <- fit_species(sp_name)

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

  message("Valid surveys for ", sp_name, ": ", length(fits_passed))

  # Create survey list from passed fits
  survey_names <- names(fits_passed)
  survey_fits <- purrr::map(survey_names, ~ {
    list(
      species = sp_name,
      fit = fits_passed[[.x]],
      abbrev = unique(fits_passed[[.x]]$data$survey_abbrev),
      tag_prefix = tolower(gsub("fit_", "", .x))
    )
  })

  return(survey_fits)
}

#' Run simulation for species × survey × parameter combination
#'
#' @param sp_name Species name
#' @param survey_config Survey configuration (fit, abbrev, tag_prefix)
#' @param param_combo Single row from parameter combinations
#' @param param_grid Full parameter grid (for getting replicates)
#' @param restricted_df Restricted dataframe for simulation
#' @param hbll_grid HBLL grid for spatial joins
#' @param hbll_last_sampled_year Last sampled year by survey
#' @param hbll_allocations HBLL allocations
#' @param sim_dir Directory to save simulated data
#' @param check_cache Check for cached simulations
#'
#' @return List with simulation results and metadata
run_survey_simulation <- function(sp_name,
                                  survey_config,
                                  param_combo,
                                  param_grid,
                                  restricted_df,
                                  hbll_grid,
                                  hbll_last_sampled_year,
                                  hbll_allocations,
                                  sim_dir,
                                  check_cache = TRUE) {

  survey_abbrev <- survey_config$abbrev

  message("\n[", sp_name, " | ", survey_abbrev, "] ",
          "MPA:", param_combo$mpa_trend, " AR1:", param_combo$ar1_scenario,
          " Time:", param_combo$time_scenario)

  # Get replicates for this parameter combination
  combo_reps <- param_grid |>
    filter(
      mpa_trend == param_combo$mpa_trend,
      ar1_scenario == param_combo$ar1_scenario,
      time_scenario == param_combo$time_scenario,
      formula_scenario == param_combo$formula_scenario
    )

  # Create hash and check cache BEFORE running simulations
  # Use string-based hash for consistency across platforms
  hash_string <- paste(
    sp_name,
    survey_abbrev,
    param_combo$mpa_trend,
    param_combo$ar1_scenario,
    param_combo$time_scenario,
    deparse(param_combo$formula[[1]]),
    param_combo$rho_V,
    param_combo$sigma_V,
    param_combo$phi,
    nrow(combo_reps),
    sep = "|"
  )
  sim_hash <- digest::digest(hash_string, algo = "xxhash64")

  # Generate filename
  fname <- generate_sim_filename(sp_name, survey_abbrev, param_combo, sim_hash)
  fpath <- file.path(sim_dir, fname)

  # Check cache - return early if found
  if (check_cache && file.exists(fpath)) {
    message("  Cache hit: ", fname)
    return(list(
      species = sp_name,
      survey_abbrev = survey_abbrev,
      param_combo = param_combo,
      file = fpath,
      from_cache = TRUE
    ))
  }

  message("  Cache miss: running simulations")

  # Run all replicates for this survey
  sim_dat_all_reps <- purrr::pmap_dfr(combo_reps, function(...) {
    row <- list(...)

    message("  - Replicate ", row$replicate, " (seed: ", row$seed, ")")

    # Run simulation for this survey
    survey_sim <- simulate_hbll(
      fit = survey_config$fit,
      restricted_df = restricted_df,
      sim_dir = sim_dir,
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

    # Add replicate column and remove fyear columns
    survey_sim |>
      select(!contains("fyear")) |>
      mutate(replicate = row$replicate)
  })

  # Post-process survey data: add spatial joins and calendar years
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
    survey_abbrev = survey_abbrev,
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

  # Save survey file
  saveRDS(sim_dat_all_reps, fpath)
  message("  Saved: ", fname)

  return(list(
    species = sp_name,
    survey_abbrev = survey_abbrev,
    param_combo = param_combo,
    file = fpath,
    from_cache = FALSE
  ))
}

# =============================================================================
# Main execution
# =============================================================================

# Define species list
species_list <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish",
  "pacific halibut"
)

# Filter to species with recovery rates
missing_rates <- setdiff(species_list, unique(recovery_rates$species))
if (length(missing_rates) > 0) {
  warning("Skipping species without recovery rates: ", paste(missing_rates, collapse = ", "))
  species_list <- setdiff(species_list, missing_rates)
}

if (length(species_list) == 0) {
  stop("No species with recovery rates available. Cannot proceed.")
}

extra_rates <- setdiff(unique(recovery_rates$species), species_list)
if (length(extra_rates) > 0) {
  message("Note: Recovery rates available but not used for: ",
          paste(extra_rates, collapse = ", "))
}

message("Running simulations for ", length(species_list), " species")

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
  "no_AR1", NA_real_, NA_real_,           # No temporal AR1 variation
  "moderate_AR1", 0.3, 0.2                # Some temporal autocorrelation; 0.2 similar to sd on year effects from conditioning models
)

# Time scenarios
# - year_covariate must be a list column
# - Example with multiple: tribble(~time_scenario, ~year_covariate, "short", list(1:11), "long", list(1:21))
time_scenarios <- tribble(
  ~time_scenario, ~year_covariate,
  "twenty-five_years", list(1:25)                      # 25 years of simulated data
)

# Formula scenarios
# - formula must be a list column
# - Example with depth: tribble(~formula_scenario, ~formula, "with_depth", list(~ 1 + restricted * year_covariate + poly(log_depth, 2)))
formula_scenarios <- tribble(
  ~formula_scenario, ~formula,
  "standard", list(~ 1 + restricted * year_covariate)  # MPA × time interaction
)

nreps <- 50

# Note: Parameter grids are now created per-species using recovery rates
# See task grid creation below for species-specific implementation

# =============================================================================
# Prepare fits and create flattened task grid
# =============================================================================

message("\n=== Loading Species Fits ===")
# Load fits for all species
all_species_fits <- purrr::map(species_list, prepare_species_fits)
names(all_species_fits) <- species_list

# Remove species with no valid fits
all_species_fits <- purrr::compact(all_species_fits)

if (length(all_species_fits) == 0) {
  stop("No valid fits for any species. Stopping.")
}

# Create species-specific task grid
message("\n=== Creating Species-Specific Parameter Grids ===")

task_grid <- purrr::map_dfr(names(all_species_fits), function(sp_name) {
  # Get species-specific recovery rates
  sp_rates <- recovery_rates |>
    filter(species == sp_name) |>
    pull(lambda)

  message("Species: ", sp_name, " - Rates: ", paste(round(sp_rates, 4), collapse = ", "))

  # Create parameter grid for this species
  sp_param_grid <- create_sim_param_grid(
    mpa_trend = sp_rates,
    ar1_scenarios = ar1_scenarios,
    time_scenarios = time_scenarios,
    formula_scenarios = formula_scenarios,
    nreps = nreps
  )

  # Get unique parameter combinations (excluding replicate)
  sp_param_combos <- sp_param_grid |>
    distinct(mpa_trend, ar1_scenario, time_scenario, formula_scenario,
             phi, rho_V, sigma_V, formula, year_covariate)

  # Get survey fits for this species
  survey_fits <- all_species_fits[[sp_name]]

  # Create tasks for each survey × param combination
  purrr::map_dfr(survey_fits, function(survey_config) {
    purrr::pmap_dfr(sp_param_combos, function(...) {
      param_combo <- tibble(...)
      tibble(
        species = sp_name,
        survey_config = list(survey_config),
        param_combo = list(param_combo),
        param_grid = list(sp_param_grid)  # Store full grid for this species
      )
    })
  })
})

message("\n=== Task Grid Summary ===")
message("Total tasks: ", nrow(task_grid))
message("  Species: ", length(unique(task_grid$species)))
message("  Average tasks per species: ", round(nrow(task_grid) / length(unique(task_grid$species)), 1))

# Setup parallel processing
map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)

# Run simulations across all tasks in parallel
message("\n=== Running Simulations ===")
if (USE_PARALLEL) {
  all_results <- furrr::future_pmap(
    task_grid,
    function(species, survey_config, param_combo, param_grid) {
      run_survey_simulation(
        sp_name = species,
        survey_config = survey_config,
        param_combo = param_combo,
        param_grid = param_grid,  # Now species-specific from task_grid
        restricted_df = restricted_df,
        hbll_grid = hbll_grid,
        hbll_last_sampled_year = hbll_last_sampled_year,
        hbll_allocations = hbll_allocations,
        sim_dir = sim_dir,
        check_cache = TRUE
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
} else {
  all_results <- purrr::pmap(
    task_grid,
    function(species, survey_config, param_combo, param_grid) {
      run_survey_simulation(
        sp_name = species,
        survey_config = survey_config,
        param_combo = param_combo,
        param_grid = param_grid,  # Now species-specific from task_grid
        restricted_df = restricted_df,
        hbll_grid = hbll_grid,
        hbll_last_sampled_year = hbll_last_sampled_year,
        hbll_allocations = hbll_allocations,
        sim_dir = sim_dir,
        check_cache = TRUE
      )
    }
  )
}

# Create summary
sim_summary <- purrr::map_dfr(all_results, function(x) {
  tibble(
    species = x$species,
    survey_abbrev = x$survey_abbrev,
    mpa_trend = round(x$param_combo$mpa_trend, digits = 3),
    ar1_scenario = x$param_combo$ar1_scenario,
    time_scenario = x$param_combo$time_scenario,
    file = basename(x$file),
    from_cache = x$from_cache
  )
})

# Save summary
saveRDS(sim_summary, file.path(sim_dir, "simulation-summary.rds"))
message("\n=== Simulation Summary ===")
print(sim_summary)
message("\nSummary saved to: ", file.path(sim_dir, "simulation-summary.rds"))

# Reset to sequential processing
future::plan(future::sequential)


test <- readRDS(file.path(sim_dir, "simulation-summary.rds"))
