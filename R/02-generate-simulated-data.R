# =============================================================================
# Generate simulated data for power analysis
# =============================================================================
# This script generates simulated data across multiple species and parameter
# combinations. The simulations are cached for reuse in sampling experiments.

source(here::here("R", "00-setup.R"))
# source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# Make sure we are using latest recovery rates
# knitr::purl(here::here("R", "recovery-rates-clean.Rmd"), output = here::here("R", "recovery-rates-clean.R"))
# # Source the newly created R script to run all code in the console
# source(here::here("R", "recovery-rates-clean.R"))
#
# =============================================================================
# Configuration
# =============================================================================
USE_PARALLEL <- TRUE
N_WORKERS <- 8L

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 80 #NULL
}

if (Sys.info()['user'] %in% c("jillian", "jilliandunic")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- ifelse(Sys.info()['user'] == "jillian", 10, 5)
}

# Output directory
sim_dir <- here::here("data-generated", "sim-data")
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

# Load and validate recovery rates
# ---------------------------------
#recovery_rates <- readRDS(here::here("data-generated", "recovery-rates-lambda.rds"))
recovery_rates <- tidyr::expand_grid(species = c("canary rockfish", "lingcod", "pacific halibut", "pacific spiny dogfish",
  "silvergray rockfish", "yelloweye rockfish", "quillback rockfish", "north pacific spiny dogfish"),
                              lambda = c(exp(log(c(1.05, 1.10, 1.25, 1.5)) / 25)))
                              # lambda = c(exp(log(c(1.25)) / 25)))


message("Loaded recovery rates for ", length(unique(recovery_rates$species)), " species")

ar1_parameters <- readRDS(here::here("data-generated", "ar1-parameters.rds"))
message("Loaded AR1 parameters for ", nrow(ar1_parameters), " species × survey combinations")

# Grid for data simulation
# ------------------------
restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  st_drop_geometry() |>
  mutate(log_depth = log(depth_m))

# Load allocations for status quo sampling
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
  param_grid <- tidyr::expand_grid(
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
  mpa_str <- paste0("mpa", round(param_row$mpa_trend, digits = 3))
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
#' @param fit_dir Directory containing cached model files
#'
#' @return List with survey_fits (list of survey configs) or NULL if no valid fits
prepare_species_fits <- function(sp_name, fit_dir = here::here("data-generated", "fits")) {
  message("Loading fits for species: ", sp_name)

  sp <- sp_to_hyphens(sp_name)

  # Pattern for each survey's betabinomial models
  patterns <- c(
    fit_ON = paste0(sp, "-HBLL-OUT-N-betabinomial-restricted-on-iid-"),
    fit_OS = paste0(sp, "-HBLL-OUT-S-betabinomial-restricted-on-iid-"),
    fit_IN = paste0(sp, "-HBLL-INS-N-betabinomial-restricted-on-iid-")
  )

  # Find fit files and check sanity
  fit_info <- purrr::map(patterns, function(pattern) {
    files <- list.files(fit_dir, pattern = paste0("^", pattern), full.names = TRUE)

    if (length(files) == 0) return(NULL)

    if (length(files) > 1) {
      warning("Multiple models found for ", pattern, ". Using most recent.")
      files <- files[which.max(file.mtime(files))]
    }

    # Load fit to check sanity (but don't keep in memory)
    fit <- readRDS(files)
    passed <- isTRUE(fit$sanity_check$passed)

    if (passed) {
      list(
        fit_path = files,
        abbrev = unique(fit$data$survey_abbrev)
      )
    } else {
      NULL
    }
  })

  # Remove failed fits
  fit_info <- purrr::compact(fit_info)

  if (length(fit_info) == 0) {
    warning("No valid fits for ", sp_name, ". Skipping species.")
    return(NULL)
  }

  message("Valid surveys for ", sp_name, ": ", length(fit_info))

  # Create survey list with fit paths instead of fit objects
  survey_names <- names(fit_info)
  survey_fits <- purrr::map(survey_names, ~ {
    list(
      species = sp_name,
      fit_path = fit_info[[.x]]$fit_path,  # Store path, not object
      abbrev = fit_info[[.x]]$abbrev,
      tag_prefix = tolower(gsub("fit_", "", .x))
    )
  })

  return(survey_fits)
}

#' Check cache and prepare micro-tasks for missing replicates
#'
#' Pre-computes hashes for all parameter combos, checks which replicate files
#' already exist, and creates micro-tasks only for missing replicates.
#'
#' @param task_grid Task grid with species, survey_config, param_combo, param_grid
#' @param sim_dir Directory containing simulation files
#'
#' @return Tibble with micro-tasks (one row per replicate to simulate)
check_cache_and_prepare_tasks <- function(task_grid, sim_dir) {
  message("Checking cache for ", nrow(task_grid), " parameter combinations...")

  micro_tasks <- purrr::pmap_dfr(task_grid, function(species, survey_config, param_combo, param_grid) {
    survey_abbrev <- survey_config$abbrev
    fit_id <- if (!is.null(survey_config$fit_path)) basename(survey_config$fit_path) else "fit-in-memory"

    # Get replicates for this parameter combination
    combo_reps <- param_grid |>
      filter(
        mpa_trend == param_combo$mpa_trend,
        ar1_scenario == param_combo$ar1_scenario,
        time_scenario == param_combo$time_scenario,
        formula_scenario == param_combo$formula_scenario
      )

    # Calculate hash (must match original formula exactly)
    hash_string <- paste(
      species,
      survey_abbrev,
      param_combo$mpa_trend,
      param_combo$ar1_scenario,
      param_combo$time_scenario,
      deparse(param_combo$formula[[1]]),
      param_combo$rho_V,
      param_combo$sigma_V,
      param_combo$phi,
      fit_id,
      nrow(combo_reps),
      sep = "|"
    )
    sim_hash <- digest::digest(hash_string, algo = "xxhash64")

    # Generate base filename pattern (without rep number)
    sp <- sp_to_hyphens(species)
    mpa_str <- paste0("mpa", round(param_combo$mpa_trend, digits = 3))
    ar1_str <- param_combo$ar1_scenario
    time_str <- param_combo$time_scenario
    formula_str <- if (param_combo$formula_scenario != "standard") {
      param_combo$formula_scenario
    } else {
      NULL
    }
    name_parts <- c(sp, survey_abbrev, mpa_str, ar1_str, time_str, formula_str)
    base_name <- paste(name_parts, collapse = "-")
    base_name <- gsub("[^a-zA-Z0-9_.-]", "-", base_name)

    # Check which replicate files are missing
    missing_reps <- purrr::map_dfr(1:nrow(combo_reps), function(i) {
      rep_num <- combo_reps$replicate[i]
      rep_seed <- combo_reps$seed[i]

      # Generate replicate filename
      rep_fname <- paste0(base_name, "-rep", sprintf("%03d", rep_num), "-",
                          substr(sim_hash, 1, 8), ".rds")
      rep_fpath <- file.path(sim_dir, rep_fname)

      # Only add to micro_tasks if file doesn't exist
      if (!file.exists(rep_fpath)) {
        tibble(
          species = species,
          survey_config = list(survey_config),
          param_combo = list(param_combo),
          replicate = rep_num,
          seed = rep_seed,
          replicate_file_path = rep_fpath,
          sim_hash = sim_hash
        )
      } else {
        NULL
      }
    })

    return(missing_reps)
  })

  if (nrow(micro_tasks) == 0) {
    message("All replicates already cached!")
  } else {
    message("Found ", nrow(micro_tasks), " missing replicates to simulate")
  }

  return(micro_tasks)
}

#' Run single replicate simulation and save immediately
#'
#' Simulates ONE replicate and saves it to individual file for fault tolerance.
#'
#' @param species Species name
#' @param survey_config Survey configuration (fit, abbrev, tag_prefix)
#' @param param_combo Parameter combination (single row with mpa_trend, ar1_scenario, etc.)
#' @param replicate Replicate number
#' @param seed Random seed for this replicate
#' @param replicate_file_path Full path where to save this replicate
#' @param sim_hash Hash for this parameter combination
#' @param restricted_df Restricted dataframe for simulation
#' @param hbll_grid HBLL grid for spatial joins
#' @param hbll_last_sampled_year Last sampled year by survey
#'
#' @return TRUE if successful, FALSE otherwise
run_single_replicate_simulation <- function(species,
                                           survey_config,
                                           param_combo,
                                           replicate,
                                           seed,
                                           replicate_file_path,
                                           sim_hash,
                                           restricted_df,
                                           hbll_grid,
                                           hbll_last_sampled_year) {

  survey_abbrev <- survey_config$abbrev

  # Load fit from disk - each worker gets independent TMB object
  fit <- readRDS(survey_config$fit_path)

  # Run simulation for this single replicate
  survey_sim <- simulate_hbll(
    fit = fit,
    restricted_df = restricted_df,
    sim_dir = dirname(replicate_file_path),
    check_cache = FALSE,
    save_sim = FALSE,
    formula = param_combo$formula[[1]],
    seed = seed,
    year_covariate = param_combo$year_covariate[[1]],
    mpa_trend = log(param_combo$mpa_trend),
    rho_V = if (is.na(param_combo$rho_V)) NULL else param_combo$rho_V,
    sigma_V = if (is.na(param_combo$sigma_V)) NULL else param_combo$sigma_V,
    use_fixed_spatial_field = TRUE,
    sigma_E = NULL,
    phi = if (is.na(param_combo$phi)) NULL else param_combo$phi,
    tag = paste0(survey_config$tag_prefix, "-rep", replicate)
  )

  # Add block_id and convert year to calendar year
  sim_dat <- survey_sim |>
    select(!contains("fyear")) |>
    left_join(hbll_grid |> select(X, Y, block_id, grouping_code), by = c("X", "Y")) |>
    left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
    mutate(
      year_counter = year,
      year = last_sampled_year + year,
      d = "simulated",
      replicate = replicate
    ) |>
    select(-last_sampled_year)

  # Add parameter metadata as attributes
  attr(sim_dat, "sim_params") <- list(
    species = species,
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
    replicate = replicate,
    sim_hash = sim_hash,
    created_date = Sys.time()
  )

  # Save immediately
  saveRDS(sim_dat, replicate_file_path)

  return(TRUE)
}

#' Create simulation summary from replicate files
#'
#' Scans directory for replicate files and builds summary tibble.
#'
#' @param sim_dir Directory containing replicate files
#'
#' @return Tibble with simulation summary
create_summary_from_replicate_files <- function(sim_dir) {
  message("Scanning ", sim_dir, " for replicate files...")

  # Find all replicate files
  replicate_files <- list.files(sim_dir, pattern = "-rep[0-9]{3}-.*\\.rds$", full.names = TRUE)

  if (length(replicate_files) == 0) {
    warning("No replicate files found in ", sim_dir)
    return(tibble())
  }

  message("Found ", length(replicate_files), " replicate files")

  # Parse filenames and load attributes to build summary
  summary <- purrr::map_dfr(replicate_files, function(fpath) {
    # Load attributes without loading full data
    params <- attr(readRDS(fpath), "sim_params")

    tibble(
      species = params$species,
      survey_abbrev = params$survey_abbrev,
      mpa_trend = round(params$mpa_trend, 3),
      ar1_scenario = params$ar1_scenario,
      time_scenario = params$time_scenario,
      replicate = params$replicate,
      file = basename(fpath),
      sim_hash = params$sim_hash,
      created_date = params$created_date
    )
  })

  # Group by combo and count replicates
  summary_by_combo <- summary |>
    group_by(species, survey_abbrev, mpa_trend, ar1_scenario, time_scenario, sim_hash) |>
    summarise(
      n_replicates = n(),
      created_date = max(created_date),
      file_pattern = paste0(
        gsub("-rep[0-9]{3}-", "-rep*-", first(file))
      ),
      .groups = "drop"
    ) |>
    rename(file = file_pattern)

  return(summary_by_combo)
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
                                  sim_dir,
                                  check_cache = TRUE) {

  survey_abbrev <- survey_config$abbrev
  fit_id <- if (!is.null(survey_config$fit_path)) basename(survey_config$fit_path) else "fit-in-memory"

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
    fit_id,
    nrow(combo_reps),
    sep = "|"
  )
  sim_hash <- digest::digest(hash_string, algo = "xxhash64")

  # Generate filename
  fname <- generate_sim_filename(sp_name, survey_abbrev, param_combo, sim_hash)
  fpath <- file.path(sim_dir, fname)

  # Return early if cache found
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
  # Process each replicate individually to reduce memory usage
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
      use_fixed_spatial_field = TRUE,
      sigma_E = NULL,  # NULL = use fitted sigma_E (includes spatiotemporal variation)
      phi = if (is.na(row$phi)) NULL else row$phi,
      tag = paste0(survey_config$tag_prefix, "-rep", row$replicate)
    )

    # Add block_id for spatial joins in sampling script
    survey_sim |>
      select(!contains("fyear")) |>
      left_join(hbll_grid |> select(X, Y, block_id, grouping_code), by = c("X", "Y")) |>
      left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
      mutate(
        year_counter = year,  # Store original simulation year
        year = last_sampled_year + year + 1,  # Convert to calendar year with future trend starting at 0
        d = "simulated",
        replicate = row$replicate
      ) |>
      select(-last_sampled_year)  # Don't need to store this in every row
  })

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
# Defensive check: Test simulation validity
# =============================================================================

sp <- sp_to_hyphens("yelloweye rockfish")
# sp <- sp_to_hyphens("lingcod")\
fit_files <- list.files(here::here("data-generated", "fits"),
                        pattern = paste0("^", sp, "-HBLL-INS-N-betabinomial-restricted-on-iid-"),
                        full.names = TRUE)
if (length(fit_files) > 0) {
  fit <- readRDS(fit_files[1])
  fit$sanity_check$passed

  test_sim <- simulate_hbll(
    fit = readRDS(fit_files[1]), restricted_df = restricted_df, sim_dir = sim_dir,
    check_cache = FALSE, save_sim = FALSE, formula = ~ 1 + restricted * year_covariate,
    seed = 999, year_covariate = 0:4, mpa_trend = log(1.01), use_fixed_spatial_field = TRUE
  )

  catch_prop <- test_sim$observed / test_sim$hook_count

  checks <- c(
    `NaN` = sum(is.nan(test_sim$observed)),
    `Inf` = sum(is.infinite(test_sim$observed)),
    all_zero = all(test_sim$observed == 0, na.rm = TRUE),
    bad_range = any(catch_prop < 0 | catch_prop > 1, na.rm = TRUE)
  )

  if (any(checks > 0)) stop("Simulation check failed: ", paste(names(checks)[checks > 0], collapse = ", "))
  message("✓ Simulation check passed")
}

# =============================================================================
# Main execution
# =============================================================================

# Define species list
sp_list <- c(
  # "yelloweye rockfish",
  # "north pacific spiny dogfish",
  "yelloweye rockfish"
  # "quillback rockfish",
  # "pacific halibut",
  # "canary rockfish",
  # "silvergray rockfish"
)

# Filter to species with recovery rates
missing_rates <- setdiff(sp_list, unique(recovery_rates$species))
if (length(missing_rates) > 0) {
  warning("Skipping species without recovery rates: ", paste(missing_rates, collapse = ", "))
  sp_list <- setdiff(sp_list, missing_rates)
}

if (length(sp_list) == 0) {
  stop("No species with recovery rates available. Cannot proceed.")
}

extra_rates <- setdiff(unique(recovery_rates$species), sp_list)
if (length(extra_rates) > 0) {
  message("Note: Recovery rates available but not used for: ",
          paste(extra_rates, collapse = ", "))
}

message("Running simulations for ", length(sp_list), " species")

# sp <- "yelloweye rockfish"

# =============================================================================
# Define simulation parameter scenarios
# =============================================================================

# AR1 temporal variation scenarios
# Note: rho_V and sigma_V are added per-survey in task grid creation
# - NA_real_ means no AR1 process (converted to NULL when passed to simulate_hbll)
# - rho_V: AR(1) correlation parameter
# - sigma_V: marginal standard deviation for AR(1) process
ar1_scenarios <- tribble(
  ~ar1_scenario,
 # "no_AR1",           # No temporal AR1 variation
  "fitted_AR1"        # Use species-survey-specific values from ar1_parameters.rds
)

# Time scenarios
# - year_covariate must be a list column
# - Example with multiple: tribble(~time_scenario, ~year_covariate, "short", list(0:10), "long", list(0:20))
time_scenarios <- tribble(
  ~time_scenario, ~year_covariate,
  "twenty-five_years", list(0:24)                      # 25 future years with the first year at t = 0
)

# Formula scenarios
# - formula must be a list column
# - Example with depth: tribble(~formula_scenario, ~formula, "with_depth", list(~ 1 + restricted * year_covariate + poly(log_depth, 2)))
formula_scenarios <- tribble(
  ~formula_scenario, ~formula,
  "standard", list(~ 1 + restricted * year_covariate)  # MPA × time interaction
)

nreps <- 120
nreps <- 11

# Note: Parameter grids are now created per-species using recovery rates
# See task grid creation below for species-specific implementation

# =============================================================================
# Prepare fits and create flattened task grid
# =============================================================================

message("\n=== Loading Species Fits ===")
# Load fits for all species
all_species_fits <- purrr::map(sp_list, prepare_species_fits)
names(all_species_fits) <- sp_list

# Remove species with no valid fits
all_species_fits <- purrr::compact(all_species_fits)

if (length(all_species_fits) == 0) {
  stop("No valid fits for any species. Stopping.")
}

# Create species × survey-specific task grid
message("\n=== Creating Species × Survey-Specific Parameter Grids ===")

task_grid <- purrr::map_dfr(names(all_species_fits), function(sp_name) {
  # Get species-specific recovery rates
  sp_rates <- recovery_rates |>
    filter(species == sp_name) |>
    pull(lambda)

  message("Species: ", sp_name, " - Rates: ", paste(round(sp_rates, 4), collapse = ", "))

  # Get survey fits for this species
  survey_fits <- all_species_fits[[sp_name]]

  # Create tasks for each survey × param combination
  purrr::map_dfr(survey_fits, function(survey_config) {
    survey_abbrev <- survey_config$abbrev

    # Get species-survey-specific AR1 parameters
    sp_survey_ar1 <- ar1_parameters |>
      filter(species == sp_name, survey_abbrev == !!survey_abbrev)

    if (nrow(sp_survey_ar1) != 1) {
      warning("Expected 1 AR1 parameter row for ", sp_name, " × ", survey_abbrev,
              ", found ", nrow(sp_survey_ar1))
    }

    # Extract AR1 values safely (NA if missing)
    fitted_rho_V <- if (nrow(sp_survey_ar1) == 1) sp_survey_ar1$rho_V[1] else NA_real_
    fitted_sigma_V <- if (nrow(sp_survey_ar1) == 1) sp_survey_ar1$sigma_V[1] else NA_real_

    # Create AR1 scenarios with survey-specific values
    ar1_scenarios_survey <- ar1_scenarios |>
      mutate(
        rho_V = if_else(ar1_scenario == "no_AR1", NA_real_, fitted_rho_V),
        sigma_V = if_else(ar1_scenario == "no_AR1", NA_real_, fitted_sigma_V)
      )

    # Create parameter grid for this species × survey combination
    sp_survey_param_grid <- create_sim_param_grid(
      mpa_trend = sp_rates,
      ar1_scenarios = ar1_scenarios_survey,  # Survey-specific AR1 values
      time_scenarios = time_scenarios,
      formula_scenarios = formula_scenarios,
      nreps = nreps
    )

    # Get unique parameter combinations (excluding replicate)
    sp_survey_param_combos <- sp_survey_param_grid |>
      distinct(mpa_trend, ar1_scenario, time_scenario, formula_scenario,
               phi, rho_V, sigma_V, formula, year_covariate)

    # Create task rows
    purrr::pmap_dfr(sp_survey_param_combos, function(...) {
      param_combo <- tibble(...)
      tibble(
        species = sp_name,
        survey_config = list(survey_config),
        param_combo = list(param_combo),
        param_grid = list(sp_survey_param_grid)
      )
    })
  })
})

message("\n=== Task Grid Summary ===")
message("Total tasks: ", nrow(task_grid))
message("  Species: ", length(unique(task_grid$species)))
message("  Average tasks per species: ", round(nrow(task_grid) / length(unique(task_grid$species)), 1))

# Setup parallel processing
tictoc::tic("Starting micro-task simulations")
map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)

# Pre-computation: check which replicate files are missing
message("\n=== Checking Cache and Preparing Micro-Tasks ===")
micro_tasks <- check_cache_and_prepare_tasks(task_grid, sim_dir)

# Parallel execution (only if missing replicates exist)
if (nrow(micro_tasks) > 0) {
  message("\n=== Running ", nrow(micro_tasks), " Micro-Tasks ===")
  message("Missing replicates will be simulated and saved individually")


  if (USE_PARALLEL) {
    results <- furrr::future_pmap(
      micro_tasks,
      run_single_replicate_simulation,
      restricted_df = restricted_df,
      hbll_grid = hbll_grid,
      hbll_last_sampled_year = hbll_last_sampled_year,
      .options = furrr::furrr_options(seed = TRUE, globals = TRUE)
    )
  } else {
    results <- purrr::pmap(
      micro_tasks,
      run_single_replicate_simulation,
      restricted_df = restricted_df,
      hbll_grid = hbll_grid,
      hbll_last_sampled_year = hbll_last_sampled_year
    )
  }


  # Count successful saves
  n_success <- sum(unlist(results))
  message("Successfully saved ", n_success, " / ", nrow(micro_tasks), " replicates")

} else {
  message("\n=== All Replicates Cached ===")
  message("No simulations needed!")
}

# Create summary by scanning all replicate files
message("\n=== Creating Simulation Summary ===")
sim_summary <- create_summary_from_replicate_files(sim_dir)
saveRDS(sim_summary, file.path(sim_dir, "simulation-summary.rds"))

message("\n=== Simulation Summary ===")
print(sim_summary)
message("\nSummary saved to: ", file.path(sim_dir, "simulation-summary.rds"))

# Reset to sequential processing
future::plan(future::sequential)
tictoc::toc()
meep()

# test <- readRDS(file.path(sim_dir, "simulation-summary.rds"))

# test2 <- readRDS(file.path(sim_dir, test$file[1]))
# glimpse(test2)
# max(test2$replicate)

# fit <- readRDS(file.path("data-generated", "fits", "yelloweye-rockfish-HBLL-OUT-N-betabinomial-restricted-on-iid-31cebb495143a6d5.rds" ))
# recovery_rates |> filter(species == "yelloweye rockfish")
# ar1_parameters |> filter(species == "yelloweye rockfish", survey_abbrev == "HBLL OUT N")

# test <- simulate_hbll(
#   fit = fit,
#   restricted_df = restricted_df,
#   sim_dir = sim_dir,
#   check_cache = FALSE,
#   save_sim = FALSE,
#   formula = ~ 1 + restricted * year_covariate,
#   seed = 1,
#   year_covariate = 0:24,
#   mpa_trend = log(1.02),
#   rho_V = -0.126,
#   sigma_V = 0.106
# )

# test_w <- sample(fit$data$hook_count, size = nrow(input_dat), replace = TRUE)
#  test <- sdmTMB::simulate_new(
#     formula = formula,
#     data = input_dat,
#     mesh = input_mesh,
#     family = family,
#     time = "year",
#     time_varying = ~ 1,
#     time_varying_type = "ar1",
#     sigma_V = sigma_V,
#     rho_time = rho_V,
#     # rho = ar1_rho, # affects AR1 deviations of the GMRF
#     sigma_E = sigma_E_sim,
#     phi = phi,
#     range = range_val,
#     fixed_re = fixed_re,
#     B = B,
#     # offset = offset,
#     # weights = weights,
#     weights = test_w,
#     seed = seed,
#     silent = FALSE,
#     ...
#   ) |>
#     as_tibble()
