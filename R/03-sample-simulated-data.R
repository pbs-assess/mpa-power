# =============================================================================
# Sample simulated data using different sampling designs
# =============================================================================
# This script loads simulated data and applies various sampling designs to
# test different survey strategies for power analysis.

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

library(purrr)

# =============================================================================
# Configuration
# =============================================================================

sim_dir <- here::here("data-generated", "sim-data")
sample_dir <- here::here("data-generated", "sampled-data")
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

# Load simulation summary
sim_summary <- readRDS(file.path(sim_dir, "simulation-summary.rds"))

# Load required spatial data
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL"))
grid_allocations <- left_join(hbll_grid, hbll_allocations) |>
  XY_to_sf(crs_to = st_crs(simple_mpa) ) |>
  st_join(simple_mpa, join = st_within) |>
  mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
  st_drop_geometry()

historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
  drop_na(block_id) # needed because there are lat/lon locations that were surveyed
  # but that are not in the simulation grid.


# =============================================================================
# Helper Functions
# =============================================================================
#' Load simulated data for a specific parameter combination
#'
#' @param species Species name
#' @param mpa_trend MPA trend value
#' @param ar1_scenario AR1 scenario name
#' @param time_scenario Time scenario name
#' @param sim_summary Simulation summary tibble
#' @param sim_dir Directory containing simulated data
#'
#' @return Simulated data tibble
load_sim_data <- function(species, mpa_trend, ar1_scenario, time_scenario,
                         sim_summary, sim_dir) {

  # Find matching file
  file_info <- sim_summary |>
    filter(
      species == !!species,
      mpa_trend == !!mpa_trend,
      ar1_scenario == !!ar1_scenario,
      time_scenario == !!time_scenario
    )

  if (nrow(file_info) == 0) {
    stop("No simulation found for: ", species, ", mpa_trend=", mpa_trend,
         ", ar1=", ar1_scenario, ", time=", time_scenario)
  }

  if (nrow(file_info) > 1) {
    warning("Multiple simulations found, using first")
    file_info <- file_info[1, ]
  }

  # Load data
  fpath <- file.path(sim_dir, file_info$file)
  sim_dat <- readRDS(fpath)

  message("Loaded: ", file_info$file)

  return(sim_dat)
}

filter_hbll_survey_years <- function(sim_dat) {
  sim_dat |>
    mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
    filter(
    (survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
    (survey_abbrev == "HBLL OUT S" & odd_even == "even")
  )
}

run_sampling <- function(sim_dat) {

  # Get replicates
  replicates <- unique(sim_dat$replicate)

  # Sample each replicate with its own seed
  purrr::map_dfr(replicates, function(rep) {

    # Filter to this replicate
    sim_rep <- sim_dat |> filter(replicate == rep)

    # Historical location sampling plan ------------------------
    sample_effort_historical <- sim_rep |>
      distinct(survey_abbrev, year) |>
      left_join(grid_allocations, by = "survey_abbrev", relationship = "many-to-many") |>
      filter(paste(survey_abbrev, block_id) %in% paste(historical_locations$survey_abbrev, historical_locations$block_id)) |>
      mutate(n_samps = allocation) |>
      filter_hbll_survey_years()

    sampled_historical <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_historical,
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep  # Use replicate number as seed
    ) |>
      mutate(plan = "historical locations only")

    # Status quo sampling plan ------------------------
    sample_effort_status_quo <- sim_rep |>
      distinct(survey_abbrev, year) |>
      left_join(grid_allocations, by = "survey_abbrev", relationship = "many-to-many") |>
      mutate(n_samps = allocation) |>
      filter_hbll_survey_years()

    sampled_status_quo <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo,
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep
    ) |>
      mutate(plan = "status quo")

    sampled_status_quo_1.1 <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.1)),
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 1000  # Offset seed for different plan
    ) |>
      mutate(plan = "status quo + 10% effort")

    sampled_status_quo_1.2 <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.2)),
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 2000
    ) |>
      mutate(plan = "status quo + 20% effort")

    sampled_status_quo_1.4 <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.4)),
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 3000
    ) |>
      mutate(plan = "status quo + 40% effort")

    # Increased sampling effort every 5 years
    sample_effort_status_quo_5 <- sample_effort_status_quo |>
      mutate(n_samps = case_when(
        survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & year %% 5 == 0 ~ round(n_samps * 1.4),
        survey_abbrev == "HBLL OUT S" & (year - 1) %% 5 == 0 ~ round(n_samps * 1.4),
        TRUE ~ round(n_samps)
      ))

    sampled_status_quo_5 <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo_5,
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 4000
    ) |>
      mutate(plan = "status quo + 40% effort every 5 years")

    # Combine all plans for this replicate
    bind_rows(
      sampled_historical,
      sampled_status_quo,
      sampled_status_quo_1.1,
      sampled_status_quo_1.2,
      sampled_status_quo_1.4,
      sampled_status_quo_5
    ) |>
      mutate(replicate = rep)
  })
}

# =============================================================================
# Main execution
# =============================================================================

USE_PARALLEL <- FALSE  # Set to TRUE for HPC
N_WORKERS <- NULL

# Setup parallel processing
if (USE_PARALLEL) {
  if (is.null(N_WORKERS)) N_WORKERS <- floor(parallel::detectCores() / 2)
  future::plan(future::multisession, workers = N_WORKERS)
  map_fn <- furrr::future_pmap_dfr
  message("Using ", N_WORKERS, " parallel workers")
} else {
  future::plan(future::sequential)
  map_fn <- purrr::pmap_dfr
  message("Using sequential processing")
}

# Get unique species
species_list <- unique(sim_summary$species)
message("\nProcessing ", length(species_list), " species")

# Process each species
purrr::walk(species_list, function(sp, check_cache = FALSE) {

  # Generate output filename
  sp_clean <- sp_to_hyphens(sp)
  fname <- paste0(sp_clean, "-all-sampled.rds")
  fpath <- file.path(sample_dir, fname)

  # Check cache
  if (check_cache) {
    if (file.exists(fpath)) {
      message("\n=== Cache hit: ", fname, " ===")
      return(invisible(NULL))
      }
  }

  message("\n========================================")
  message("Sampling simulated data for species: ", sp)
  message("========================================")

  # Get all parameter combinations for this species
  sp_sims <- sim_summary |> filter(species == sp)
  message("  ", nrow(sp_sims), " parameter combinations")

  # Process each parameter combination
  sp_sampled <- map_fn(sp_sims, function(...) {
    row <- list(...)

    message("  - mpa=", row$mpa_trend, ", ar1=", row$ar1_scenario,
            ", time=", row$time_scenario)

    # Load simulation
    sim_dat <- load_sim_data(
      species = row$species,
      mpa_trend = row$mpa_trend,
      ar1_scenario = row$ar1_scenario,
      time_scenario = row$time_scenario,
      sim_summary = sim_summary,
      sim_dir = sim_dir
    )

    # Apply all sampling designs
    sampled <- run_sampling(sim_dat)

    # Add simulation metadata
    sampled |>
      mutate(
        sim_mpa_trend = row$mpa_trend,
        sim_ar1_scenario = row$ar1_scenario,
        sim_time_scenario = row$time_scenario
      )
  }, .options = if (USE_PARALLEL) furrr::furrr_options(seed = TRUE) else list())

  # Save species file
  saveRDS(sp_sampled, fpath)

  message("  Saved: ", fname)
})

message("\n=== All sampling complete ===")
message("Files saved to: ", sample_dir)

ye_samps <- readRDS(file.path(sample_dir, "yelloweye-rockfish-all-sampled.rds"))
glimpse(ye_samps)
distinct(ye_samps, plan, sim_mpa_trend, sim_ar1_scenario, sim_time_scenario)
