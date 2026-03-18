# =============================================================================
# Sample simulated data using different sampling designs
# =============================================================================
# This script loads simulated data and applies various sampling designs to
# test different survey strategies for power analysis.

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

library(purrr)
library(tidyr)

# =============================================================================
# Configuration
# =============================================================================

sim_dir <- here::here("data-generated", "sim-data")
sample_dir <- here::here("data-generated", "sampled-data")
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

# Not used, but sometimes usefulf or me to see
# hbll_last_sampled_years <- tribble(
#   ~survey_abbrev, ~last_sampled_year, ~start_year,
#   "HBLL OUT N", 2024, 2026,
#   "HBLL OUT S", 2025, 2027,
#   "HBLL INS N", 2024, 2026
# )

# Load simulation summary
sim_summary0 <- readRDS(file.path(sim_dir, "simulation-summary.rds"))
sim_summary <- sim_summary0 |>
  mutate(mpa_trend = round(mpa_trend, digits = 3)) #|>
  # left_join(hbll_last_sampled_years, by = "survey_abbrev")

# Load required spatial data
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL"))


# Not needed because it is included in the simulated data
# if (!file.exists(file.path("data-generated", "grid-allocations.rds"))) {
#   grid_allocations <- left_join(hbll_grid, hbll_allocations) |>
#     XY_to_sf(crs_to = st_crs(simple_mpa) ) |>
#     st_join(simple_mpa, join = st_within) |>
#     mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
#     st_drop_geometry()
#   saveRDS(grid_allocations, file.path("data-generated", "grid-allocations.rds"))
# } else {
#   grid_allocations <- readRDS(file.path("data-generated", "grid-allocations.rds"))
# }

historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
  tidyr::drop_na(block_id) # needed because there are lat/lon locations that were surveyed
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
load_sim_data <- function(species, survey_abbrev, mpa_trend, ar1_scenario,
                          time_scenario, sim_summary, sim_dir) {

  # Find matching file
  file_info <- sim_summary |>
    filter(
      species == .env$species,
      survey_abbrev == .env$survey_abbrev,
      mpa_trend == .env$mpa_trend,
      ar1_scenario == .env$ar1_scenario,
      time_scenario == .env$time_scenario
    )

  if (nrow(file_info) == 0) {
    stop("No simulation found for: ", species, ", survey_abbrev=", survey_abbrev,
         ", mpa_trend=", mpa_trend,
         ", ar1=", ar1_scenario, ", time=", time_scenario)
  }

  if (nrow(file_info) > 1) {
    warning("Multiple simulations found, using first - meaning you haven't accounted for some kind of parameter combo")
    file_info <- file_info[1, ]
  }

  # Load data - handle both old combined files and new replicate files
  file_pattern <- file_info$file

  # Check if this is a pattern for replicate files (contains rep*)
  if (grepl("-rep\\*-", file_pattern)) {
    # New approach: glob for all replicate files
    replicate_files <- list.files(
      path = sim_dir,
      pattern = glob2rx(file_pattern),
      full.names = TRUE
    )

    if (length(replicate_files) == 0) {
      stop("No replicate files found matching pattern: ", file_pattern)
    }

    message("Loading ", length(replicate_files), " replicate files matching: ", file_pattern)

    # Load and combine all replicate files
    sim_dat <- purrr::map_dfr(replicate_files, readRDS)

    # Verify we have all expected replicates
    n_reps <- length(unique(sim_dat$replicate))
    expected_reps <- file_info$n_replicates

    if (n_reps != expected_reps) {
      warning("Expected ", expected_reps, " replicates, found ", n_reps,
              " for ", file_pattern)
    } else {
      message("Successfully loaded all ", n_reps, " replicates")
    }

  } else {
    # Old approach: single combined file (for backward compatibility)
    fpath <- file.path(sim_dir, file_pattern)

    if (!file.exists(fpath)) {
      stop("Simulation file not found: ", fpath)
    }

    message("Loading combined file: ", basename(fpath))
    sim_dat <- readRDS(fpath)
  }

  # Add joins that were moved from script 02 to reduce file size
  sim_dat <- sim_dat |>
    left_join(hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
    mutate(spatial_grouping_id = ifelse(pfma %in% c("5A", "4B"), "5A4B", pfma))

  historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
    tidyr::drop_na(block_id) |>
    mutate(historical_location = 1) |>
    select(survey_abbrev, block_id, historical_location)

  sim_dat <- sim_dat |>
    mutate(historical_location = ifelse(block_id %in% historical_locations$block_id, 1, 0))

  return(sim_dat)
}

#' Generate clean filename for sampled data
#'
#' @param species Species name
#' @param survey_abbrev Survey abbreviation
#' @param mpa_trend MPA trend value
#' @param ar1_scenario AR1 scenario name
#' @param time_scenario Time scenario name
#' @param plan Sampling plan name
#'
#' @return Character string with filename
generate_sample_filename <- function(species, survey_abbrev, mpa_trend,
                                     ar1_scenario, time_scenario, plan) {
  # Clean plan name for filename
  plan_slug <- gsub("[^a-zA-Z0-9]", "-", plan) |>
    gsub("-+", "-", x = _) |>
    tolower()

  # Build filename
  paste0(
    survey_abbrev, "_",
    "mpa", round(mpa_trend, digits = 3), "_",
    ar1_scenario, "_",
    time_scenario, "_",
    plan_slug,
    ".rds"
  ) |>
    gsub(" ", "-", x = _)
}

filter_hbll_survey_years <- function(sim_dat) {
  sim_dat |>
    mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
    filter(
    (survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
    (survey_abbrev == "HBLL OUT S" & odd_even == "even")
  )
}

run_sampling <- function(sim_dat, replicates = NULL) {
  # Get replicates
  replicates <- unique(sim_dat$replicate)

  # Sample each replicate with its own seed
  purrr::map_dfr(replicates, function(rep) {

    # Filter to this replicate
    sim_rep <- sim_dat |> filter(replicate == rep)

    # Historical location sampling plan ------------------------
    # This was useful for testing accuracy of simulated data
    # sample_effort_historical <- sim_rep |>
    #   filter(historical_location == 1) |>
    #   mutate(n_samps = allocation) |>
    #   select(survey_series_id, survey_abbrev,
    #          year, X, Y, block_id, grouping_code, pfma, strata_depth,
    #          restricted, allocation, n_samps) |>
    #   filter_hbll_survey_years()

    # sampled_historical <- sample_by_plan(
    #   sim_dat = sim_rep,
    #   sampling_effort = sample_effort_historical,
    #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    #   seed = rep  # Use replicate number as seed
    # ) |>
    #   mutate(plan = "historical locations only")

    # Case 1: Status quo sampling plan ------------------------
    sample_effort_status_quo <- sim_rep |>
      mutate(n_samps = allocation) |>
      select(survey_series_id, survey_abbrev,
        year, X, Y, block_id, grouping_code, pfma, strata_depth,
        restricted, allocation, n_samps) |>
      filter_hbll_survey_years()

    sampled_status_quo <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo,
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep
    ) |>
      mutate(plan = "status quo")

    # Case 2: MPA sampling every 5 years; Status quo reallocated outside MPAs in off years
    sampled_mpas_5_years <- sample_by_plan(
      sim_dat = sim_rep,
      sampling_effort = sample_effort_status_quo |>
        group_by(survey_abbrev) |>
        mutate(first_year = min(year)) |>
        filter(restricted == 0 | (year - first_year) %% 4 == 0) |>
        select(-first_year) |>
        ungroup(),
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 6000
    ) |>
      mutate(plan = "MPAs at 5 year intervals")

    # Case 3: Status quo + 20% effort ------------------------
    # For low power species, does increasing sampling make a difference to power?
    # sampled_status_quo_1.2 <- sample_by_plan(
    #   sim_dat = sim_rep,
    #   sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.2)),
    #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    #   seed = rep + 2000 # Offset seed for different plan
    # ) |>
    #   mutate(plan = "status quo + 20% effort")

    # # Case X: Status quo - no sampling in MPAs ------------------------
    # # Would show nothing because all MPA values are empty. Could do an extrapolation example
    # sampled_status_quo_0_mpa <- sample_by_plan(
    #   sim_dat = sim_rep,
    #   sampling_effort = sample_effort_status_quo |> filter(restricted == 0), # no sampling in MPAs
    #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    #   seed = rep + 5000
    # ) |>
    #   mutate(plan = "status quo - no sampling in MPAs")

    # sampled_status_quo_1.1 <- sample_by_plan(
    #   sim_dat = sim_rep,
    #   sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.1)),
    #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    #   seed = rep + 1000  # Offset seed for different plan
    # ) |>
    #   mutate(plan = "status quo + 10% effort")


    # Increased sampling effort every 5 years ------------------------
    # sample_effort_status_quo_5 <- sample_effort_status_quo |>
    #   mutate(n_samps = case_when(
    #     survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & year %% 5 == 0 ~ round(n_samps * 1.4),
    #     survey_abbrev == "HBLL OUT S" & (year - 1) %% 5 == 0 ~ round(n_samps * 1.4),
    #     TRUE ~ round(n_samps)
    #   ))

    # sampled_status_quo_5_year <- sample_by_plan(
    #   sim_dat = sim_rep,
    #   sampling_effort = sample_effort_status_quo_5,
    #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    #   seed = rep + 4000
    # ) |>
    #   mutate(plan = "status quo + 40% effort every 5 years")

    # Combine all plans for this replicate
    bind_rows(
      # sampled_historical,
      sampled_status_quo,
      # sampled_status_quo_1.1,
      # sampled_status_quo_1.2,
      # sampled_status_quo_1.4,
      # sampled_status_quo_5_year,
      sampled_mpas_5_years
    ) |>
      mutate(replicate = rep)
  })
}

# =============================================================================
# Defensive check: test sampling on yelloweye before main execution
# =============================================================================
ye_sim <- purrr::map_dfr(c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"), ~{
  load_sim_data("yelloweye rockfish", .x, 1.011, "fitted_AR1", "twenty-five_years", sim_summary, sim_dir) |> filter(replicate == 1)
}) |> filter(replicate == 1)

unique(ye_sim$survey_abbrev)

ye_sampled <- run_sampling(ye_sim)

ggplot() +
  aes(X, Y, colour = observed, shape = factor(restricted)) +
  geom_point(data = ye_sampled |> filter(restricted == 0), shape = 21) +
  geom_point(data = ye_sampled |> filter(restricted == 1), shape = 19) +
  scale_colour_viridis_c(trans = "log10") +
  facet_grid(cols = vars(plan), rows = vars(year))

janitor::tabyl(ye_sampled, year, restricted, survey_abbrev)

# =============================================================================
# Main execution
# =============================================================================

USE_PARALLEL <- TRUE
N_WORKERS <- 6

# Setup parallel processing
if (USE_PARALLEL) {
  if (is.null(N_WORKERS)) N_WORKERS <- floor(future::availableCores() / 2)

  if (Sys.info()['user'] %in% c("dunic", "anderson")) {
    future::plan(future::multicore, workers = N_WORKERS)
    message("Using ", N_WORKERS, " parallel workers (multicore)")
  } else {
    future::plan(future::multisession, workers = N_WORKERS)
    message("Using ", N_WORKERS, " parallel workers (multisession)")
  }
  map_fn <- furrr::future_pmap_dfr
} else {
  future::plan(future::sequential)
  map_fn <- purrr::pmap_dfr
  message("Using sequential processing")
}

# Get unique species
species_list <- unique(sim_summary$species)
message("\nProcessing ", length(species_list), " species")

# Process each species
sampling_summary <- purrr::map_dfr(species_list, function(sp) {

  message("\n========================================")
  message("Sampling simulated data for species: ", sp)
  message("========================================")

  # Create species subdirectory
  sp_clean <- sp_to_hyphens(sp)
  sp_dir <- file.path(sample_dir, sp_clean)
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

  # Get all parameter combinations for this species
  sp_sims <- sim_summary |> filter(species == sp)
  message("  ", nrow(sp_sims), " parameter combinations")

  # Process each parameter combination
  sp_metadata <- map_fn(sp_sims, function(...) {

    row <- list(...)

    # Check if all plan files already exist for this parameter combination
    plan_names <- c(
      # "historical locations only",
      "status quo",
      "MPAs at 5 year intervals"#,
      # "status quo + 20% effort"
      # "status quo - no sampling in MPAs"
    )

    expected_files <- sapply(plan_names, function(plan) {
      fname <- generate_sample_filename(
        species = sp_clean,
        survey_abbrev = row$survey_abbrev,
        mpa_trend = row$mpa_trend,
        ar1_scenario = row$ar1_scenario,
        time_scenario = row$time_scenario,
        plan = plan
      )
      file.path(sp_dir, fname)
    })

    # If all files exist, skip sampling and load metadata
    if (all(file.exists(expected_files))) {
      message("  Cache hit: survey=", row$survey_abbrev, ", mpa=", row$mpa_trend,
              ", ar1=", row$ar1_scenario, ", time=", row$time_scenario)

      # Load metadata from existing files
      file_metadata <- purrr::map_dfr(seq_along(expected_files), function(i) {
        existing_data <- readRDS(expected_files[i])
        tibble(
          species = row$species,
          survey_abbrev = row$survey_abbrev,
          mpa_trend = row$mpa_trend,
          ar1_scenario = row$ar1_scenario,
          time_scenario = row$time_scenario,
          plan = unique(existing_data$plan),
          file = file.path(sp_clean, basename(expected_files[i])),
          n_replicates = length(unique(existing_data$replicate))
        )
      })

      return(file_metadata)
    }

    # At least one file missing - proceed with sampling
    message("  - survey=", row$survey_abbrev, ", mpa=", row$mpa_trend,
            ", ar1=", row$ar1_scenario, ", time=", row$time_scenario)

    # Load simulation
    sim_dat <- load_sim_data(
      species = row$species,
      survey_abbrev = row$survey_abbrev,
      mpa_trend = row$mpa_trend,
      ar1_scenario = row$ar1_scenario,
      time_scenario = row$time_scenario,
      sim_summary = sim_summary,
      sim_dir = sim_dir
    )

    # Apply all sampling designs
    sampled <- run_sampling(sim_dat = sim_dat)

    # Add simulation metadata
    sampled <- sampled |>
      mutate(
        sim_species = row$species,
        sim_survey_abbrev = row$survey_abbrev,
        sim_mpa_trend = row$mpa_trend,
        sim_ar1_scenario = row$ar1_scenario,
        sim_time_scenario = row$time_scenario
      )
    # Split by plan and save separately
    file_metadata <- sampled |>
      group_by(plan) |>
      group_split() |>
      purrr::map_dfr(function(plan_data) {
        plan_name <- unique(plan_data$plan)
        n_reps <- length(unique(plan_data$replicate))

        # Generate filename
        fname <- generate_sample_filename(
          species = sp_clean,
          survey_abbrev = row$survey_abbrev,
          mpa_trend = row$mpa_trend,
          ar1_scenario = row$ar1_scenario,
          time_scenario = row$time_scenario,
          plan = plan_name
        )

        fpath <- file.path(sp_dir, fname)

        # Save this plan
        saveRDS(plan_data, fpath)
        message("    Saved: ", fname, " (", n_reps, " replicates)")

        # Return metadata for summary
        tibble(
          species = row$species,
          survey_abbrev = row$survey_abbrev,
          mpa_trend = row$mpa_trend,
          ar1_scenario = row$ar1_scenario,
          time_scenario = row$time_scenario,
          plan = plan_name,
          file = file.path(sp_clean, fname),
          n_replicates = n_reps
        )
      })

    return(file_metadata)
  }, .options = if (USE_PARALLEL) furrr::furrr_options(seed = TRUE) else list())
})

# Save sampling summary catalog
summary_file <- file.path(sample_dir, "sampling-summary.rds")
saveRDS(sampling_summary, summary_file)
message("\n=== All sampling complete ===")
message("Files saved to: ", sample_dir)
message("Summary saved to: ", summary_file)
message("Total files created: ", nrow(sampling_summary))

# Reset to sequential processing
future::plan(future::sequential)
meep()

# Example: inspect sampling summary
# glimpse(sampling_summary)
# head(sampling_summary)

# # Example: load a specific sampling scenario using the new helper function
# ye_sample <- load_sampled_data(
#   species = "yelloweye rockfish",
#   survey_abbrev = "HBLL OUT N",
#   plan = "status quo",
#   mpa_trend = 1.011,
#   ar1_scenario = "no_AR1",
#   time_scenario = "twenty-five_years",
#   sampling_summary = sampling_summary,
#   sample_dir = sample_dir
# )
# glimpse(ye_sample)
# max(ye_sample$replicate)

# Test scenario buildng:

# test <- load_sim_data(
#   species = "yelloweye rockfish",
#   survey_abbrev = "HBLL OUT N",
#   mpa_trend = 1.011,
#   ar1_scenario = "no_AR1",
#   time_scenario = "twenty_years",
#   sim_summary = sim_summary,
#   sim_dir = sim_dir
# ) |> filter(replicate == 1)

# display_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa-500m.rds"))
# test2 <- run_sampling(test) |> XY_to_sf()

# test2 |>
#   ggplot(data = _) +
#   geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
#   geom_sf(data = display_mpa, fill = "grey85", colour = "grey85") +
#   geom_sf(aes(colour = factor(restricted)), shape = 21, size = 1.2) +
#   # scale_colour_manual(name = "Restricted", values = c(`0` = 21, `1` = 19)) +
#   facet_grid(plan ~ year) +
#   gfplot::coord_sf_auto(test2)

# ggsave(here::here("figures", "test2.png"), test2, width = 10, height = 10)
