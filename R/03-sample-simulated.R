# =============================================================================
# Sample simulated data using different sampling designs
# =============================================================================
# This script loads simulated data and applies various sampling designs to
# test different survey strategies for power analysis.

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

message("Using simulation directory: ", sim_dir)
message("Sample directory: ", sample_dir)

library(purrr)
library(tidyr)

# =============================================================================
# Configuration
# =============================================================================

dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

# Load simulation summary
sim_summary0 <- readRDS(file.path(sim_dir, "simulation-summary.rds"))
sim_summary <- sim_summary0 |>
  mutate(mpa_trend = round(mpa_trend, digits = 4)) #|>
  # left_join(hbll_last_sampled_years, by = "survey_abbrev")

# Load required spatial data
hbll_allocations <- readRDS(hbll_allocations_file) |>
  as_tibble()
# hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
#   filter(stringr::str_detect(survey_abbrev, "HBLL"))


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

historical_locations <- readRDS(historical_locations_file) |>
  tidyr::drop_na(block_id) # needed because there are lat/lon locations that were surveyed
  # but that are not in the simulation grid.

fit_characteristics <- readRDS(fit_characteristics_file) |>
  rename(survey_abbrev = survey) |>
  select(species, survey_abbrev, phi) |>
  distinct()


# =============================================================================
# Helper Functions
# =============================================================================
#' Resolve multiple cache versions for the same simulation combo
#'
#' When legacy cache versions exist for the same parameter combination, keep
#' the version with the most replicates, breaking ties by newest creation time.
resolve_sim_summary_versions <- function(sim_summary) {
  if (nrow(sim_summary) == 0) {
    return(sim_summary)
  }

  if (!"created_date" %in% names(sim_summary)) {
    sim_summary <- sim_summary |>
      mutate(created_date = as.POSIXct(NA))
  } else {
    sim_summary <- sim_summary |>
      mutate(created_date = as.POSIXct(created_date, origin = "1970-01-01", tz = "UTC"))
  }

  resolved <- sim_summary |>
    group_by(species, survey_abbrev, mpa_trend, ar1_scenario, time_scenario) |>
    arrange(desc(n_replicates), desc(created_date), .by_group = TRUE) |>
    slice(1) |>
    ungroup()

  n_dropped <- nrow(sim_summary) - nrow(resolved)
  if (n_dropped > 0) {
    message(
      "Resolved ", n_dropped, " duplicate simulation summary row(s) by keeping ",
      "the cache version with the most replicates."
    )
  }

  resolved
}

#' Load simulated data for a specific parameter combination
#'
#' @param species Species name
#' @param mpa_trend MPA trend value
#' @param ar1_scenario AR1 scenario name
#' @param time_scenario Time scenario name
#' @param sim_summary Simulation summary tibble
#' @param sim_dir Directory containing simulated data
#' @param replicates Integer vector of replicate numbers to load, or NULL to load all
#' @param sim_hash Optional simulation hash to require
#' @param file_pattern Optional exact file pattern from the simulation summary
#'
#' @return Simulated data tibble
load_sim_data <- function(species, survey_abbrev, mpa_trend, ar1_scenario,
                          time_scenario, sim_summary, sim_dir, replicates = NULL,
                          sim_hash = NULL, file_pattern = NULL) {

  # Find matching file
  file_info <- sim_summary |>
    filter(
      species == .env$species,
      survey_abbrev == .env$survey_abbrev,
      mpa_trend == .env$mpa_trend,
      ar1_scenario == .env$ar1_scenario,
      time_scenario == .env$time_scenario
    )

  if (!is.null(sim_hash)) {
    file_info <- file_info |>
      filter(sim_hash == .env$sim_hash)
  }

  if (!is.null(file_pattern)) {
    file_info <- file_info |>
      filter(file == .env$file_pattern)
  }

  if (nrow(file_info) == 0) {
    stop("No simulation found for: ", species, ", survey_abbrev=", survey_abbrev,
         ", mpa_trend=", mpa_trend,
         ", ar1=", ar1_scenario, ", time=", time_scenario)
  }

  if (nrow(file_info) > 1) {
    if ("created_date" %in% names(file_info)) {
      file_info <- file_info |>
        arrange(desc(created_date))
      warning("Multiple simulations found; using the newest created simulation set.")
    } else {
      warning("Multiple simulations found, using first - meaning you haven't accounted for some kind of parameter combo")
    }
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

    # Filter to requested replicates if specified
    if (!is.null(replicates)) {
      rep_patterns <- sprintf("rep%03d", replicates)
      replicate_files <- replicate_files[grepl(paste(rep_patterns, collapse = "|"), replicate_files)]

      if (length(replicate_files) == 0) {
        stop("No files found for requested replicates: ", paste(replicates, collapse = ", "))
      }

      message("Loading ", length(replicate_files), " requested replicate file(s)")
    } else {
      message("Loading ", length(replicate_files), " replicate files matching: ", file_pattern)
    }

    # Load and combine replicate files
    sim_dat <- purrr::map_dfr(replicate_files, readRDS)

    # Verify we have all expected replicates
    n_reps <- length(unique(sim_dat$replicate))

    if (is.null(replicates)) {
      expected_reps <- file_info$n_replicates
      if (n_reps != expected_reps) {
        warning("Expected ", expected_reps, " replicates, found ", n_reps,
                " for ", file_pattern)
      } else {
        message("Successfully loaded all ", n_reps, " replicates")
      }
    } else {
      message("Successfully loaded ", n_reps, " replicate(s)")
    }

  } else {
    # Old approach: single combined file (for backward compatibility)
    fpath <- file.path(sim_dir, file_pattern)

    if (!file.exists(fpath)) {
      stop("Simulation file not found: ", fpath)
    }

    message("Loading combined file: ", basename(fpath))
    sim_dat <- readRDS(fpath)

    # Filter to requested replicates if specified
    if (!is.null(replicates)) {
      sim_dat <- sim_dat |> filter(replicate %in% replicates)
      message("Filtered to ", length(unique(sim_dat$replicate)), " replicate(s)")
    }
  }

  # Add joins that were moved from script 02 to reduce file size
  sim_dat <- sim_dat |>
    left_join(hbll_allocations, by = c("survey_abbrev", "grouping_code"))

  historical_locations <- readRDS(historical_locations_file) |>
    tidyr::drop_na(block_id) |>
    mutate(historical_location = 1) |>
    select(survey_abbrev, block_id, historical_location)

  sim_dat <- sim_dat |>
    mutate(historical_location = ifelse(block_id %in% historical_locations$block_id, 1, 0))

  return(sim_dat)
}

#' Get available replicate IDs for a simulation file pattern
#'
#' @param file_pattern File pattern from simulation summary
#' @param sim_dir Directory containing simulated data
#'
#' @return Sorted integer vector of available replicate IDs
get_available_sim_replicates <- function(file_pattern, sim_dir) {
  if (!grepl("-rep\\*-", file_pattern)) {
    return(integer())
  }

  replicate_files <- list.files(
    path = sim_dir,
    pattern = glob2rx(file_pattern),
    full.names = FALSE
  )

  if (length(replicate_files) == 0) {
    return(integer())
  }

  replicate_ids <- stringr::str_match(replicate_files, "rep([0-9]{3})")[, 2] |>
    as.integer()

  sort(unique(replicate_ids[!is.na(replicate_ids)]))
}

#' Generate clean filename for sampled data
#'
#' @param species Species name
#' @param survey_abbrev Survey abbreviation
#' @param mpa_trend MPA trend value
#' @param ar1_scenario AR1 scenario name
#' @param time_scenario Time scenario name
#' @param plan Sampling plan name
#' @param replicate Replicate number
#'
#' @return Character string with filename
generate_sample_filename <- function(species, survey_abbrev, mpa_trend,
                                     ar1_scenario, time_scenario, plan, replicate) {
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
    plan_slug, "_",
    sprintf("rep%03d", replicate),
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

draw_betabinomial_observed <- function(mu, hook_count, phi) {
  mu <- pmin(pmax(mu, .Machine$double.eps), 1 - .Machine$double.eps)
  p <- stats::rbeta(length(mu), shape1 = mu * phi, shape2 = (1 - mu) * phi)
  stats::rbinom(length(mu), size = hook_count, prob = p)
}

reallocate_bootstrap_locations_outside_mpa <- function(bootstrap_locations, sim_grid,
                                                       reallocation_years) {
  locations_to_keep <- bootstrap_locations |>
    anti_join(reallocation_years, by = c("survey_abbrev", "year")) |>
    mutate(replacement_draw = FALSE)

  bootstrap_locations_outside <- bootstrap_locations |>
    semi_join(reallocation_years, by = c("survey_abbrev", "year")) |>
    filter(restricted == 0) |>
    mutate(replacement_draw = FALSE)

  replacement_locations <- bootstrap_locations |>
    semi_join(reallocation_years, by = c("survey_abbrev", "year")) |>
    group_by(survey_abbrev, year) |>
    summarise(n_restricted = sum(restricted == 1), .groups = "drop") |>
    filter(n_restricted > 0) |>
    purrr::pmap_dfr(function(survey_abbrev, year, n_restricted) {
      available_outside_locations <- bootstrap_locations_outside |>
        filter(survey_abbrev == .env$survey_abbrev, year == .env$year)

      if (nrow(available_outside_locations) == 0) {
        available_outside_locations <- sim_grid |>
          filter(survey_abbrev == .env$survey_abbrev, restricted == 0) |>
          transmute(survey_abbrev, year = .env$year, X, Y, grouping_code, restricted)
      }

      available_outside_locations |>
        slice_sample(n = n_restricted, replace = TRUE) |>
        mutate(replacement_draw = TRUE)
    })

  bind_rows(
    locations_to_keep,
    bootstrap_locations_outside,
    replacement_locations
  )
}

bootstrap_historical_survey_years <- function(sim_dat, species, hist_clean_dir,
                                              seed = NULL, drop_restricted = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  future_survey_years <- filter_hbll_survey_years(sim_dat) |>
    distinct(survey_abbrev, year)

  sim_grid <- sim_dat |>
    distinct(survey_abbrev, X, Y, grouping_code, restricted)

  historical_templates <- purrr::map_dfr(unique(future_survey_years$survey_abbrev), function(survey_abbrev) {
    hist_file <- file.path(
      hist_clean_dir,
      paste0(sp_to_hyphens(species), "-", sp_to_hyphens(survey_abbrev), ".rds")
    )

    if (!file.exists(hist_file)) {
      stop("Historical cleaned data not found: ", hist_file, call. = FALSE)
    }

    readRDS(hist_file) |>
      select(survey_abbrev, hist_year = year, X, Y, restricted)
  })

  mapped_templates <- purrr::map_dfr(unique(historical_templates$survey_abbrev), function(survey_abbrev) {
    purrr::map_dfr(sort(unique(historical_templates$restricted[historical_templates$survey_abbrev == survey_abbrev])), function(restricted_value) {
      hist_subset <- historical_templates |>
        filter(survey_abbrev == .env$survey_abbrev, restricted == .env$restricted_value)

      grid_subset <- sim_grid |>
        filter(survey_abbrev == .env$survey_abbrev, restricted == .env$restricted_value)

      if (nrow(hist_subset) == 0 || nrow(grid_subset) == 0) {
        return(tibble())
      }

      hist_sf <- hist_subset |>
        mutate(x_m = X * 1000, y_m = Y * 1000) |>
        sf::st_as_sf(coords = c("x_m", "y_m"), crs = 32609)

      grid_sf <- grid_subset |>
        mutate(x_m = X * 1000, y_m = Y * 1000) |>
        sf::st_as_sf(coords = c("x_m", "y_m"), crs = 32609)

      nearest_idx <- sf::st_nearest_feature(hist_sf, grid_sf)

      tibble(
        survey_abbrev = hist_subset$survey_abbrev,
        hist_year = hist_subset$hist_year,
        X = grid_subset$X[nearest_idx],
        Y = grid_subset$Y[nearest_idx],
        grouping_code = grid_subset$grouping_code[nearest_idx],
        restricted = grid_subset$restricted[nearest_idx]
      ) |>
        distinct()
    })
  })

  sampled_year_map <- future_survey_years |>
    group_by(survey_abbrev) |>
    mutate(
      hist_year = sample(
        unique(mapped_templates$hist_year[mapped_templates$survey_abbrev == first(survey_abbrev)]),
        size = n(),
        replace = TRUE
      )
    ) |>
    ungroup()

  bootstrap_locations <- purrr::pmap_dfr(sampled_year_map, function(survey_abbrev, year, hist_year) {
    mapped_templates |>
      filter(survey_abbrev == .env$survey_abbrev, hist_year == .env$hist_year) |>
      transmute(survey_abbrev, year = year, X, Y, grouping_code, restricted)
  })

  if (drop_restricted) {
    reallocation_years <- future_survey_years |>
      group_by(survey_abbrev) |>
      mutate(first_year = min(year)) |>
      filter((year - first_year) %% 4 != 0) |>
      ungroup() |>
      select(survey_abbrev, year)

    bootstrap_locations <- reallocate_bootstrap_locations_outside_mpa(
      bootstrap_locations = bootstrap_locations,
      sim_grid = sim_grid,
      reallocation_years = reallocation_years
    )
  } else {
    bootstrap_locations <- bootstrap_locations |>
      mutate(replacement_draw = FALSE)
  }

  sampled_dat <- sim_dat |>
    inner_join(bootstrap_locations, by = c("survey_abbrev", "year", "X", "Y", "grouping_code", "restricted")) |>
    left_join(
      fit_characteristics |>
        filter(species == .env$species) |>
        select(survey_abbrev, phi),
      by = "survey_abbrev"
    )

  replacement_idx <- which(sampled_dat$replacement_draw)
  if (length(replacement_idx) > 0) {
    if (anyNA(sampled_dat$phi[replacement_idx])) {
      stop("Missing phi value for replacement draws in species=", species, call. = FALSE)
    }

    sampled_dat$observed[replacement_idx] <- draw_betabinomial_observed(
      mu = sampled_dat$mu[replacement_idx],
      hook_count = sampled_dat$hook_count[replacement_idx],
      phi = sampled_dat$phi[replacement_idx]
    )
  }

  sampled_dat |>
    select(-replacement_draw, -phi)
}

run_sampling <- function(sim_dat, species) {
  # Verify sim_dat contains single replicate
  rep <- unique(sim_dat$replicate)
  if (length(rep) != 1) {
    stop("sim_dat must contain exactly one replicate, found: ", length(rep))
  }

  # Historical location sampling plan ------------------------
  # This was useful for testing accuracy of simulated data
  # sample_effort_historical <- sim_dat |>
  #   filter(historical_location == 1) |>
  #   mutate(n_samps = allocation) |>
  #   select(survey_series_id, survey_abbrev,
  #          year, X, Y, block_id, grouping_code, pfma, strata_depth,
  #          restricted, allocation, n_samps) |>
  #   filter_hbll_survey_years()

  # sampled_historical <- sample_by_plan(
  #   sim_dat = sim_dat,
  #   sampling_effort = sample_effort_historical,
  #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
  #   seed = rep  # Use replicate number as seed
  # ) |>
  #   mutate(plan = "historical locations only")

  if (RUN_NON_BOOTSTRAP_PLANS) {
    # Case 1: Status quo sampling plan ------------------------
# browser()
    sample_effort_status_quo <- sim_dat |>
      mutate(n_samps = allocation) |>
      select(survey_series_id, survey_abbrev,
        year, X, Y, block_id, grouping_code, grouping_spatial_id, grouping_depth_id, strata_depth,
        restricted, allocation, n_samps) |>
      filter_hbll_survey_years()

    sampled_status_quo <- sample_by_plan(
      sim_dat = sim_dat,
      sampling_effort = sample_effort_status_quo,
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep
    ) |>
      mutate(plan = "status quo", replicate = rep)
    # ggplot(sampled_status_quo, aes(X, Y, colour = factor(restricted))) + geom_point() + facet_wrap(~year)
    # Case 2: MPA sampling every 5 years; Status quo reallocated outside MPAs in off years
    sampled_mpas_5_years <- sample_by_plan(
      sim_dat = sim_dat,
      sampling_effort = sample_effort_status_quo |>
        group_by(survey_abbrev) |>
        mutate(first_year = min(year)) |>
        filter(restricted == 0 | (year - first_year) %% 4 == 0) |>
        select(-first_year) |>
        ungroup(),
      grouping_vars = c("survey_abbrev", "year", "grouping_code"),
      seed = rep + 6000
    ) |>
      mutate(plan = "MPAs every 4 years", replicate = rep)
  } else {
    sampled_status_quo <- tibble()
    sampled_mpas_5_years <- tibble()
  }

  sampled_historical_bootstrap <- bootstrap_historical_survey_years(
    sim_dat = sim_dat,
    species = species,
    hist_clean_dir = cleaned_data_dir,
    seed = rep + 12000
  ) |>
    mutate(plan = "historical survey-year bootstrap", replicate = rep)

  sampled_historical_bootstrap_outside_only <- bootstrap_historical_survey_years(
    sim_dat = sim_dat,
    species = species,
    hist_clean_dir = cleaned_data_dir,
    seed = rep + 13000,
    drop_restricted = TRUE
  ) |>
    mutate(plan = "historical survey-year bootstrap - no MPA every 2nd survey", replicate = rep)

  # Case 3: Status quo + 20% effort ------------------------
  # For low power species, does increasing sampling make a difference to power?
  # sampled_status_quo_1.2 <- sample_by_plan(
  #   sim_dat = sim_dat,
  #   sampling_effort = sample_effort_status_quo |> mutate(n_samps = round(n_samps * 1.2)),
  #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
  #   seed = rep + 2000 # Offset seed for different plan
  # ) |>
  #   mutate(plan = "status quo + 20% effort")

  # # Case X: Status quo - no sampling in MPAs ------------------------
  # # Would show nothing because all MPA values are empty. Could do an extrapolation example
  # sampled_status_quo_0_mpa <- sample_by_plan(
  #   sim_dat = sim_dat,
  #   sampling_effort = sample_effort_status_quo |> filter(restricted == 0), # no sampling in MPAs
  #   grouping_vars = c("survey_abbrev", "year", "grouping_code"),
  #   seed = rep + 5000
  # ) |>
  #   mutate(plan = "status quo - no sampling in MPAs")

  # sampled_status_quo_1.1 <- sample_by_plan(
  #   sim_dat = sim_dat,
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
  #   sim_dat = sim_dat,
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
    sampled_mpas_5_years,
    sampled_historical_bootstrap,
    sampled_historical_bootstrap_outside_only
  )
}

#' Load sampled data files for selected parameter combinations
#'
#' @param species Species name(s)
#' @param survey_abbrev Survey abbreviation(s)
#' @param mpa_trend MPA trend value(s) (rate)
#' @param ar1_scenario AR1 scenario name(s)
#' @param time_scenario Time scenario name(s)
#' @param plan Optional sampling plan name(s)
#' @param replicates Optional integer vector of replicates to load
#' @param sampling_summary Sampling summary tibble
#' @param sample_dir Directory containing sampled data files
#'
#' @return Combined sampled data tibble
load_sampled_data <- function(species, survey_abbrev, mpa_trend, ar1_scenario,
                              time_scenario, plan = NULL, replicates = NULL,
                              sampling_summary, sample_dir) {

  selected_files <- sampling_summary |>
    filter(
      species %in% .env$species,
      survey_abbrev %in% .env$survey_abbrev,
      mpa_trend %in% .env$mpa_trend,
      ar1_scenario %in% .env$ar1_scenario,
      time_scenario %in% .env$time_scenario
    )

  if (!is.null(plan)) {
    selected_files <- selected_files |>
      filter(plan %in% .env$plan)
  }

  if (!is.null(replicates)) {
    selected_files <- selected_files |>
      filter(replicate %in% .env$replicates)
  }

  if (nrow(selected_files) == 0) {
    stop(
      "No sampled files found for species=", species,
      ", survey_abbrev=", survey_abbrev,
      ", mpa_trend=", mpa_trend,
      ", ar1_scenario=", ar1_scenario,
      ", time_scenario=", time_scenario
    )
  }

  message("Loading ", nrow(selected_files), " sampled file(s)")

  selected_files |>
    mutate(file_path = file.path(sample_dir, file)) |>
    pull(file_path) |>
    purrr::map_dfr(readRDS)
}
# stop()
# source(here::here("R", "03-sample-simulated.R"))
# =============================================================================
# Defensive check: test sampling on yelloweye before main execution
# =============================================================================
# Get first available MPA trend for yelloweye rockfish
# ye_mpa_trend <- sim_summary |>
#   filter(species == "yelloweye rockfish",
#          ar1_scenario == "fitted_AR1",
#          time_scenario == "thirty_years") |>
#   pull(mpa_trend) |>
#   unique() |>
#   first()

# if (is.na(ye_mpa_trend)) {
#   message("ℹ Skipping yelloweye defensive checks (no data available)")
# } else {
#   message("Running defensive checks with yelloweye MPA trend: ", ye_mpa_trend)

#   ye_sim <- purrr::map_dfr(c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"), ~{
#     load_sim_data(
#       "yelloweye rockfish", .x, ye_mpa_trend, "fitted_AR1", "thirty_years",
#       sim_summary, sim_dir, replicates = 1
#     )
#   })

#   unique(ye_sim$survey_abbrev)

#   ye_sampled <- run_sampling(ye_sim)

#   ggplot() +
#     aes(X, Y, colour = observed, shape = factor(restricted)) +
#     geom_point(data = ye_sampled |> filter(restricted == 0), shape = 21) +
#     geom_point(data = ye_sampled |> filter(restricted == 1), shape = 19) +
#     scale_colour_viridis_c(trans = "log10") +
#     facet_grid(cols = vars(plan), rows = vars(year))

#   ye_sampled |>
#     filter(plan == "status quo") |>
#     janitor::tabyl(year, restricted, survey_abbrev)

#   ye_sampled |>
#     filter(plan == "MPAs every 4 years") |>
#     janitor::tabyl(year, restricted, survey_abbrev)

#   # Quick checks
#   stopifnot("Odd years: HBLL INS N/OUT N" =
#     all((ye_sampled |> filter(survey_abbrev %in% c("HBLL INS N", "HBLL OUT N")) |> pull(year)) %% 2 == 1))
#   stopifnot("Even years: HBLL OUT S" =
#     all((ye_sampled |> filter(survey_abbrev == "HBLL OUT S") |> pull(year)) %% 2 == 0))
#   stopifnot("MPA 4-yr: restricted every 4 years" =
#     all(diff(unique((ye_sampled |> filter(plan == "MPAs every 4 years", survey_abbrev %in% c("HBLL INS N", "HBLL OUT N"), restricted == 1) |> pull(year)))) %in% c(4, 5)))
#   stopifnot("MPA 4-yr: restricted every 4 years" =
#     all(diff(unique((ye_sampled |> filter(plan == "MPAs every 4 years", survey_abbrev == "HBLL OUT S", restricted == 1) |> pull(year)))) %in% c(4, 5)))
#   message("✓ Sampling checks passed")
# }

# =============================================================================
# Main execution
# =============================================================================

### SETTINGS
source(here::here("R", "sample-fit-config.R"))

# Setup parallel processing
if (USE_PARALLEL) {
  if (is.null(N_WORKERS)) N_WORKERS <- max(1L, floor(future::availableCores() / 2))

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

# =============================================================================
# Task grid filtering (optional - set filters to process subset)
# =============================================================================
# Define which combinations to process. Set to NULL to process all.
#
# Usage examples:
#
# 1. Process single species:
#    FILTER_SPECIES <- "yelloweye rockfish"
#
# 2. Process multiple surveys:
#    FILTER_SURVEY <- c("HBLL OUT N", "HBLL OUT S")
#
# 3. Process specific MPA trends:
#    FILTER_MPA_TREND <- c(1.021, 1.024)
#
# 4. Process first 10 replicates only:
#    FILTER_REPLICATES <- 1:10
#
# 5. Combine filters (e.g., yelloweye, HBLL OUT N, first 5 replicates):
#    FILTER_SPECIES <- "yelloweye rockfish"
#    FILTER_SURVEY <- "HBLL OUT N"
#    FILTER_REPLICATES <- 1:5
#
# Available filter options:
#   FILTER_SPECIES: Character vector of species names
#   FILTER_SURVEY: Character vector of survey abbreviations
#   FILTER_MPA_TREND: Numeric vector of MPA trend values
#   FILTER_AR1_SCENARIO: Character vector of AR1 scenario names
#   FILTER_TIME_SCENARIO: Character vector of time scenario names
#   FILTER_REPLICATES: Integer vector of replicate numbers

# Apply filters to simulation summary
filter_config <- list(
  list(
    values = FILTER_SPECIES,
    column = "species",
    label = "species"
  ),
  list(
    values = FILTER_SURVEY,
    column = "survey_abbrev",
    label = "surveys"
  ),
  list(
    values = FILTER_MPA_TREND,
    column = "mpa_trend",
    label = "MPA trends"
  ),
  list(
    values = FILTER_AR1_SCENARIO,
    column = "ar1_scenario",
    label = "AR1 scenarios"
  ),
  list(
    values = FILTER_TIME_SCENARIO,
    column = "time_scenario",
    label = "time scenarios"
  )
)

sim_summary_filtered <- purrr::reduce(filter_config, function(dat, cfg) {
  if (is.null(cfg$values)) {
    return(dat)
  }

  message("Filtering to ", cfg$label, ": ", paste(cfg$values, collapse = ", "))
  dat |>
    filter(.data[[cfg$column]] %in% cfg$values)
}, .init = sim_summary)

sim_summary_filtered <- resolve_sim_summary_versions(sim_summary_filtered)

# Get unique species from filtered summary
species_list <- unique(sim_summary_filtered$species)
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

  # Get all parameter combinations for this species (from filtered summary)
  sp_sims <- sim_summary_filtered |> filter(species == sp)
  message("  ", nrow(sp_sims), " parameter combinations")

  # Capture replicate filter for use in parallel workers
  # Must be captured here (inside function(sp)) so it's in scope for map_fn
  replicate_filter <- FILTER_REPLICATES

  # Process each parameter combination
  sp_metadata <- map_fn(sp_sims, function(...) {

    row <- list(...)

    # Determine which replicate IDs actually exist on disk for this combo.
    available_replicates <- get_available_sim_replicates(row$file, sim_dir)

    if (length(available_replicates) == 0) {
      warning(
        "No simulation replicate files found for species=", row$species,
        ", survey=", row$survey_abbrev,
        ", mpa=", row$mpa_trend,
        ", ar1=", row$ar1_scenario,
        ", time=", row$time_scenario,
        ". Skipping this parameter combination."
      )
      return(tibble())
    }

    replicates_to_process <- if (!is.null(replicate_filter)) {
      intersect(replicate_filter, available_replicates)
    } else {
      available_replicates
    }

    if (length(replicates_to_process) == 0) {
      message(
        "  No requested replicates available for survey=", row$survey_abbrev,
        ", mpa=", row$mpa_trend,
        ", requested=", paste(replicate_filter, collapse = ", "),
        ", available=", paste(available_replicates, collapse = ", "),
        ". Skipping."
      )
      return(tibble())
    }

    # Define sampling plan names
    plan_names <- c(
      if (RUN_NON_BOOTSTRAP_PLANS) c(
        # "historical locations only",
        "status quo",
        "MPAs every 4 years"
      ),
      "historical survey-year bootstrap",
      "historical survey-year bootstrap - no MPA every 2nd survey"
      # "status quo + 20% effort"
      # "status quo - no sampling in MPAs"
    )

    # Process each replicate (filtered if applicable)
    replicate_metadata <- purrr::map_dfr(replicates_to_process, function(rep_num) {

      # Check if all plan files exist for this replicate
      expected_files <- sapply(plan_names, function(plan) {
        fname <- generate_sample_filename(
          species = sp_clean,
          survey_abbrev = row$survey_abbrev,
          mpa_trend = row$mpa_trend,
          ar1_scenario = row$ar1_scenario,
          time_scenario = row$time_scenario,
          plan = plan,
          replicate = rep_num
        )
        file.path(sp_dir, fname)
      })

      sim_created_date <- if ("created_date" %in% names(row)) as.POSIXct(row$created_date) else NA
      sampled_files_current <- all(file.exists(expected_files)) &&
        (is.na(sim_created_date) || all(file.mtime(expected_files) >= sim_created_date))

      # If all files exist for this replicate and are newer than the source simulations, skip and load metadata
      if (sampled_files_current) {
        message("  Cache hit: survey=", row$survey_abbrev, ", mpa=", row$mpa_trend,
                ", ar1=", row$ar1_scenario, ", time=", row$time_scenario, ", rep=", rep_num)

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
            replicate = rep_num,
            file = file.path(sp_clean, basename(expected_files[i]))
          )
        })

        return(file_metadata)
      }

      if (all(file.exists(expected_files)) && !sampled_files_current) {
        message("  Cache stale: resampling survey=", row$survey_abbrev, ", mpa=", row$mpa_trend,
                ", ar1=", row$ar1_scenario, ", time=", row$time_scenario, ", rep=", rep_num)
      }

      # At least one file missing - proceed with sampling this replicate
      message("  - survey=", row$survey_abbrev, ", mpa=", row$mpa_trend,
              ", ar1=", row$ar1_scenario, ", time=", row$time_scenario, ", rep=", rep_num)

      # Load single replicate
      sim_dat <- load_sim_data(
        species = row$species,
        survey_abbrev = row$survey_abbrev,
        mpa_trend = row$mpa_trend,
        ar1_scenario = row$ar1_scenario,
        time_scenario = row$time_scenario,
        sim_summary = sim_summary,
        sim_dir = sim_dir,
        replicates = rep_num,
        sim_hash = row$sim_hash,
        file_pattern = row$file
      )

      # Apply all sampling designs to this replicate
      sampled <- run_sampling(sim_dat = sim_dat, species = row$species)

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

          # Generate filename
          fname <- generate_sample_filename(
            species = sp_clean,
            survey_abbrev = row$survey_abbrev,
            mpa_trend = row$mpa_trend,
            ar1_scenario = row$ar1_scenario,
            time_scenario = row$time_scenario,
            plan = plan_name,
            replicate = rep_num
          )

          fpath <- file.path(sp_dir, fname)

          # Save this plan × replicate
          saveRDS(plan_data, fpath)
          message("    Saved: ", fname)

          # Return metadata for summary
          tibble(
            species = row$species,
            survey_abbrev = row$survey_abbrev,
            mpa_trend = row$mpa_trend,
            ar1_scenario = row$ar1_scenario,
            time_scenario = row$time_scenario,
            plan = plan_name,
            replicate = rep_num,
            file = file.path(sp_clean, fname)
          )
        })

      return(file_metadata)
    })

    return(replicate_metadata)
  }, .options = if (USE_PARALLEL) furrr::furrr_options(seed = TRUE) else list())
})


# Save sampling summary catalog - merge with existing to preserve prior batches
summary_file <- file.path(sample_dir, "sampling-summary.rds")
existing_summary <- if (file.exists(summary_file)) readRDS(summary_file) else tibble()
combined_summary <- bind_rows(existing_summary, sampling_summary) |>
  distinct(species, survey_abbrev, mpa_trend, ar1_scenario, time_scenario, plan, replicate, .keep_all = TRUE)
saveRDS(combined_summary, summary_file)
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
#   survey_abbrev = "HBLL OUT S",
#   mpa_trend = c(1.021, 1.024),
#   ar1_scenario = "fitted_AR1",
#   time_scenario = "twenty-five_years",
#   plan = "status quo",
#   replicates = 1,
#   sampling_summary = sampling_summary,
#   sample_dir = sample_dir
# )
# glimpse(ye_sample)
# max(ye_sample$replicate)

# Test scenario buildng:

# test <- load_sim_data(
#   species = "yelloweye rockfish",
#   survey_abbrev = "HBLL INS N",
#   mpa_trend = 1.009,
#   ar1_scenario = "fitted_AR1",
#   time_scenario = "twenty-five_years",
#   sim_summary = sim_summary,
#   sim_dir = sim_dir,
#   replicates = 1
# )

# # display_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa-500m.rds"))
# test2 <- run_sampling(test, species = "yelloweye rockfish") #|> XY_to_sf()

# wide_test2 <- test2 |>
#   group_by(plan, year, restricted) |>
#   summarise(n = n(), .groups = "drop") |>
#   mutate(
#     plan_group = ifelse(plan %in% c("historical survey-year bootstrap", "status quo"), "status quo", "every 4 years"),
#     strategy = ifelse(grepl("bootstrap", plan), "bootstrap", "allocation")
#   ) |>
#   summarise(n = sum(n), .by = c(plan_group, year, restricted, strategy)) |>
#   pivot_wider(names_from = strategy, values_from = n)

# wide_test2 |>
#   mutate(diff = bootstrap - allocation) |>
#   arrange(plan_group, restricted, year) |>
#   print(n = Inf)

# # on average realised sampling is about 97% of the bootstrapped allocations
# wide_test2 |>
#   mutate(ratio = bootstrap / allocation) |>
#   pull(ratio) |>
#   mean(na.rm = TRUE)

# glimpse(test2)
# test2 |> filter(year %in% 2025:2030) |>
# ggplot(aes(X, Y, colour = factor(restricted))) +
#   geom_point() +
#   facet_grid(year~plan)

# test2 |>
#   group_by(plan, year, restricted) |>
#   summarise(n = n()) |>
#   mutate(plan_group = ifelse(plan %in% c("historical survey-year bootstrap", "status quo"), "status quo", "every 4 years")) |>
#   mutate(strategy = ifelse(grepl("bootstrap", plan), "boostrap", "allocation")) |>
# ggplot(aes(year, n, colour = factor(strategy))) +
#   geom_jitter(width = 0.1) +
#   facet_wrap(factor(restricted)~plan_group)

# test2 |>
#   ggplot(data = _) +
#   geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
#   geom_sf(data = display_mpa, fill = "grey85", colour = "grey85") +
#   geom_sf(aes(colour = factor(restricted)), shape = 21, size = 1.2) +
#   # scale_colour_manual(name = "Restricted", values = c(`0` = 21, `1` = 19)) +
#   facet_grid(plan ~ year) +
#   gfplot::coord_sf_auto(test2)

# ggsave(here::here("figures", "test2.png"), test2, width = 10, height = 10)

