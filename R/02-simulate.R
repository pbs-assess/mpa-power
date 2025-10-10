source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

library(purrr)

# Grid for data simulation
# ------------------------
restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  st_drop_geometry() |>
  mutate(log_depth = log(depth_m))

# Load allocations for status quo sampling
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
grid_allocations <- left_join(hbll_grid, hbll_allocations)
hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL"))

#' Run single simulation replicate
#'
#' @param replicate Replicate number
#' @param formula Simulation formula
#' @param year_covariate Year covariate (1:20)
#' @param rho_V AR(1) correlation parameter
#' @param sigma_V AR(1) marginal standard deviation
#' @param mpa_trend MPA recovery trend (multiplicative)
#' @param phi Dispersion parameter
#' @param seed Random seed
#' @param save_sim Whether to save individual simulation files
#' @param surveys List of survey configurations. Default runs all three surveys.
#'   To run specific surveys, pass a list like:
#'   list(list(fit = fit_OS, abbrev = "HBLL OUT S", tag_prefix = "out-s"))
#'
#' @return Tibble with simulation data including replicate column
run_single_simulation <- function(replicate,
                                 formula = ~ 1 + restricted * year_covariate,
                                 year_covariate = 1:20,
                                 rho_V = NULL,
                                 sigma_V = NULL,
                                 mpa_trend = 1,
                                 phi = NULL,
                                 seed = 42,
                                 weights = NULL,
                                 check_cache = FALSE,
                                 save_sim = FALSE,
                                 surveys = list(
                                   list(fit = fit_IN, abbrev = "HBLL INS N", tag_prefix = "ins-n"),
                                   list(fit = fit_ON, abbrev = "HBLL OUT N", tag_prefix = "out-n"),
                                   list(fit = fit_OS, abbrev = "HBLL OUT S", tag_prefix = "out-s")
                                 ),
                                 ...) {

  message("  - Running replicate ", replicate, " (seed: ", seed, ")")

  # Run simulation for each survey using DRY approach
  sim_results <- purrr::map(surveys, function(survey_config) {

    # Get model parameters for sigma_V default
    b <- get_model_pars(survey_config$fit)
    sigma_V_use <- if (!is.null(sigma_V)) sigma_V else sd(b$estimate[grepl("fyear", b$term)])

    # Run simulation
    simulate_hbll(
      fit = survey_config$fit,
      restricted_df = restricted_df,
      sim_dir = "data-generated/sim-dat",
      check_cache = check_cache,
      save_sim = save_sim,
      formula = formula,
      seed = seed,
      year_covariate = year_covariate,
      mpa_trend = log(mpa_trend),
      rho_V = rho_V,
      sigma_V = sigma_V_use,
      fixed_spatial_re = TRUE,
      fixed_spatiotemporal_re = FALSE,
      phi = phi,
      tag = paste0(survey_config$tag_prefix, "-rep", replicate),
      ...
    ) |>
    mutate(survey_abbrev = survey_config$abbrev)
  })
# TODO: add some of this to the main simulate_hbll function?
  # Get last sampled year for calendar year conversion

  # Combine and process results
  bind_rows(sim_results) |>
    left_join(hbll_grid |> select(X, Y, block_id, grouping_code), by = c("X", "Y")) |>
    select(!contains("fyear")) |>
    left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
    mutate(
      year_counter = year_covariate,
      year = last_sampled_year + year,
      d = "simulated",
      replicate = replicate
    ) |>
    left_join(hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
    mutate(spatial_grouping_id = ifelse(pfma %in% c("5A", "4B"), "5A4B", pfma))
}

# ------------------------------------------------------------
# Run simulations
# ------------------------------------------------------------
# Loaded from 01-fit-conditioning-models.R
sp_name <- "yelloweye rockfish"

run_species_sims <- function(sp_name) {
  sp <- sp_to_hyphens(sp_name)
  fits <- fit_species(sp)

  # fits_passed <- purrr::keep(fits, ~ isTRUE(.x$sanity_check$passed))
  # fits_failed <- purrr::keep(fits, ~ isFALSE(.x$sanity_check$passed)) |>
  #   map(~ paste(unique(.x$data$species_common_name), unique(.x$data$survey_abbrev), sep = ", "))
  surveys <- list(
    list(fit = fits$fit_IN, abbrev = "HBLL INS N", tag_prefix = "ins-n"),
    list(fit = fits$fit_ON, abbrev = "HBLL OUT N", tag_prefix = "out-n"),
    list(fit = fits$fit_OS, abbrev = "HBLL OUT S", tag_prefix = "out-s")
  )

# Simple replicates with same parameters
# ------------------------------------------------------------
nreps <- 20

# Helper function to run simulations for a given MPA trend
run_sims <- function(mpa_trend_value, nreps, surveys) {
  sim_params <- tibble(
    replicate = 1:nreps,
    mpa_trend = mpa_trend_value,
    check_cache = TRUE,
    save_sim = TRUE,
    seed = 42 + (1:nreps) - 1
  ) |>
    mutate(
      formula = purrr::map(1:nreps, ~ ~ 1 + restricted * year_covariate),
      year_covariate = purrr::map(1:nreps, ~ 1:21),
      rho_V = purrr::map(1:nreps, ~ NULL),
      sigma_V = purrr::map(1:nreps, ~ NULL),
      phi = purrr::map(1:nreps, ~ NULL)
    )

  purrr::pmap_dfr(sim_params, run_single_simulation, surveys = surveys) |>
    select(!contains("fyear"), -last_sampled_year) |>
    mutate(mpa_trend = mpa_trend_value)

}
# Run simulations for different MPA trends
sim_dat_no_change <- run_sims(mpa_trend_value = 1, 5, surveys)
sim_dat_1.05 <- run_sims(mpa_trend_value = 1.05, nreps, surveys)
sim_dat_1.02 <- run_sims(mpa_trend_value = 1.02, nreps, surveys)
sim_dat_1.01 <- run_sims(mpa_trend_value = 1.01, nreps, surveys)
sim_dat_1.005 <- run_sims(mpa_trend_value = 1.005, nreps, surveys)
}
meep()

# run_species_sims("yelloweye rockfish")

