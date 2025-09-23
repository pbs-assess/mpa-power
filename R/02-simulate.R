source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

library(purrr)

# Grid for data simulation
# ------------------------
restricted_df <- hbll_grid |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) %>%
  st_join(., comm_ll_activity_status |> st_transform(crs = st_crs(.)), join = st_within) |>
  mutate(activity_status_label = ifelse(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label != "outside", 1, 0)) |>
  st_drop_geometry() |>
  mutate(log_depth = log(depth_m))

# Load allocations for status quo sampling
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
grid_allocations <- left_join(hbll_grid, hbll_allocations)

# Prepare historical data for comparison and future modelling
historical <- sp_dat |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  st_join(comm_ll_activity_status |> st_transform(crs = 3156), join = st_within) |>
  mutate(activity_status_label = if_else(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label == "outside", 0, 1)) |>
  st_join(hbll_grid_poly |> select(block_id, grouping_code) |> st_transform(crs = 3156), join = st_within) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, fishing_event_id, latitude, longitude, X, Y,
    block_id, fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y,
    depth_m, offset, hook_count,
    catch_count, restricted) |>
  mutate(after = 0) |>
  left_join(hbll_allocations,
    by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id")) |>
  mutate(observed = catch_count, replicate = 0, d = "historical")

hbll_last_sampled_year <- historical |>
  group_by(survey_abbrev) |>
  slice(which.max(year)) |>
  select(ssid, survey_abbrev, last_sampled_year = year)


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
  hbll_last_sampled_year <- sp_dat |>
    group_by(survey_abbrev) |>
    summarise(last_sampled_year = max(year))

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

sample_historical_locations <- function(sim_dat, allocation, grouping_vars) {
  hist_locations <- unique(historical$block_id)
  hist_locations <- hist_locations[!is.na(hist_locations)]
  sim_dat |>
    filter(block_id %in% hist_locations) |>
  sample_by_plan(
    sim_dat = _,
    sampling_effort = allocation |> mutate(n_samps = allocation),
    grouping_vars = grouping_vars) |>
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
  # filter out years that we don't sample
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "even")) |>
  mutate(plan = "historical locations only - status quo allocations")
}

status_quo_sampling <- function(sim_dat) {
  sample_by_plan(
    sim_dat = sim_dat,
    sampling_effort = allocation_lu |> mutate(n_samps = allocation),
    grouping_vars = c("survey_abbrev", "year", "spatial_grouping_id", "strata_depth")) |>
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
  # filter out years that we don't sample
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "even")) |>
  mutate(plan = "status quo")
}

# ------------------------------------------------------------
# Run simulations
# ------------------------------------------------------------
# Loaded from 01-fit-conditioning-models.R
surveys <- list(
  list(fit = fit_IN, abbrev = "HBLL INS N", tag_prefix = "ins-n"),
  list(fit = fit_ON, abbrev = "HBLL OUT N", tag_prefix = "out-n"),
  list(fit = fit_OS, abbrev = "HBLL OUT S", tag_prefix = "out-s")
)

# Simple replicates with same parameters
# ------------------------------------------------------------
nreps <- 3
sim_params <- tibble(
  replicate = 1:nreps,
  mpa_trend = 1,
  check_cache = TRUE,
  save_sim = TRUE,
  seed = 42 + (1:nreps) - 1
) |>
mutate(
  formula = map(1:nreps, ~ ~ 1 + restricted * year_covariate),
  year_covariate = map(1:nreps, ~ 1:21),
  rho_V = map(1:nreps, ~ NULL),
  sigma_V = map(1:nreps, ~ NULL),
  phi = map(1:nreps, ~ NULL)
)
sim_dat <- purrr::pmap_dfr(sim_params, run_single_simulation, surveys = surveys) |>
  select(!contains("fyear"))
meep()

test <- historical |>
  group_by(survey_abbrev) |>
  mutate(year_counter = year - min(year) + 1) |>
  drop_na(block_id) |>
  ungroup()
samp <- distinct(test, year_counter, ssid, survey_abbrev, block_id)

# # The old blocks continue to haunt me...
# missing_blocks <- unique(test$block_id)[which(!unique(test$block_id) %in% unique(sim_dat_OS$block_id))]
# missing_blocks
# historical |>
#   filter(block_id %in% missing_blocks)
#
comb <- sim_dat |>
  mutate(year_counter = year - last_sampled_year) |>
  right_join(samp) |>
  bind_rows(test) |>
  drop_na(X)

comb <- sim_dat |>
  sample_status_quo(
    allocation = grid_allocations,
    grouping_vars = c("survey_abbrev", "year", "grouping_code", "replicate")) |>
  bind_rows(test) |>
  drop_na(X)

comb_sf <- comb |> XY_to_sf()

pdat <- comb_sf |>
  filter(survey_abbrev == "HBLL INS N") |>
  filter(year_counter < 21)
hist_years <- unique(filter(pdat, replicate == 0)$year_counter)
pdat |>
  ggplot() +
  geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
  geom_sf(aes(colour = observed)) +
  scale_colour_viridis_c(trans = ggsidekick::fourth_root_power_trans()) +
  gfplot::coord_sf_auto(pdat) +
  facet_wrap(replicate ~ year_counter, ncol = length(hist_years))
meep()

zdf <- comb |>
  group_by(d, survey_abbrev, replicate) |>
  summarise(mean_catch_count = mean(observed),
    zero_count = sum(observed == 0),
    non_zero_count = sum(observed > 0)) |>
  arrange(survey_abbrev, replicate) |>
  mutate(pzero = zero_count / (zero_count + non_zero_count))

ggplot() +
  geom_histogram(data = zdf |> filter(d != "historical"),
    aes(x = mean_catch_count, fill = d), position = "identity", alpha = 0.5) +
  geom_vline(data = zdf |> filter(d == "historical"),
    aes(xintercept = mean_catch_count), color = "black") +
  facet_wrap(~ survey_abbrev, scales = "free")

ggplot() +
  geom_histogram(data = zdf |> filter(d != "historical"),
    aes(x = pzero, fill = d), position = "identity", alpha = 0.5) +
  geom_vline(data = zdf |> filter(d == "historical"),
    aes(xintercept = pzero), color = "black") +
   facet_wrap(~ survey_abbrev, scales = "free")

ggplot(data = comb |> filter(replicate %in% 0:19)) +
  geom_point(mapping = aes(x = year, y = observed, color = d), shape = 21) +
  scale_colour_manual(values = c(historical = "steelblue", simulated = "orange")) +
  facet_wrap(~ survey_abbrev, scales = "free")
  # facet_wrap(~ replicate, scales = "fixed")

# TODO: explain why I am using this metric instead of the raw cumulative catch - rank plots
cdf <- comb |>
  group_by(d, survey_abbrev, replicate, restricted) |>
  arrange(observed) |>
  mutate(
    cumul_catch = cumsum(observed),
    prop_cumul = cumul_catch / sum(observed),  # Proportion of total
    rank_prop = rank(observed) / n()          # Proportion of samples
  )
ggplot(cdf, aes(x = rank_prop, y = prop_cumul, colour = d)) +
  geom_line(aes(group = replicate), alpha = 0.5) +
  scale_colour_manual(values = c(historical = "steelblue", simulated = "orange")) +
  labs(x = "Proportion of samples (ranked)",
        y = "Proportion of total catch") +
  facet_wrap(restricted ~ survey_abbrev)
