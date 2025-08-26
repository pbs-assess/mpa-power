source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

sp_dat_block_id <- sp_dat |>
  select(ssid, survey_abbrev, year, latitude, longitude, grouping_code, hook_count, catch_count) |>
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326, crs_to = 32609) |>
  # XY_to_sf(crs_to = st_crs(hbll_grid_poly)) |>
  st_join(hbll_grid_poly |>
    rename(survey_abbrev_grid = survey_abbrev, grouping_code_grid = grouping_code),
    join = st_within) |>
  st_drop_geometry()

# Grid for data simulation
restricted_df <- hbll_grid |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) %>%
  st_join(., comm_ll_activity_status |> st_transform(crs = st_crs(.)), join = st_within) |>
  mutate(activity_status_label = ifelse(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label != "outside", 1, 0)) |>
  st_drop_geometry() |>
  mutate(log_depth = log(depth_m))

# formula = ~ 0 + as.factor(year) + restricted * year_covariate
update_sim_dat <- function(formula = ~ 0 + as.factor(year) + restricted * year_covariate,
                           rho_V = 0, sigma_V = 0, mpa_trend = 1, phi = 0.3, check_cache = FALSE, seed = 42) {

b_IN <- get_model_pars(fit_IN)
b_ON <- get_model_pars(fit_ON)
b_OS <- get_model_pars(fit_OS)

sim_IN <- simulate_hbll(fit_IN, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = check_cache,
  formula = formula,
  seed = seed,
  year_covariate = 1:20,
  mpa_trend = log(mpa_trend),
  rho_V = rho_V,
  sigma_V = if (!is.null(sigma_V)) sigma_V else sd(b_IN$estimate[grepl("fyear", b_IN$term)]),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  phi = phi,
  tag = "ins-n"
) |>
  mutate(survey_abbrev = "HBLL INS N")

sim_ON <- simulate_hbll(fit_ON, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = check_cache,
  formula = formula,
  seed = seed,
  year_covariate = 1:20,
  mpa_trend = log(mpa_trend),
  rho_V = rho_V,
  sigma_V = if (!is.null(sigma_V)) sigma_V else sd(b_ON$estimate[grepl("fyear", b_ON$term)]),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  phi = phi,
  tag = "out-n"
) |>
  mutate(survey_abbrev = "HBLL OUT N")

sim_OS <- simulate_hbll(fit_OS, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = check_cache,
  formula = formula,
  seed = seed,
  year_covariate = 1:20,
  mpa_trend = log(mpa_trend),
  rho_V = rho_V,
  sigma_V = if (!is.null(sigma_V)) sigma_V else sd(b_OS$estimate[grepl("fyear", b_OS$term)]),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  phi = phi,
  tag = "out-s"
) |>
  mutate(survey_abbrev = "HBLL OUT S")

# Make years calendar year
hbll_last_sampled_year <- sp_dat |>
  group_by(survey_abbrev) |>
  summarise(last_sampled_year = max(year))

# Full simulated dataset (all years/grid cells)
sim_dat0 <- bind_rows(sim_IN, sim_ON, sim_OS) |>
  left_join(hbll_grid |> select(X, Y, block_id, grouping_code)) |>
  select(!contains("as.factor(year)")) |>
  left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
  mutate(year = last_sampled_year + year,
         hook_count = round(exp(offset)),
         observed_capped = pmin(observed, hook_count)) # Question: what to do about observed values > hook count?

# Extract simulation parameters from all surveys
sim_params <- list(
  IN = attr(sim_IN, "simulation_params"),
  ON = attr(sim_ON, "simulation_params"),
  OS = attr(sim_OS, "simulation_params")
)

# Spatial version of simulated data
sim_dat_sf <- XY_to_sf(sim_dat0)

# Return list and assign to global environment)
result_list <- list(
  sim_IN = sim_IN,
  sim_ON = sim_ON,
  sim_OS = sim_OS,
  hbll_allocations = hbll_allocations,
  allocation_lu = allocation_lu,
  sim_dat0 = sim_dat0,
  sim_dat_sf = sim_dat_sf,
  sim_dat = sim_dat,
  sim_params = sim_params
)

if (return_result_list) list2env(result_list, envir = .GlobalEnv)

sim_dat0
}

# -----------------------------------------------------------------------------
# Sample simulated data to be used in monitoring models
# -----------------------------------------------------------------------------
sim_dir <- here::here("data-generated", "simulated-sampling-data")
dir.create(sim_dir, showWarnings = FALSE)

# Status quo sampling effort
# --------------------------
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))

sim_dat <- left_join(sim_dat0, hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
  mutate(spatial_grouping_id = ifelse(pfma %in% c("5A", "4B"), "5A4B", pfma)) # Group 5A and 4B together for sampling purposes

allocation_lu <- sim_dat |>
  distinct(survey_abbrev, spatial_grouping_id, strata_depth, allocation)

# HBLL INS N odd years
# HBLL OUT N odd years
# HBLL OUT S even years
sp_dat |>
  group_by(survey_abbrev) |>
  summarise(last_sampled_year = max(year)) |>
  mutate(odd_even = ifelse(last_sampled_year %% 2 == 0, "even", "odd"))

hbll_year_filter <- function(sim_df) {
  sim_df |>
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "even"))
}


# Status quo sampling, based on historical allocations:
# ------------------------------------------------------------


# Note: status quo for outside sampling assumes perfect sampling - based on real data
# this does not seem to be the case, especially in 2022 and 2024, which had 172 and 174
# blocks sampled respectively, rather than the 198 in the historical allocations.
# samp1 |> group_by(survey_abbrev, year) |>
#   summarise(n = n())
# sp_dat |>
#   group_by(survey_abbrev, year) |>
#   summarise(n = n())
status_quo_sampling <- function() {
samp1 <- sample_by_plan(
  sim_dat = sim_dat,
  sampling_effort = allocation_lu |> mutate(n_samps = allocation),
  grouping_vars = c("survey_abbrev", "year", "spatial_grouping_id", "strata_depth")) |>
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
  # maybe it is simplest to sample every year and just filter out years we don't actually sample?
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
      (survey_abbrev %in% c("HBLL OUT S") & odd_even == "even")) |>
  mutate(plan = "status quo")
# samp1_summary <- samp1 |>
#   group_by(year, restricted) |>
#   summarise(n = n()) |>
#   filter(restricted == 1)
  # tidyr::pivot_wider(names_from = restricted, values_from = n, values_fill = 0)
samp1
}

if (FALSE) { # don't plot right now
  plan_name <- "Status quo"
  p1 <- local(plot_sampling_plan(samp1, plan_name))
  ggsave(here::here("draft-figures", paste0("sim-dat-", plan_name, ".png")), width = 9, height = 18)
}

# Compare simulated and historical data
# ------------------------------------------------------------
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
    by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id"))

update_sim_dat(
  formula = ~ 1 + restricted * year_covariate,
  rho_V = NULL,
  sigma_V = NULL,
  mpa_trend = 1,
  phi = NULL,
  check_cache = FALSE,
  seed = 42
)
samp1 <- status_quo_sampling()

sim_params$IN
sim_params$ON
sim_params$OS

simulated <-
  samp1 |>
  select(ssid = "survey_series_id", survey_abbrev, year, year_covariate, X, Y,
    # depth_m,
    offset, hook_count,
    catch_count = "observed", restricted)

# 'Posterior check' - compare simulated and historical data
test <- bind_rows(
  historical |> mutate(d = "historical"),
  simulated |> mutate(d = "simulated")
  ) |>
  mutate(catch_count = ifelse(catch_count > hook_count, hook_count, catch_count))

test |> group_by(d, survey_abbrev) |>
  summarise(mean_catch_count = mean(catch_count),
    zero_count = sum(catch_count == 0), non_zero_count = sum(catch_count > 0))

test |>
  ggplot(aes(x = year, y = catch_count, colour = pfma)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  # not offset problem
  # geom_point(data = test |> group_by(survey_abbrev, year) |> summarise(n = n(), hooks = mean(hook_count)), aes(x = year, y = hooks), color = "red") +
  facet_wrap(~ survey_abbrev, scale = "free_y")

# Posterior check on AR1 deviations
# ---------------------------------
b <- get_model_pars(fit_IN)
year_effects <- b$estimate[grepl("fyear", b$term)]
last_year_mean <- year_effects[length(year_effects)]
last_year <- as.numeric(gsub("fyear", "", b$term[grepl("fyear", b$term)]))[length(year_effects)]
sigma_V <- sd(year_effects)
rho_V <- 0.2

# Create 50 replicates of AR1 simulation
n_replicates <- 50
sim_ar1_replicates <- purrr::map_dfr(1:n_replicates, function(rep_id) {
  tibble(
    replicate = rep_id,
    year = last_year + 1:20,
    estimate = sim_ar1_deviations(
      rho = rho_V, sigma = sigma_V, years = 1:20
    )
  )
}, .id = "replicate") |>
 mutate(estimate = estimate + last_year_mean)

# Plot historical data with all 50 replicates
historical_data <- b |>
  filter(grepl("fyear", term)) |>
  mutate(year = as.numeric(gsub("fyear", "", term)),
         data_type = "Historical") |>
  select(year, estimate, data_type)

ggplot() +
  geom_line(data = sim_ar1_replicates,
            aes(x = year - 1, y = estimate, group = replicate),
            alpha = 0.15, color = "steelblue") +
  geom_line(data = historical_data, aes(x = year, y = estimate),
            color = "black") +
  geom_vline(xintercept = last_year, linetype = "dotted", alpha = 0.7) +
  labs(title = "Historical Year Effects + 50 AR1 Realizations",
        subtitle = paste0("σ_V = ", round(sigma_V, 3), ", ρ = ", rho_V),
        x = "Year", y = "Log abundance (year effect)")

posterior_check_data <- bind_rows(
  # Historical data as "replicate 0"
  historical_data |>
    mutate(replicate = "Historical Data",
            series_id = "Historical Data"),

  sim_ar1_replicates |>
    filter(replicate %in% sample(1:n_replicates, 10)) |>
    select(year, estimate, replicate) |>
    mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
    mutate(survey_abbrev = "HBLL OUT N") |>
    hbll_year_filter()
)

  # Posterior predictive check plot
posterior_check_data |>
    ggplot(aes(x = year, y = estimate)) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.8) +
    facet_grid(~ replicate, scales = "free_x") +
    labs(title = "Posterior Predictive Check: Historical vs AR1 Simulations",
         x = "Year", y = "Log abundance (year effect)") +
    theme(#strip.text = element_text(size = 9),
          # strip.text.x = element_blank(),
          axis.text = element_text(size = 8))





# Test fit to monitoring data and see if we can recover the trends
# -----------------------------------------------------------------------------
# Function does not currently deal with multiple surveys at once, maybe it does just work but I haven't double checked
fit_monitoring <- function(historical, simulated, .formula, mpa_start_year = 2026, ...) {
  last_sampled_year <- max(historical$year)

  d <- bind_rows(historical, simulated) |>
    mutate(mpa_start_year = mpa_start_year,
          last_sampled_year = last_sampled_year,
          year = ifelse(is.na(year), last_sampled_year + year_covariate, year),
          year_collapsed = ifelse(year < mpa_start_year,
                                    year,
                                    last_sampled_year),
          year_since_mpa = ifelse(year < mpa_start_year, 0, year_covariate),
          after_mpa = ifelse(year_since_mpa > 0, 1, 0))

    mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 10)
    fm <- sdmTMB::sdmTMB(
      formula = .formula,
      data = d,
      mesh = mesh,
      family = nbinom2(link = "log"),
      offset = d$offset,
      ...
    )
}




# Question about offset and observations of more fish than hooks..
# I wouldn't mind better understanding is using the censored NB2 before april.
# mean expectation = exp(eta + offset) = mu
# eta = linear predictor (at location s, time t)
# offset = mean log effort of historical data (incorporates hook adjustment corrected)
# observed = NB2(mu, phi)


# samp1 |>
#   mutate(hook_count = round(exp(offset)),
#          observed_capped = pmin(observed, hook_count)) |> # Cap values at hook count
#   mutate(method = "samp") |>
# ggplot(data = _) +
#   geom_point(aes(x = block_id, y = observed), colour = "black", shape = 19, alpha = 0.6) +
#   geom_point(aes(x = block_id, y = observed_capped), colour = "red", shape = 21, alpha = 0.6) +
#   geom_point(data = sp_dat_block_id |> mutate(method = "raw"), aes(x = block_id, y = catch_count), colour = "blue", shape = 21, alpha = 0.6)  +
#   facet_wrap(~ method)


# # distribution of catch counts
# bind_rows(
#   sim_dat |> mutate(d = "simulated", catch_count = observed_hard_cap),
#   sp_dat_block_id |> mutate(d = "historical", catch_count = catch_count)
# ) |>
#   ggplot(aes(x = catch_count, fill = d)) +
#   geom_density(alpha = 0.7) +
#   facet_wrap(~ d, scales = "free") +
#   labs(title = "Distribution of catch counts by method",
#        subtitle = "Compare simulation methods with historical data")

