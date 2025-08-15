library(sdmTMB)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(purrr)
library(furrr)
library(progress)
theme_set(gfplot::theme_pbs())


# From: https://github.com/pbs-assess/gfmpa/blob/9429210b9da5b5044f3afddcf6eaa9cffaec4d40/analysis/sim.Rmd
# - Simulate recovery in restricted areas to assess whether current survey effort is sufficient to detect population recovery

# - Approach:
#     1. fit geostatistical models to observed data
#     2. use parameters from that model to simulate new data with observations at the actual historically observed locations
#     3. when simulating, simulate recovery at some rate within closed areas and a stationary abundance/density outside of closed areas

# - Dimensions that will likely affect the answer:
#     1. species (therefore estimated spatial and spatiotemporal SD, spatial correlation range, observation error)
#     2. rate of 'recovery'
#     3. number of years observed
#     4. whether one assesses all restricted areas together or individually
#     5. **level of fishing (or activity of concern) that occurred before management actions taken**
#     6. size of restricted area?

# Notes:
# - to start no depth because I think there is a lot of uncertainty in the grid depth
#   - TODO: add depth to prediction grid (part of gfdata updates I have going)

# Housekeeping
# ------------------------------------------------------------
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-utils.R"))

fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S")
  # depth not currently used in sim - I think I need to fix the depth data before including
  # filter(depth_m > 0) |>
  # mutate(depth_mean = mean(log(depth_m), na.rm = TRUE),
  #        depth_sd = sd(log(depth_m), na.rm = TRUE),
  #        depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1],
        #  depth_scaled2 = depth_scaled^2)
hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S")

sp <- sp_to_hyphens("yelloweye rockfish")
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets
comm_ll_activity_status <- readRDS(here::here("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))

mpa_shape_simplified <- comm_ll_activity_status |> st_simplify(dTolerance = 100)

sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) #|>
  # mutate(x = X * 1000, y = Y * 1000) |>
  # st_as_sf(coords = c("x", "y"), crs = 3156) |>
  # st_join(comm_ll_activity_status |> st_transform(crs = 3156), join = st_within) |>
  # mutate(activity_status_label = if_else(is.na(activity_status_label), "outside", activity_status_label)) |>
  # mutate(restricted = ifelse(activity_status_label %in% c("not allowed", "conditional"), 1, 0)) |>
  # st_drop_geometry()

combined <- st_intersection(
    st_as_sfc(st_bbox(hbll_grid_poly |> st_transform(crs = st_crs(mpa_shape_simplified)))),
    st_as_sfc(st_bbox(mpa_shape_simplified))
  ) |>
  st_as_sf()

plot_limits_combined <- get_plot_limits(combined, buffer = 1000)

# Testing
# # fit conditioning model
# fit_ins <- fit_hbll(dat = sp_dat, survey_type = "HBLL INS",
#                     species = "yelloweye-rockfish",
#                     fit_dir = fit_dir)
# fit_out <- fit_hbll(dat = sp_dat, survey_type = "HBLL OUT",
#                     species = "yelloweye-rockfish",
#                     fit_dir = fit_dir)

# # get predictions
# pred_ins <- predict_hbll(fit_ins, hbll_grid, re_form = NULL)
# pred_out <- predict_hbll(fit_out, hbll_grid, re_form = NULL)

# # plot predictions
# plot_hbll_predictions(pred_ins |> filter(year %in% 2024),
#   rotation = 90) + facet_wrap(~ year)
# plot_hbll_predictions(pred_out |> filter(year %in% 2022), rotation = NULL) + facet_wrap(~ year)

# pred_p <- bind_rows(pred_ins, pred_out) |>
#   filter(year == 2024) |>
#   plot_hbll_predictions(rotation = NULL, crs = 4326, buffer = 0) +
#   ggtitle("Yelloweye rockfish HBLL predictions 2024")
# pred_p
# ggsave(here::here("draft-figures", "yelloweye-rockfish-hbll-predictions.pdf"), width = 10, height = 10)

# run simulation
# ------------------------------------------------------------
# start with just simulating coastwide data into the future
# - no depth effect
# - fit each survey region separately and simulate from that in the same temporal
# pattern as sampling occurs

# restricted_df <- sp_dat_sf |> select(fishing_event_id, year, restricted)

# TODO: make it easier to keep track of what should be restricted and what should be outside.
# e.g., to allow quick switch to using national or existing as full protection status.
restricted_df <- hbll_grid |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) %>%
  st_join(., comm_ll_activity_status |> st_transform(crs = st_crs(.)), join = st_within) |>
  mutate(activity_status_label = ifelse(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label != "outside", 1, 0)) |>
  # mutate(restricted = ifelse(activity_status_label %in% c("not allowed", "conditional"), 1, 0)) |>
  # For now assume all restricted areas will be full protection
  # mutate(restricted = ifelse(!is.na(activity_status_label), 1, 0)) |>
  # select(survey_abbrev, grouping_code, restricted, X, Y) |>
  st_drop_geometry()

fit_OS <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT S",
  formula = catch_count ~ 0 + fyear,
  species = sp,
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir
)
fit_ON <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT N",
  formula = catch_count ~ 0 + fyear,
  species = sp,
  spatiotemporal = "iid",
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir
)
fit_IN <- fit_hbll(dat = sp_dat, # didn't converge with spatiotemporal = "iid"
  survey_type = "HBLL INS N",
  formula = catch_count ~ 0 + fyear,
  species = sp,
  spatiotemporal = "off",
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir
)
meep()

# Just want to see what rho looks like for AR1
fit_OS_ar1 <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT S",
  formula = catch_count ~ 0 + fyear,
  species = sp,
  spatiotemporal = "ar1",
  use_extra_time = TRUE,
  time = "year",
  fit_dir = fit_dir,
  tag = "ar1"
)
fit_ON_ar1 <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT N",
  formula = catch_count ~ 0 + fyear,
  species = sp,
  spatiotemporal = "ar1",
  use_extra_time = TRUE,
  time = "year",
  fit_dir = fit_dir,
  tag = "ar1"
)

# TODO: evaluate and compare conditioning models: see - https://github.com/mis-assess/shrimp_surveydesign_csas/blob/794abdf0d4657dff5ed3316fe876b58afab0dd83/Reproducible_Examples/coastwide-density.R#L157

# -----------------------------------------------------------------------------
# Simulate data on HBLL grid for all three surveys
# ------------------------------------------------------------
# Do you ever make the starting value of a simulation the values of the last
# observed year so that things match up? Does this even matter or change things?

# Setup deviations on a trend, following an AR1 process

# This comment block is not worth anyone else reading. I am making better method notes elsewhere.
# # Me walking through all this... will move out later
# b_conditioning <- get_model_pars(fit_ON)
# # Use the sigma_E (marginal standard deviation of the spatiotemporal random field)
# # from the conditioning model as the target marginal sd of the AR1 process used
# # in simulation (marginal sd = sd of all the observations in a time-series)
# sigma_E_marginal_target <- b_conditioning$estimate[b_conditioning$term == "sigma_E"]
# rho <- 0.7 # Choose something for now but I'm not sure how we should
# # parameterize this. Look into values from literature? try fits to other species?
# # get from trawl fits?
# n_years <- 10

# # Marginal Variance = Conditional Variance / (1-ρ²)
# # Conditional Variance = (1-ρ²) × Marginal Variance

# # TMB/sdmTMB will then compute the innovation sd (conditional sd)
# # see: https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html
# # note to self - I still don't quite understand/follow that TMB documentation
# # note - phi in the TMB documentation is the AR1 rho
# # conditional sd: sd of innovations (annual deviations, i.e., 'conditioned' on previous year's value)
# # $Conditional Variance = (1 - \rho^2) × Marginal Variance$
# # $Conditional sd = \sqrt{(1 - \rho^2) × Marginal variance)}$
# # $Conditional sd = \sqrt{(1 - \rho^2)} × Marginal sd$
# sigma_E_innovation <- sqrt(1 - rho^2) * sigma_E_marginal_target  # AR1 innovation sd (conditional sd)
# # wikipedia:
# # "innovation is the difference between the observed value of a variable at time
# # t and the optimal forecast of that value based on information available prior
# # to time t."
# innovations <- rnorm(n_years) * sigma_E_innovation
# # TIL... you can just do this, which is why we get the innovations by doing rnorm(1) * innovation_sd
# # d1 <- rnorm(1000000, sd = 3)
# # d2 <- rnorm(1000000, sd = 1) * 3
# # sd(d1)
# # sd(d2)

imputed_sigma_E <- max(c(get_marginal_sigma_E(fit_ON), get_marginal_sigma_E(fit_OS)))

sim_IN <- simulate_hbll(fit_IN, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 0 + as.factor(year_covariate) + restricted * year_covariate,
  seed = 42,
  year_covariate = 1:20,
  mpa_trend = log(1.05),
  ar1_rho = 0.7,
  ar1_sigma_E = imputed_sigma_E,
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "ins-n"
) |>
  mutate(survey_abbrev = "HBLL INS N")

sim_ON <- simulate_hbll(fit_ON, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 0 + as.factor(year_covariate) + restricted * year_covariate,
  seed = 42,
  year_covariate = 1:20,
  mpa_trend = log(1.05),
  ar1_rho = 0.7,
  ar1_sigma_E = get_marginal_sigma_E(fit_ON),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "out-n"
) |>
  mutate(survey_abbrev = "HBLL OUT N")

sim_OS <- simulate_hbll(fit_OS, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 0 + as.factor(year_covariate) + restricted * year_covariate,
  seed = 42,
  year_covariate = 1:20,
  mpa_trend = log(1.05),
  ar1_rho = 0.7,
  ar1_sigma_E = get_marginal_sigma_E(fit_OS),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "out-s"
) |>
  mutate(survey_abbrev = "HBLL OUT S")

sim_dat0 <- bind_rows(sim_IN, sim_ON, sim_OS) |>
  left_join(hbll_grid |> select(X, Y, block_id, grouping_code)) |>
  select(!contains("as.factor(year_covariate)")) |>
  mutate(offset = 0) # Question: Yes??
sim_dat_sf <- XY_to_sf(sim_dat0)

ggplot(data = sim_dat_sf) +
  geom_sf(aes(colour = eta)) +
  geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
  scale_colour_viridis_c(option = "A") +
  facet_wrap(~ year) +
  get_plot_limits(sim_dat_sf, buffer = 5000)

# Look at diff between years
# diff_df <- sim_dat |>
#   group_by(X, Y, block_id) |>
#   arrange(year, .by_group = TRUE) |>
#   mutate(
#     mu_diff = mu - lag(mu),
#     eta_diff = eta - lag(eta),
#     percent_change = (exp(eta_diff) - 1) * 100,
#     observed_diff = observed - lag(observed)
#   ) |>
#   ungroup()

# p1 <- left_join(hbll_grid_poly, diff_df, by = "block_id")  |>
#   ggplot(data = _) +
#     geom_sf(data = pacea::bc_coast, fill = "grey90", colour = "grey90") +
#     geom_sf(data = mpa_shape_simplified, aes(fill = activity_allowance_label),
#       colour = NA, alpha = 0.8) +
#     scale_fill_manual(name = "activity_allowance_label",
#       values = c("allowed" = "#0072B2",
#         "not allowed" = "#D55E00", "conditional" = "#F0E442",
#         "not applicable" = "#999999"), na.value = "grey90") +
#     geom_sf(aes(colour = percent_change), alpha = 0.5) +
#     scale_colour_viridis_c(end = 0.8) +
#     theme(legend.position = "inside",
#           legend.position.inside = c(0.9, 0.2)) +
#     facet_wrap(~ year, nrow = 3) +
#     plot_limits_combined
# # p1
# p1 %+% (left_join(hbll_grid_poly, diff_df, by = "block_id") |> filter(year == 20)) +
#   theme(legend.position.inside = c(0.1, 0.2))

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

# HBLL INS N even years
# HBLL OUT N even years
# HBLL INS S odd years
sp_dat |> distinct(survey_abbrev, year)
hbll_year_filter <- function(df) {
  df |>
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "even") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "odd"))
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

samp1 <- sample_by_plan(
  sim_dat = sim_dat,
  sampling_effort = allocation_lu |> mutate(n_samps = allocation),
  grouping_vars = c("survey_abbrev", "year", "spatial_grouping_id", "strata_depth")) |>
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
  # maybe it is simplest to sample every year and just filter out years we don't actually sample?
  hbll_year_filter() |>
  mutate(plan = "status quo")
# samp1_summary <- samp1 |>
#   group_by(year, restricted) |>
#   summarise(n = n()) |>
#   filter(restricted == 1)
  # tidyr::pivot_wider(names_from = restricted, values_from = n, values_fill = 0)

plan_name <- "Status quo"
p1 <- local(plot_sampling_plan(samp1, plan_name))
dir.create(here::here("draft-figures"), showWarnings = FALSE)
ggsave(here::here("draft-figures", paste0("sim-dat-", plan_name, ".png")), width = 9, height = 18)

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

# Prep data for monitoring model
historical <- sp_dat |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  st_join(comm_ll_activity_status |> st_transform(crs = 3156), join = st_within) |>
  mutate(activity_status_label = if_else(is.na(activity_status_label), "outside", activity_status_label)) |>
  mutate(restricted = ifelse(activity_status_label == "outside", 0, 1)) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, fishing_event_id, latitude, longitude, X, Y,
    depth_m, offset,
    catch_count, restricted) |>
  mutate(after = 0)

simulated <- samp1 |>
  mutate(offset = 0) |>
  select(ssid = "survey_series_id", survey_abbrev, year_covariate, X, Y,
    # depth_m,
    offset,
    catch_count = "observed", restricted)


# Case 1 - fit historical with as.factor(year) and then use the last factor year
# as the intercept for the simulated data
test1 <- fit_monitoring(
  historical = historical |> filter(survey_abbrev == "HBLL OUT S"),
  simulated = simulated |> filter(survey_abbrev == "HBLL OUT S"),
  .formula = catch_count ~ as.factor(year_collapsed) + restricted * year_since_mpa,
  spatial = "on",
  spatiotemporal = "iid",
  time = "year_collapsed"
)

# Case 2 - use continuous time with breakpoint at stat of MPA
test2 <- fit_monitoring(
  historical = historical |> filter(survey_abbrev == "HBLL OUT S"),
  simulated = simulated |> filter(survey_abbrev == "HBLL OUT S"),
  .formula = catch_count ~ restricted * year_since_mpa,
  spatial = "on",
  spatiotemporal = "AR1",
  time = "year"
)

# Case 3 - use AR1 only post-MPA --> this didn't work;
test3 <- fit_monitoring(
  historical = historical |> filter(survey_abbrev == "HBLL OUT S"),
  simulated = simulated |> filter(survey_abbrev == "HBLL OUT S"),
  .formula = catch_count ~ 0 + after_mpa + restricted * year_since_mpa,
  spatial = "on",
  spatiotemporal = "AR1",
  time = "year",
  time_varying = ~ 0 + after_mpa,
  # extra_time = seq(2007:2042),
  time_varying_type = "ar1"
)

meep()


stop()

# Strategy 2: status quo + 10% additional effort
# ------------------------------------------------------------
samp2 <- sample_by_plan(
  sim_dat = sim_dat,
  sampling_effort = allocation_lu |> mutate(n_samps = round(allocation * 1.1)),
  grouping_vars = c("survey_abbrev", "year", "spatial_grouping_id", "strata_depth")) |>
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
  # maybe it is simplest to sample every year and just filter out years we don't actually sample?
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "even") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "odd")) |>
  mutate(plan = "status quo + 10% effort")

plan_name_2 <- "Status quo + 10% effort"
p2 <- local(plot_sampling_plan(samp2, plan_name_2))

# Look at what we have for baseline data in MPAs
# -----------------------------------------------------------------------------
# One question - how many unique blocks have been sampled in an MPA over the years?
# Welp, turns out that there are HBLL INS N locations that are not in the grid but
# that have sampled within MPAs, so... for now I will work without these, but I
# do think this is worth revisiting or at least making note of so that these
# locations could be resampled or included in future resampling for the MPAs.
# Most of the values that do not overlap the grid are the 47 and 48 grouping_code
# values. The southern ones are for the most part ones that should be part of the
# HBLL INS S grid.
# I guess none of this really matters but should be aware that that these 47 and
# and 48 are outside of the grid? Should we increase grid for the simulation?
# I guess this matters depending on how much we care about being able to account
# for the couple of relatively well-sampled MPA locations that fall outside of the
# existing grid. We would need to expand the prediction grid to be able to
# include these locations because right now we do not generate predictions at these
# locations.

# Side bar: showing previously sampled locations that are outside of the HBLL grid
sp_dat |>
  select(ssid, survey_abbrev, year, latitude, longitude, grouping_code) |>
  # select(ssid, year, latitude, longitude, grouping_code) |>
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326, crs_to = 32609) |>
  st_join(hbll_grid_poly, join = st_within) |>
  st_drop_geometry() |>
  # filter(is.na(block_id))
  # distinct(ssid, grouping_code.x, grouping_code.y, block_id) |>
  filter(is.na(block_id)) |>
  # filter(grouping_code.x < 317) |>
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326, crs_to = 32609) |>
  ggplot() +
    geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
    geom_sf(data = hbll_grid_poly, aes(fill = survey_abbrev), colour = NA, alpha = 0.3) +
    geom_sf(aes(colour = factor(grouping_code.x), shape = survey_abbrev.x), size = 2) +
    geom_sf(data = mpa_shape_simplified, fill = "black", colour = "black", alpha = 0.1) +
    geom_sf(data = mpa_shape_simplified |> filter(stringr::str_detect(category_simple, "RCA")), fill = "green", colour = "green", alpha = 0.1) +
    scale_colour_viridis_d() +
    # plot_limits_combined
    get_plot_limits(XY_to_sf(sim_dat |> filter(survey_abbrev %in% c("HBLL INS N"))), buffer = 50000)
ggsave(here::here("draft-figures", "hbll-ins-n-NA-grid-locations.pdf"), width = 10, height = 10)

# Identify which sampled locations are inside/outside MPA polygons
sp_dat_mpa_status <- sp_dat |>
  select(ssid, survey_abbrev, year, latitude, longitude, grouping_code) |>
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326,
    crs_to = st_crs(comm_ll_activity_status)) |>
  st_join(comm_ll_activity_status, join = st_within) |>
  mutate(
    in_mpa = ifelse(is.na(activity_status_label), 0, 1),
    activity_status_label = ifelse(is.na(activity_status_label), "outside", activity_status_label)
  ) |>
  st_drop_geometry()

janitor::tabyl(sp_dat_mpa_status, year, in_mpa) |>
  mutate(prop = round(`1` / (`1` + `0`), 2)) |>
  reframe(prop = mean(prop), n = mean(`1`))

sp_dat_mpa_status |>
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326, crs_to = 32609) |>
  ggplot() +
    geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
    geom_sf(data = mpa_shape_simplified, fill = "grey85", colour = "grey85") +
    geom_sf(aes(colour = factor(in_mpa)), shape = 21, size = 2) +
    # scale_shape_manual(name = "In MPA", values = c("0" = 21, "1" = 19)) +
    scale_colour_manual(name = "In MPA", values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                        labels = c("0" = "Outside", "1" = "Inside")) +
    plot_limits_combined
ggsave(here::here("draft-figures", "HBLL-sampled-inside-outside-mpa.pdf"), width = 10, height = 10)

# Get unique MPA locations with grid information
sampled_mpa_locations_with_grid <- sp_dat_mpa_status |>
  # select(-survey_abbrev, -grouping_code) |> # simplest to remove from one or the other dataframes
  XY_to_sf(x_col = "longitude", y_col = "latitude", mult = 1, crs_from = 4326, crs_to = 32609) |>
  st_join(hbll_grid_poly, join = st_within)

unique_samples_in_mpa <- bind_rows(
  sampled_mpa_locations_with_grid |>
    filter(!is.na(block_id)) |>
    st_drop_geometry() |>
    distinct(ssid, block_id, in_mpa, .keep_all = TRUE),
  sampled_mpa_locations_with_grid |>
    filter(is.na(block_id)) |>
    st_drop_geometry() |>
    distinct(ssid, in_mpa, .keep_all = TRUE)
)

unique_samples_in_mpa |>
  janitor::tabyl(in_mpa) |>
  filter(in_mpa == 1) |>
  pull(n) %>%
  cat("Number of unique samples in MPAs:", .)


# Strategy 3: Status quo + n additional MPA samples every 5 years - this could result in bias?
# ------------------------------------------------------------
# Get previously sampled MPA sites from historical data
historically_sampled_mpa_sites <- unique_samples_in_mpa |>
  filter(in_mpa == 1) |>
  # Because the strata information is in a messy state, let's make the strata for
  # these based on the recorded depths:
  mutate(strata_depth = case_when(depth_m >= 40 & depth_m < 70 ~ "40-70",
                                  depth_m >= 70 & depth_m < 150 ~ "70-150",
                                  depth_m >= 150 & depth_m < 259 ~ "150-259",
                                  TRUE ~ NA))

# Strategy 2:
# What if we double the number of locations sampled in MPAs every 5 years:
# * ~60 * 2 sites in MPAs every 5 years (almost equivalent to a whole doubling of HBLL survey effort)
# * stratify the sampling by depth
# * sample only in previously sampled locations
# * will want to try out more stratification but let's start with this
# * will want to probably target the MPAs likely to have full protection status
# * will want to try targetting MPAs that are restricting the most pre-implementation
#   fishing pressure

