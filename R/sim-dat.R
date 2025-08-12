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
  left_join(hbll_grid |> select(X, Y, block_id, grouping_code))
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

# Status quo sampling effort
# ------------------------------------------------------------
# FIXME - these are broken right now - they are empty files
sampled_cpue_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "sampled-cpue-mpa-overlap-yelloweye.rds"))
ll_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "ll-mpa-overlap-yelloweye.rds"))


hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))

sim_dat <- left_join(sim_dat0, hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
  mutate(spatial_grouping_id = ifelse(pfma %in% c("5A", "4B"), "5A4B", pfma)) # Group 5A and 4B together for sampling purposes


allocation_lu <- sim_dat |>
  distinct(survey_abbrev, spatial_grouping_id, strata_depth, allocation)

# HBLL INS N even years
# HBLL OUT N even years
# HBLL INS S odd years
sp_dat |> distinct(survey_abbrev, year)

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
  filter((survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "even") |
        (survey_abbrev %in% c("HBLL OUT S") & odd_even == "odd")) |>
  mutate(plan = "status quo")
# samp1_summary <- samp1 |>
#   group_by(year, restricted) |>
#   summarise(n = n()) |>
#   filter(restricted == 1)
  # tidyr::pivot_wider(names_from = restricted, values_from = n, values_fill = 0)

plan_name <- "Status quo"
p1 <- local(plot_sampling_plan(samp1, plan_name))
ggsave(here::here("draft-figures", paste0("sim-dat-", plan_name, ".png")), width = 9, height = 18)

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
