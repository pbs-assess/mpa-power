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
  filter(survey_abbrev != "HBLL INS S") |>
  mutate(cell_id = row_number())
  # depth not currently used in sim - I think I need to fix the depth data before including
  # filter(depth_m > 0) |>
  # mutate(depth_mean = mean(log(depth_m), na.rm = TRUE),
  #        depth_sd = sd(log(depth_m), na.rm = TRUE),
  #        depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1],
        #  depth_scaled2 = depth_scaled^2)
hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |>
  mutate(cell_id = row_number())

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

# TODO: evaluate and compare conditioning models: see - https://github.com/mis-assess/shrimp_surveydesign_csas/blob/794abdf0d4657dff5ed3316fe876b58afab0dd83/Reproducible_Examples/coastwide-density.R#L157

# -----------------------------------------------------------------------------
# Simulate data on HBLL grid for all three surveys
# ------------------------------------------------------------
sim_IN <- simulate_hbll(fit_IN, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 1 + restricted * year_covariate,
  seed = 42,
  year_covariate = seq(from = 1, to = 21, by = 2),
  mpa_trend = log(1.05),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "ins-n"
) |>
  mutate(survey_abbrev = "HBLL INS N")

sim_ON <- simulate_hbll(fit_ON, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 1 + restricted * year_covariate,
  seed = 42,
  year_covariate = seq(from = 0, to = 20, by = 2),
  mpa_trend = log(1.05),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "out-n"
) |>
  mutate(survey_abbrev = "HBLL OUT N")

sim_OS <- simulate_hbll(fit_OS, restricted_df,
  sim_dir = "data-generated/sim-dat",
  check_cache = TRUE,
  formula = ~ 1 + restricted * year_covariate,
  seed = 42,
  year_covariate = seq(from = 1, to = 21, by = 2),
  mpa_trend = log(1.05),
  fixed_spatial_re = TRUE,
  fixed_spatiotemporal_re = FALSE,
  tag = "out-s"
) |>
  mutate(survey_abbrev = "HBLL OUT S")

sim_dat <- bind_rows(sim_IN, sim_ON, sim_OS) |>
  left_join(hbll_grid |> select(X, Y, block_id, grouping_code))
sim_dat_sf <- XY_to_sf(sim_dat)

# TODO: Add deviations from trend (need if we want to simulate more years than sampled no?)
# - Simulate future years with a linear effect of time.
# - But we could get fancier and have deviations from a trend,
#   in which case you’d pre-simulate those and pass them in a coefficients on future factor years.
# - Having deviations from a mean trend would be most realistic
#   but adds more variables to pick like the level of autocorrelation and SD of
#   those deviations.

ggplot(data = sim_dat_sf) +
  geom_sf(aes(colour = eta)) +
  geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
  scale_colour_viridis_c(option = "A") +
  facet_wrap(~ year) +
  get_plot_limits(sim_dat_sf, buffer = 5000)

# Look at diff between years
# diff_df <- sim_dat |>
#   group_by(X, Y, cell_id) |>
#   arrange(year, .by_group = TRUE) |>
#   mutate(
#     mu_diff = mu - lag(mu),
#     eta_diff = eta - lag(eta),
#     percent_change = (exp(eta_diff) - 1) * 100,
#     observed_diff = observed - lag(observed)
#   ) |>
#   ungroup()

# p1 <- left_join(hbll_grid_poly, diff_df, by = "cell_id")  |>
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
# p1 %+% (left_join(hbll_grid_poly, diff_df, by = "cell_id") |> filter(year == 20)) +
#   theme(legend.position.inside = c(0.1, 0.2))

# Status quo sampling effort
# ------------------------------------------------------------
sampled_cpue_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "sampled-cpue-mpa-overlap-yelloweye.rds"))
ll_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "ll-mpa-overlap-yelloweye.rds"))

hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
# NOTE: Need to group by ssid, pfma, and strata_depth because of the 5A4B shenanigans

left_join(sim_dat, hbll_allocations, by = c("survey_abbrev", "grouping_code")) |>
glimpse()

# - random sampling of blocks from year to year
samp1 <- sp_dat |>
  group_by(survey_abbrev, year) |>
  summarise(n = n()) |>
  group_by(survey_abbrev) |>
  summarise(n_samps = round(mean(n))) |>
  mutate(plan = "status quo")

sampled_sim_dat1 <- sample_by_plan(sim_dat, samp1,
  grouping_vars = c("survey_abbrev", "year"))

# - resample previously sampled MPA cells if sampling within MPA
# - resample at current rate within MPAs

samp2 <- sp_dat |>
  group_by(survey_abbrev, year, restricted) |>
  summarise(n = n()) |>
  group_by(survey_abbrev, restricted) |>
  summarise(n_samps = round(mean(n))) |>
  mutate(plan = "resample MPA blocks at current rate")

sampled_sim_dat2 <- sample_by_plan(sim_dat, samp2,
  grouping_vars = c("survey_abbrev", "year", "restricted"))

# - resample every previously sampled MPA cell every five years?
# (this method probably doesn't make sense for every MPA site because some have
# quite good coverage, so would need to come up with some kind of rule perhaps)
all_mpa_blocks <- sp_dat |>
  distinct(survey_abbrev, restricted, X, Y) |>
  filter(restricted == 1) |>
  mutate(plan = "resample every previously sampled MPA cell every five years")




restricted_df |>
  distinct(survey_abbrev, restricted) |>
  mutate(schedule = ifelse(restricted == 1, year - (year %% 5)) |>
  distinct(survey_abbrev, census_year, restricted)

sp_dat |> distinct(survey_abbrev, restricted, block_id)
samp3 <- sp_dat |>

sp_dat |>
distinct(survey_abbrev, year, restricted)

census_schedule <- sim_dat |>
  distinct(year, survey_abbrev) |>
  arrange(year) |>
  slice(2:4) |>
  rename(first_census_year = year)

test <- sample_by_plan2(sim_dat,
  sampling_effort = samp3,
  grouping_vars = c("survey_abbrev", "year", "restricted"),
  mpa_schedule = 5
)


sampled_sim_dat <- bind_rows(sampled_sim_dat1, sampled_sim_dat2)

# # So few actual fish observed is this right??
# sampled_sim_dat <- sim_dat |>
#   left_join(sampling_effort, by = "survey_abbrev") |> glimpse()
#   group_by(year) |>
#   sample_n(size = sampling_effort, replace = FALSE)

# sampling_plan <- "status quo"
sampling_plan <- "resample MPA blocks at current rate"

n_mpa_samps <- sampled_sim_dat |>
  filter(restricted == 1) |>
  filter(plan == sampling_plan) |>
  group_by(year) |>
  summarise(n = n())

# plot_limits <- get_plot_limits(XY_to_sf(sampled_sim_dat), buffer = 5000)
plot_years <- c(0, 1, 5, 6, 10, 11, 15, 16, 20, 21)
plot_dat <- sampled_sim_dat |>
  filter(plan == sampling_plan) |>
  XY_to_sf() |>
  filter(year %in% plot_years)

ggplot(data = plot_dat) +
  geom_sf(data = mpa_shape_simplified, aes(fill = activity_allowance_label),
    colour = NA, alpha = 0.3) +
  scale_fill_manual(name = "MPA allowance",
    values = c("allowed" = "#0072B2",
      "not allowed" = "#D55E00", "conditional" = "#F0E442",
      "not applicable" = "#999999"), na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
  ggnewscale::new_scale_fill() +
  # geom_sf(data = plot_dat |> filter(restricted == 1), colour = "black", size = 1.1) +
  geom_sf(aes(colour = eta, shape = factor(restricted)), size = 1.2) +
  scale_shape_manual(name = "Restricted",values = c(`0` = 21, `1` = 19)) +
  scale_colour_viridis_c(name = "eta", option = "A", end = 0.8) +
  facet_wrap(~ year, nrow = 5) +
  plot_limits_combined +
  # theme(legend.position = "inside",
  #   legend.position.inside = c(0.87, 0.1),
  #   legend.box = "horizontal") +
  geom_sf_text(data = n_mpa_samps |> mutate(X = -126, Y = 54) |> filter(year %in% plot_years) |>
    XY_to_sf(crs_from = 4326, crs_to = 4326),
    aes(x = X, y = Y, label = paste0("n = ", n)), size = 3) +
  ggtitle(sampling_plan)
# ggsave(here::here("draft-figures", paste0("sim-dat-", sampling_plan, ".pdf")), width = 14.5, height = 12)
ggsave(here::here("draft-figures", paste0("sim-dat-", sampling_plan, ".pdf")), width = 9, height = 18)
ggsave(here::here("draft-figures", paste0("sim-dat-", sampling_plan, ".png")), width = 9, height = 18)

# Refit to detect change:
sim_mesh <- make_mesh(sampled_sim_dat, xy_cols = c("X", "Y"), cutoff = 10)
sim_fit <- sdmTMB::sdmTMB(
  formula = observed ~ 1 + as.factor(year) * restricted,
  data = sampled_sim_dat,
  mesh = sim_mesh,
  family = nbinom2(link = "log"),
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  offset = NULL, # Question: Does the offset matter for sampling of the simulated data?
  anisotropy = TRUE
)
meep()
summary(sim_fit)

test <- glm_fit <- MASS::glm.nb(
  formula = observed ~ 1 + as.factor(year) * restricted,
  data = sampled_sim_dat,
  link = "log"
)
summary(test)

glm(est ~ type*year, data = ind, family = Gamma(link = "log"))

glm(observed ~ 1 + year * restricted, data = sampled_sim_dat) |> summary()

# Compare estimates to simulated data



# Simulate data on HBLL grid - single area test
# ------------------------------------------------------------
fit <- fit_IN

# get the model parameters
b <- get_model_pars(fit)

# fixed random effects (get single draw from rf distributions)
osp <- one_sample_posterior(fit)
omega_s <- osp[grepl("omega_s", names(osp))] |> matrix()

# random effect SDs
omega_s_sd <- b$estimate[b$term == "sigma_O"]
epsilon_st_sd <- b$estimate[b$term == "sigma_E"]
if (fit$spatiotemporal == "off") epsilon_st_sd <- 0

# input df to simulate with
# Question: if using simulate and fyear, this is not going into future?
# year_covariate <- get_model_years(fit) - min(get_model_years(fit))
year_covariate <- seq(from = 0, to = 20, by = 2) # Question: I don't fully understand how to include time, does time index from 0?

input_dat <- restricted_df |>
  filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
  select(X, Y, restricted) |>
  sdmTMB::replicate_df(
    time_name = "year_covariate", # why doesn't this need to be coded as a factor?
    time_values = year_covariate) |>
  mutate(restricted = restricted,
         year = as.numeric(year_covariate),
  )

input_mesh <- make_mesh(input_dat, xy_cols = c("X", "Y"), mesh = fit$spde$mesh)

mpa_trend <- log(1.05) # 5% increase per year
.seed <- 714
sim_dat <- sdmTMB::sdmTMB_simulate(
  formula = ~ 1 + restricted * year_covariate,
  data = input_dat,
  mesh = input_mesh,
  family = nbinom2(link = "log"),
  time = "year",
  # rho = 2 * plogis(m$model$par[['ar1_phi']]) - 1, # If AR process used: "Spatiotemporal correlation between years; between -1 and 1"
  rho = NULL,
  # sigma_O = b$estimate[b$term == "sigma_O"],
  sigma_E = epsilon_st_sd,
  phi = b$estimate[b$term == "phi"],
  range = b$estimate[b$term == "range"],
  fixed_re = list(omega_s = omega_s, epsilon_st = NULL, zeta_s = NULL),
  # (Intercept), restrictedTRUE, year_covariate, restrictedTRUE:year_covariate
  B = c(mean(b[grep("year", b$term),"estimate"]), 0, 0, mpa_trend), # every year has same mean yes? so years as factors being same as historical doesn't matter?
  seed = .seed
) |> as_tibble()

sim_dat |> glimpse()

# Start with sampling effort equivalent to HBLL OUT N (mean blocks per year)
# - random sampling of blocks from year to year
sampling_effort <- sp_dat |> filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
  group_by(year) |>
  summarise(n = n()) |>
  pull(n) |>
  mean()

# So few actual fish observed is this right??
sampled_sim_dat <- sim_dat |>
  group_by(year) |>
  sample_n(size = sampling_effort, replace = FALSE)
# plot_limits <- coord_sf(xlim = c(-134, -124), ylim = c(50, 54.5), crs = 4326)
# plot_limits <- coord_sf(xlim = c(-134, -124), ylim = c(51.5, 54.5), crs = 4326)

n_mpa_samps <- sampled_sim_dat |>
  filter(restricted == 1) |>
  group_by(year) |>
  summarise(n = n())

plot_limits <- get_plot_limits(XY_to_sf(sampled_sim_dat), buffer = 5000)
ggplot() +
  geom_sf(data = mpa_shape_simplified, aes(fill = activity_allowance_label),
    colour = NA, alpha = 0.3) +
  scale_fill_manual(name = "activity_allowance_label",
    values = c("allowed" = "#0072B2",
      "not allowed" = "#D55E00", "conditional" = "#F0E442",
      "not applicable" = "#999999"), na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = sampled_sim_dat |> XY_to_ll(), aes(colour = eta, shape = factor(restricted))) +
  scale_shape_manual(values = c(`0` = 21, `1` = 19)) +
  scale_colour_viridis_c(option = "A") +
  facet_wrap(~ year) +
  plot_limits +
  theme(legend.position = "top") +
  geom_sf_text(data = n_mpa_samps |> mutate(X = -126, Y = 50.1) |>
    XY_to_sf(crs_from = 4326, crs_to = 4326),
    aes(x = X, y = Y, label = paste0("n = ", n)), size = 3)





# Refit to detect change:
sim_mesh <- make_mesh(sampled_sim_dat, xy_cols = c("X", "Y"), cutoff = 10)
sim_fit <- sdmTMB::sdmTMB(
  formula = observed ~ 1 + as.factor(year) * restricted,
  data = sampled_sim_dat,
  mesh = sim_mesh,
  family = nbinom2(link = "log"),
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  offset = NULL, # Question: Does the offset matter for sampling of the simulated data???
  anisotropy = TRUE
)
meep()
summary(sim_fit)

# Compare estimates to simulated data



# I don't understand how to use project
# Try out using project instead of simulate, but need to add observation error after the fact.
# Not necessary since years as factors?
historical_years <- get_model_years(fit)
future_years <- seq(max(historical_years) + 2, max(historical_years) + 10, by = 2)
all_years <- c(historical_years, future_years)

proj_grid <- restricted_df |>
  filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
  select(X, Y, restricted) |>
  sdmTMB::replicate_df(
    time_name = "year", # why doesn't this need to be coded as a factor?
    time_values = all_years) |>
  mutate(restricted = restricted,
         fyear = factor(year)
  )

test <- project(
  object = fit,
  newdata = proj_grid,
  nsim = 1,
  uncertainty = "both", # both would be single draw from rf distributions?
  silent = FALSE,
  sims_var = "eta_i",
  # sp fields, st fields, svc fields, random intercepts, tvc, smoothers
  sim_re = c(0, 1, 0, 0, 1, 0),
  return_tmb_report = FALSE
)
