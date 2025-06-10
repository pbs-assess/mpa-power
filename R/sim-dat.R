library(sdmTMB)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(purrr)
theme_set(ggsidekick::theme_sleek())


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


source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-utils.R"))

# Using data from synopsis cache
synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"

fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

hbll_grid <- load_survey_blocks(type = "coords") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(depth_m > 0) |>
  mutate(depth_mean = mean(log(depth_m), na.rm = TRUE),
         depth_sd = sd(log(depth_m), na.rm = TRUE),
         depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1],
         depth_scaled2 = depth_scaled^2)

bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
sp_dat0 <- readRDS(file.path(synopsis_cache, "yelloweye-rockfish.rds"))$survey_sets

sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  rename(ssid = "survey_series_id.x") |>
  left_join(bait_counts, by = c("year", "fishing_event_id", "ssid")) |>
  mutate(count_bait_only = replace(count_bait_only, which(count_bait_only == 0), 1),
    prop_bait_hooks = count_bait_only / hook_count,
    hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
    prop_removed = 1 - prop_bait_hooks,
    offset = log(hook_count / hook_adjust_factor),
    depth_mean = mean(log(depth_m), na.rm = TRUE),
    depth_sd = sd(log(depth_m), na.rm = TRUE),
    depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1],
    depth_scaled2 = depth_scaled^2
  ) |>
  sdmTMB::add_utm_columns()

# fit conditioning model
fit_ins <- fit_hbll(dat = sp_dat, survey_type = "HBLL INS",
                    species = "yelloweye-rockfish",
                    fit_dir = fit_dir)
fit_out <- fit_hbll(dat = sp_dat, survey_type = "HBLL OUT",
                    species = "yelloweye-rockfish",
                    fit_dir = fit_dir)

# get predictions
pred_ins <- predict_hbll(fit_ins, hbll_grid, re_form = NULL)
pred_out <- predict_hbll(fit_out, hbll_grid, re_form = NULL)

plot predictions
plot_hbll_predictions(pred_ins |> filter(year %in% 2024),
  rotation = 90) + facet_wrap(~ year)
plot_hbll_predictions(pred_out |> filter(year %in% 2022), rotation = 90) + facet_wrap(~ year)

# run simulation
# ------------------------------------------------------------
# start with just simulating coastwide data into the future
# - no trend
# - no MPA effect
# - outside survey only

# Start with outside:

fit <- fit_hbll(dat = sp_dat, survey_type = "HBLL INS",
                    species = "yelloweye-rockfish",
                    fit_dir = fit_dir)


# get the model parameters
b <- get_model_pars(fit)
last_pars_best <- fit$tmb_obj$env$last.par.best

# fixed random effects
omega_s <- get_pars(fit)$omega_s # or: fixed_omega_s <- last_pars_best[names(last_pars_best) == "omega_s"] # TMB last best estimate
epsilon_st <- get_pars(fit)$epsilon_st

# random effect SDs
omega_s_sd <- b$estimate[b$term == "sigma_O"]
epsilon_st_sd <- b$estimate[b$term == "sigma_E"]


# input df to simulate with
year_covariate <- get_model_years(fit) - min(get_model_years(fit))
# dat <- sdmTMB::replicate_df(
#   hbll_grid |>
#     filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
#     mutate(survey_abbrev = "HBLL OUT") |>
#     distinct(),
#   time_name = "year_covariate",
#   time_values = year_covariate) |>
#   mutate(year = year_covariate)


dat <- fit$data
attr(dat, "version") <- NULL
attr(dat, "date") <- NULL

input_dat <- dat |> select(X, Y) |>
  sdmTMB::replicate_df(
    time_name = "year_covariate",
    time_values = year_covariate) |>
  mutate(year = year_covariate)
input_mesh <- make_mesh(input_dat, xy_cols = c("X", "Y"), cutoff = 10) # find way to get cutoff value from fit

.seed = 714
sim_dat <- sdmTMB::sdmTMB_simulate(
  # formula = ~ 1 + restricted * year_covariate,
  formula = ~ 1 + year_covariate,
  data = input_dat,
  mesh = input_mesh,
  family = nbinom2(),
  time = "year",
  # rho = 2 * plogis(m$model$par[['ar1_phi']]) - 1, # TODO!?
  rho = NULL,
  sigma_O = b$estimate[b$term == "sigma_O"], # Question - is this just ignored if fixed_re omega_s is set?
  sigma_E = b$estimate[b$term == "sigma_E"],
  phi = b$estimate[b$term == "phi"],
  range = b$estimate[b$term == "range"],
  fixed_re = list(omega_s = omega_s, epsilon_st = NULL, zeta_s = NULL),
  # (Intercept), restrictedTRUE, year_covariate, restrictedTRUE:year_covariate
  # B = c(mean(b[grep("year", b$term),"estimate"]), 0, 0, mpa_increase_per_year),
  B = c(0, 0.1), # Intercept and year effect (0.1 is a small positive trend)
  seed = .seed
)
sim_dat |> glimpse()

ggplot(sim_dat) +
  geom_point(aes(x = X, y = Y, color = omega_s)) +
  scale_color_viridis_c() +
  facet_wrap(~ year_covariate)

ggplot(sim_dat) +
  geom_point(aes(x = X, y = Y, color = eta)) +
  scale_color_viridis_c() +
  facet_wrap(~ year_covariate)