library(sdmTMB)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
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


source(here::here("R", "00-fit-predict-plot.R"))

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

# plot predictions
plot_hbll_predictions(pred_ins |> filter(year %in% 2024),
  rotation = 90) + facet_wrap(~ year)
plot_hbll_predictions(pred_out |> filter(year %in% 2022), rotation = 90) + facet_wrap(~ year)

# run simulation
