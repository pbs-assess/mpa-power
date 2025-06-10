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
#     1. use parameters from that model to simulate new data with observations at the actual historically observed locations
#     1. when simulating, simulate recovery at some rate within closed areas and a stationary abundance/density outside of closed areas

# - Dimensions that will likely affect the answer:
#     1. species (therefore estimated spatial and spatiotemporal SD, spatial correlation range, observation error)
#     1. rate of 'recovery'
#     1. number of years observed
#     1. whether one assesses all restricted areas together or individually

# sizes of restricted areas
# how much sampling within restricted areas
#

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
fit_out1 <- fit_hbll(dat = sp_dat, survey_type = "HBLL OUT",
                    species = "yelloweye-rockfish",
                    fit_dir = fit_dir)
# Was trying to understand why this plot looked so different from gfsynopsis
# I think it's because it's plotted on link scale, and I was plotting on response scale
fit_out2 <- fit_hbll(dat = sp_dat, survey_type = "HBLL OUT S",
                    species = "yelloweye-rockfish",
                    fit_dir = fit_dir,
                    tag = "gfsynopsis",
                    formula = catch_count ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                    spatial = "on",
                    spatiotemporal = "off")

# get predictions
pred_ins <- predict_hbll(fit_ins, hbll_grid, re_form = NULL)

pred_out1 <- predict_hbll(fit_out1, hbll_grid, re_form = NULL)
pred_out2 <- predict_hbll(fit_out2,
  grid = hbll_grid |>
    filter(str_detect(survey_abbrev, "HBLL OUT S")),
  re_form = NULL)

# plot predictions
plot_hbll_predictions(pred_ins |> filter(year %in% 2024),
  rotation = 90, type = "link") + facet_wrap(~ year)
plot_hbll_predictions(pred_out1 |> filter(year %in% 2022), rotation = 90) + facet_wrap(~ year)
plot_hbll_predictions(pred_out2 |> filter(year %in% 2022), rotation = 90,
  type = "link")
plot_hbll_predictions(pred_out2 |> filter(year %in% 2022), rotation = 90,
  type = "response")

# run simulation
