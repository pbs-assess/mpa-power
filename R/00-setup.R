# Core libraries
library(sdmTMB) # pak::pak("sdmTMB/sdmTMB@6e9e7de") --> I think this is what is needed for the resdoc???
# sdmTMB     0.8.1.9001 → 1.0.0.9015 👷🔧 (GitHub: 6e9e7de), unless I updated since then????
library(ggplot2)
# library(here)
library(sf)
library(dplyr)
library(future)
library(furrr)

# Handle namespace conflicts
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))
source(here::here("R", "sample-fit-config.R"))

synopsis_cache <- set_synopsis_cache()
# so that there is a place to put some of the data dependencies
dir.create(here::here("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

ms_dir <- here::here("data-generated", run_tag)
dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)

mesh_dir         <- file.path(ms_dir, "00-mesh-cache")
fit_dir          <- file.path(ms_dir, "01-fits")
cleaned_data_dir <- file.path(ms_dir, "01-cleaned-species-data")
sim_dir          <- file.path(ms_dir, "02-sim-data")
sample_dir       <- file.path(ms_dir, "03-sampled-data")
results_dir      <- file.path(ms_dir, "04-power-results")

# ms-specific shared data files
historical_locations_file   <- file.path(ms_dir, "historical-locations.rds")
hbll_last_sampled_year_file <- file.path(ms_dir, "hbll-last-sampled-year.rds")
hbll_restricted_sf_file     <- file.path(ms_dir, "hbll-restricted-sf.rds")
hbll_allocations_file       <- file.path(ms_dir, "hbll-allocations.rds")
ar1_parameters_file         <- file.path(ms_dir, "ar1-parameters.rds")
fit_characteristics_file    <- file.path(ms_dir, "fit-characteristics.rds")
# TODO: consolidate power-df-historical-sampling.rds and power-results-df.rds
# (two filenames used across 06-figures-power-plots.R, 06-figure-power-correlations*.R,
# and 07-get-values.R) and add a power_results_file path variable here

fig_dir <- here::here("figures-ms")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

hbll_ssids <- c(22, 36, 39, 40)
syn_ssids <- c(1, 3, 4, 16)

survey_lu <- tibble::tibble(
  survey_abbrev = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG",
                    "HBLL OUT N", "HBLL OUT S",
                    "HBLL INS N", "HBLL INS S"),
  survey_series_id = c(1, 3, 4, 16, 22, 36, 39, 40)
)

if (!file.exists(file.path("data-generated", "hbll-dem-grid-depths.rds"))) {
  stop("data-generated/hbll-dem-grid-depths.rds is missing; run R/dem-data.R first; requires BC DEM raw file (large)")
}
grid_dem_depths <- readRDS(file.path("data-generated", "hbll-dem-grid-depths.rds")) |>
  left_join(survey_lu, by = "survey_series_id") |>
  select(survey_abbrev, block_id, dem_depth = depth_dem_mean) |>
  mutate(log_depth = log(dem_depth))

if (!file.exists(hbll_last_sampled_year_file) |
    !file.exists(historical_locations_file)) {
  d0 <- readRDS(file.path(synopsis_cache, "yelloweye-rockfish.rds"))$survey_sets
  simple_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa.rds"))
  hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
    filter(survey_series_id %in% c(22, 36, 39))
  hbll_last_sampled_year <- d0 |>
    filter(survey_series_id.x %in% hbll_ssids) |>
    group_by(survey_abbrev) |>
    slice(which.max(year)) |>
    select(ssid = survey_series_id.x, survey_abbrev, last_sampled_year = year) |>
    filter(ssid != 40)
  saveRDS(hbll_last_sampled_year, hbll_last_sampled_year_file)

# Really slow on remote server (I think because old RGEOS and RGDAL?); save it
# on local machine to prevent doing this on the server.
  historical_locations <- d0 |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
    filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
    add_utm_columns() |>
    XY_to_sf(crs_to = 32609) |>
    st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
    mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
    st_join(hbll_grid_poly |> select(block_id, grouping_code), join = st_within) |>
    st_drop_geometry() |>
    left_join(grid_dem_depths) |>
    distinct(survey_abbrev, block_id, uid, fishing_event_id, latitude, longitude, X, Y, year, restricted, dem_depth, log_depth)
  saveRDS(historical_locations, historical_locations_file)
}

if (!file.exists(here::here("data-generated", "spatial", "simple-mpa.rds"))) {
  source(here::here("R", "01-prepare-spatial-data.R"))
} else {
  simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))
  display_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa-500m.rds"))
}



if (!file.exists(hbll_restricted_sf_file)) {
  gfdata::load_survey_blocks(type = "XY") |>
    XY_to_sf(crs_to = 32609) |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
    st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
    mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
    select(-depth_m) |> # remove original depth and use DEM depth
    left_join(grid_dem_depths) |>
    tidyr::drop_na(log_depth) |>
  saveRDS(hbll_restricted_sf_file)
}

# Setup allocations
if (!file.exists(hbll_allocations_file)) {
  source(here::here("R", "prep-survey-allocations.R"))
}

# Something like this could be helpful to add later
# make load data file - sean included this which was smart
# dir.create("data-generated", showWarnings = FALSE)
# dir.create("figs", showWarnings = FALSE)

# if (Sys.info()[["user"]] != "seananderson") {
#   stop("This file does not need to be run; all outputs have been cached and commited in Git.")
# }

