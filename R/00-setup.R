# Core libraries
library(sdmTMB) # pak::pak("tidyverse/dplyr@778e413
library(dplyr)
library(ggplot2)
# library(here)
library(sf)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))

# Species data cache
if (Sys.info()['user'] == "jilliandunic") synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"
if (Sys.info()['user'] %in% c("dunic", "anderson")) synopsis_cache <- "/srv/anderson/src/gfsynopsis-2024/report/data-cache-2025-03"
if (Sys.info()['user'] == "seananderson") synopsis_cache <- "../gfsynopsis-2024/report/data-cache-2025-03"
if (Sys.info()['user'] == "jillian") synopsis_cache <- here::here("data-generated", "cleaned-species-data")
# so that there is a place to put some of the data dependencies
dir.create(here::here("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

hbll_ssids <- c(22, 36, 39, 40)
syn_ssids <- c(1, 3, 4, 16)

if (!file.exists(file.path("data-generated", "hbll-last-sampled-year.rds")) |
    !file.exists(file.path("data-generated", "historical-locations.rds"))) {
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
    saveRDS(hbll_last_sampled_year, file.path("data-generated", "hbll-last-sampled-year.rds"))

# Really slow on remote server (I think because old RGEOS and RGDAL?)
  historical_locations <- d0 |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
    filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
    add_utm_columns() |>
    XY_to_sf(crs_to = 32609) |>
    st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
    mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
    st_join(hbll_grid_poly |> select(block_id, grouping_code), join = st_within) |>
    st_drop_geometry() |>
    # select(ssid = survey_series_id.x, survey_abbrev, year, fishing_event_id,
    #   latitude, longitude, X, Y,block_id,
    #   fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y, restricted)
    distinct(survey_abbrev, block_id, latitude, longitude, X, Y)
  saveRDS(historical_locations, file.path("data-generated", "historical-locations.rds"))
}

survey_lu <- tibble::tibble(
  survey_abbrev = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG",
                    "HBLL OUT N", "HBLL OUT S",
                    "HBLL INS N", "HBLL INS S"),
  survey_series_id = c(1, 3, 4, 16, 22, 36, 39, 40)
)

if (!file.exists(here::here("data-generated", "spatial", "simple-mpa.rds"))) {
  source(here::here("R", "01-prepare-spatial-data.R"))
} else {
  simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))
}

if (!file.exists(file.path("data-generated", "hbll-restricted-sf.rds"))) {
  gfdata::load_survey_blocks(type = "XY") |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
    XY_to_sf(crs_to = 32609) |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) %>%
    st_join(., simple_mpa |> st_transform(crs = st_crs(.)), join = st_within) |>
    mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
  saveRDS(file.path("data-generated", "hbll-restricted-sf.rds"))
}

# Setup allocations
if (!file.exists(file.path("data-generated", "hbll-allocations.rds"))) {
  source(here::here("R", "prep-survey-allocations.R"))
}

# Something like this could be helpful to add later
# make load data file - sean included this which was smart
# dir.create("data-generated", showWarnings = FALSE)
# dir.create("figs", showWarnings = FALSE)

# if (Sys.info()[["user"]] != "seananderson") {
#   stop("This file does not need to be run; all outputs have been cached and commited in Git.")
# }

