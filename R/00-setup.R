# Core libraries
library(sdmTMB) # pak::pak("tidyverse/dplyr@778e413
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

synopsis_cache <- set_synopsis_cache()
# so that there is a place to put some of the data dependencies
dir.create(here::here("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

run_tag <- "-ms"  # set to "" for resdoc outputs
mesh_dir <- here::here("data-generated", paste0("00-mesh-cache", run_tag))
fit_dir    <- here::here("data-generated", paste0("01-fits", run_tag))
cleaned_data_dir <- here::here("data-generated", paste0("01-cleaned-species-data", run_tag))
sim_dir    <- here::here("data-generated", paste0("02-sim-data", run_tag))
sample_dir <- here::here("data-generated", paste0("03-sampled-data", run_tag))
results_dir <- here::here("data-generated", paste0("04-power-results", run_tag))

fig_dir <- here::here("figures-ms")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

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
    # select(ssid = survey_series_id.x, survey_abbrev, year, fishing_event_id,
    #   latitude, longitude, X, Y,block_id,
    #   fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y, restricted)
    distinct(survey_abbrev, block_id, uid, fishing_event_id, latitude, longitude, X, Y, restricted)
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

if (!file.exists(file.path("data-generated", "hbll-dem-grid-depths.rds"))) {
  stop("data-generated/hbll-dem-grid-depths.rds is missing; run R/dem-data.R first")
}
grid_dem_depths <- readRDS(file.path("data-generated", "hbll-dem-grid-depths.rds")) |>
  select(survey_series_id, block_id, depth_m = depth_dem_mean) |>
  mutate(log_depth = log(depth_m))

if (!file.exists(file.path("data-generated", "hbll-restricted-sf.rds"))) {
  gfdata::load_survey_blocks(type = "XY") |>
    XY_to_sf(crs_to = 32609) |>
    filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
    st_join(simple_mpa |> st_transform(crs = 32609), join = st_within) |>
    mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
    select(-depth_m) |> # remove original depth and use DEM depth
    left_join(grid_dem_depths, by = c("survey_series_id", "block_id")) |>
    tidyr::drop_na(log_depth) |>
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

