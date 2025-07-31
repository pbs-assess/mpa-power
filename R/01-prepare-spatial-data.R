library(dplyr)
library(ggplot2)
library(sf)

# Prepare spatial data (MPA polygons and human use layers)
dir.create(file.path("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

# Look up table for activity status
activity_status_lu <- tibble::enframe(c(
  X = "currently restricted",
  AC = "identified as concern",
  O = "not identified as concern",
  na = "not applicable")) |>
  rename(activity_status = name, activity_status_label = value)
saveRDS(activity_status_lu, file.path("data-generated", "spatial", "activity_status_lu.rds"))

# MPA polygons (double check what layer to use - I am assuming the most up to date one)
Q1 <- st_read(here::here("data-raw", "spatial", "Spatial_Q.gdb"),
  layer = "Q1_FULL_March2023") |>
  janitor::clean_names()

Q2 <- st_read(here::here("data-raw", "spatial", "All_Network_boundaries_Q2_2024.gdb"),
  layer = "All_Network_boundaries_Q2_2024") |>
  janitor::clean_names()

shape <- Q2

comm_ll_activity_status <- shape |>
  select(hu_commercial_harvest_bottom_longline_demersal_hookand_line,
    category_detailed, category_simple) |>
  left_join(activity_status_lu,
    by = c("hu_commercial_harvest_bottom_longline_demersal_hookand_line" = "activity_status")) |>
  sf::st_cast("MULTIPOLYGON") %>%
  mutate(mpa_id = row_number(),
         mpa_area = st_area(.))

saveRDS(comm_ll_activity_status, file.path("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))

# Human use layers
# human_layers <- st_layers(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"))
# human_layers$name

# hu_ll <- st_read(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"),
#   layer = "hu_co_demersalfishing_bottomlongline_d")


# -----------------------------------------------------------------------------
# Commercial longline data - raw
# -----------------------------------------------------------------------------

synopsis_cache <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"

# Load for all species - filter as necessary for now
# NOT ANONYMISED YET
# ll_spatial <- list.files(synopsis_cache) |>
#   purrr::map_dfr(function(x) {
#     dat <- readRDS(file.path(synopsis_cache, x))
#     dat$cpue_spatial_ll
#   })



# Bathymetric data (in progress and maybe not needed)
# -----------------------------------------------------------------------------
# Load the DEM in R
# library(terra)
# library(tidyterra)

# # st_layers("~/R_DFO/gfdata/scratch/canada_west_coast_DEM_original.gdb")
# # dem <- st_read("~/R_DFO/gfdata/scratch/canada_west_coast_DEM_original.gdb", layer = "WEST_COAST_DEM") |>
# #   janitor::clean_names()
# # dem <- terra::classify(dem0, rcl = cbind(0, Inf, NA))
# # dem <- terra::mask(dem0, dem0 <= 0)
# # writeRaster(dem, here::here("data-generated/spatial/dem_excluded_land.tif"), overwrite = TRUE)
# dem <- terra::rast(here::here("data-generated/spatial/dem_excluded_land.tif"))


# ggplot() +
#   geom_sf(data = dem, aes(fill = depth))



# # benchmark test:
# mb_results <- microbenchmark(
#   classify = {
#     dem_classify <- terra::classify(dem0, rcl = cbind(0, Inf, NA))
#   },
#   mask = {
#     dem_incorrect <- terra::mask(dem0, dem0 <= 0)
#   },
#   times = 1
# )
# beepr::beep()
