# Bathymetric data (in progress and maybe not needed)
# -----------------------------------------------------------------------------
# Load the DEM in R
library(terra)
library(tidyterra)
library(sf)
library(dplyr)

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-utils.R"))
source(here::here("R", "00-overlay-functions.R"))

strata <- readRDS(here::here("data-raw", "grouping-table.rds")) |>
  select(grouping_code, min_depth_m, max_depth_m, depth_operator, strata_depth_label) |>
  mutate(max_depth_m = as.numeric(max_depth_m))

if (!file.exists(here::here("data-generated/spatial/dem_with_land_mask_positive_depths.tif"))) {
  dem0 <- terra::rast("~/R_DFO/gfdata/scratch/canada_west_coast_DEM_original.gdb",
  lyrs = "WEST_COAST_DEM") |>
  janitor::clean_names()

  dem <- terra::ifel(dem0 >= 0, NA, dem0)
  dem <- abs(dem)

  terra::writeRaster(dem, here::here("data-generated/spatial/dem_with_land_mask_positive_depths.tif"), overwrite = TRUE)
} else {
  dem <- terra::rast(here::here("data-generated/spatial/dem_with_land_mask_positive_depths.tif"))
}
meep()

library(tmap)
tmap_mode("view")
# test <- terra::ifel(dem == 99999, NA, dem)

test2 <- terra::global(dem, fun = "isNA")
meep()

tm_shape(dem) +
   tm_raster(
    col.scale = tm_scale_continuous(
      value.na = "grey75",
      limits = c(0, 300)      # or min(dem, na.rm = TRUE) if not starting at 0
    )
  )

# HBLL grid polygons ----
hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N", "HBLL INS S"))
# Transform HBLL grid to match DEM CRS
hbll_grid_poly_transformed <- hbll_grid_poly |>
  st_transform(crs = terra::crs(dem)) |>
  left_join(strata, by = "grouping_code")

# Get DEM depth data for HBLL survey grid
# ----------------------------------------
# Below chunk is unnecessary and I think led to too many NA values.
# Step 1: Filter for HBLL survey and crop DEM
# Save the cropped DEM
# if (!file.exists(here::here("data-generated/spatial/dem_hbll_cropped.tif"))) {
#   # Crop DEM to HBLL OUT N extent
#   dem_cropped <- terra::crop(dem, hbll_grid_poly_transformed)
#   meep()
#   terra::writeRaster(dem_cropped, here::here("data-generated/spatial/dem_hbll_cropped.tif"), overwrite = TRUE)
#   # Step 2: Aggregate DEM to 100m resolution
#   dem_100m <- terra::aggregate(dem_cropped, fact = 10, fun = "mean", na.rm = TRUE)
#   terra::writeRaster(dem_100m, here::here("data-generated/spatial/dem_hbll_cropped_100m.tif"), overwrite = TRUE)
# } else {
#   dem_cropped <- terra::rast(here::here("data-generated/spatial/dem_hbll_cropped.tif"))
#   dem_100m <- terra::rast(here::here("data-generated/spatial/dem_hbll_cropped_100m.tif"))
# }

# # get 100 m DEM for full coast
# if (!file.exists(here::here("data-generated/spatial/dem_full_coast_100m.tif"))) {
#   # Aggregate DEM to 100m resolution
#   dem_100m <- terra::aggregate(dem, fact = 10, fun = "mean", na.rm = TRUE)
#   terra::writeRaster(dem_100m, here::here("data-generated/spatial/dem_full_coast_100m.tif"), overwrite = TRUE)
# } else {
#   dem_100m <- terra::rast(here::here("data-generated/spatial/dem_full_coast_100m.tif"))
# }
# dem_100m <- abs(dem_100m) # go back and apply to excluded land raster

library(tmap)
tmap_mode("view")
tm_shape(dem) +
   tm_raster(
    col.scale = tm_scale_continuous(
      value.na = "grey75",
      limits = c(0, 300)      # or min(dem, na.rm = TRUE) if not starting at 0
    )
  ) +
  tm_shape(hbll_grid_poly_transformed) +
  tm_borders(
    col = "strata_depth_label"
  ) +
  tm_polygons(
    col_alpha = 0,
    fill_alpha = 0,
    popup.vars = c("block_id", "strata_depth_label")
  )


# Calculate depth summaries for each HBLL OUT N grid cell
grid_depths <- hbll_grid_poly_transformed |>
  st_drop_geometry() |>
  select(block_id, grouping_code, min_depth_m, max_depth_m, depth_operator)

# Get all depth values from the DEM that overlap with HBLL grid cells
if (!file.exists(here::here("data-generated/spatial/dem_hbll_extracts.rds"))) {
  all_extracts <- terra::extract(dem, hbll_grid_poly_transformed, df = TRUE) |>
    rename(dem_depth = west_coast_dem)
  all_extracts$dem_depth <- abs(all_extracts$dem_depth)
  saveRDS(all_extracts, here::here("data-generated/spatial/dem_hbll_extracts.rds"))
  meep()
} else {
  all_extracts <- readRDS(here::here("data-generated/spatial/dem_hbll_extracts.rds"))
}


library(data.table)

# Option 1.
# Version of depth summaries without filtering to depths only within the designated depth stratum:
if (!file.exists(here::here("data-generated", "hbll-dem-grid-depths.rds"))) {
test_polys <- hbll_grid_poly_transformed
tictoc::tic("depth_summaries")
# This takes 1.3 hours to run...
depth_summaries <- data.frame(
  block_id         = test_polys$block_id,
  survey_series_id = test_polys$survey_series_id,
  depth_dem_mean   = terra::extract(dem, test_polys, fun = mean, na.rm = TRUE)[, 2],
  depth_dem_sd     = terra::extract(dem, test_polys, fun = sd,   na.rm = TRUE)[, 2],
  depth_dem_min    = terra::extract(dem, test_polys, fun = min,  na.rm = TRUE)[, 2],
  depth_dem_max    = terra::extract(dem, test_polys, fun = max,  na.rm = TRUE)[, 2],
  n_pixels         = terra::extract(dem, test_polys, fun = function(x) sum(!is.na(x)))[, 2] # non-na pixel count
)
tictoc::toc()
meep()

  saveRDS(depth_summaries, here::here("data-generated", "hbll-dem-grid-depths.rds"))
} else {
  depth_summaries <- readRDS(here::here("data-generated", "hbll-dem-grid-depths.rds"))
}


dem_centroids <- hbll_grid_poly_transformed |> st_centroid() |> st_transform(crs = terra::crs(dem))
dem_centroids$depth_centroid <- terra::extract(dem, terra::vect(centroids))[, 2]
glimpse(dem_centroids)
saveRDS(dem_centroids, here::here("data-generated", "hbll-dem-grid-centroid-depths.rds"))

# Option 2.
# Code to filter to depths only within the designated block stratum
# dem_hbll0 <- all_extracts |>
#   left_join(
#     hbll_grid_poly_transformed |>
#       sf::st_drop_geometry() |>
#       mutate(ID = row_number()) |>
#       select(ID, survey_series_id, block_id, depth_m, min_depth_m, max_depth_m, depth_operator),
#     by = "ID"
#   )
# meep()

# op_levels <- c(
#   ">= MIN_DEPTH_M and < MAX_DEPTH_M",
#   ">= MIN_DEPTH_M and <= MAX_DEPTH_M",
#   "> MIN_DEPTH_M and < MAX_DEPTH_M",
#   "> MIN_DEPTH_M and <= MAX_DEPTH_M",
#   ">= MIN_DEPTH_M",
#   "> MIN_DEPTH_M"
# )

# dt <- as.data.table(dem_hbll0)
# dt[, op_code := match(depth_operator, op_levels)]

# dem_hbll <- dt[
#   !is.na(dem_depth) & (
#     is.na(min_depth_m) | is.na(max_depth_m) |
#     (op_code == 1L & dem_depth >= min_depth_m & dem_depth < max_depth_m) |
#     (op_code == 2L & dem_depth >= min_depth_m & dem_depth <= max_depth_m) |
#     (op_code == 3L & dem_depth > min_depth_m & dem_depth < max_depth_m) |
#     (op_code == 4L & dem_depth > min_depth_m & dem_depth <= max_depth_m) |
#     (op_code == 5L & dem_depth >= min_depth_m) |
#     (op_code == 6L & dem_depth > min_depth_m) |
#     is.na(op_code)
#   )
# ]

# glimpse(dem_hbll)
# meep()

# I think this code block would need to be updated to match the above depth_summaries before using
# depth_summaries <- #dem_hbll[, .(
#   dem_hbll0[, .(
#     depth_dem_mean = mean(dem_depth),
#     depth_dem_sd   = sd(dem_depth),
#     depth_dem_min  = min(dem_depth),
#     depth_dem_max  = max(dem_depth),
#     n_pixels       = .N
#   ), by = .(block_id, survey_series_id)]
# saveRDS(depth_summaries, here::here("data-generated/spatial/hbll_depth_summaries.rds"))