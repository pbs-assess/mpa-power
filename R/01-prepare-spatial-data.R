library(dplyr)
library(ggplot2)
library(sf)
library(stringr)

source(here::here("R", "00-setup.R"))

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
# st_layers(here::here("data-raw", "spatial", "Public", "Spatial2025_Public_Network_footprint.gdb"))
public_mpa <- st_read(here::here("data-raw", "spatial", "Public",
  "Spatial2025_Public_Network_footprint.gdb"), layer = "Spatial2025_Public_Network_footprint") |>
  janitor::clean_names()

simple_mpa <- public_mpa |>
  select(uid, map_label, common_site_name_site_profile, category_simple, name_2025) |>
    st_cast("MULTIPOLYGON")  # Convert MULTISURFACE
saveRDS(simple_mpa, file.path("data-generated", "spatial", "simple-mpa.rds"))

mpa_500 <- simple_mpa |>
  st_simplify(dTolerance = 500)
saveRDS(mpa_500, file.path("data-generated", "spatial", "simple-mpa-500m.rds"))


comm_ll_activity_status <- public_mpa |>
  select(uid, hu_commercial_harvest_bottom_longline_demersal_hookand_line,
    category_detailed, category_simple) |>
  left_join(activity_status_lu,
    by = c("hu_commercial_harvest_bottom_longline_demersal_hookand_line" = "activity_status")) |>
  sf::st_cast("MULTIPOLYGON") %>%
  mutate(mpa_id = row_number(),
         mpa_area = st_area(.))

saveRDS(comm_ll_activity_status, file.path("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))

stop()
# Human use layers
# human_layers <- st_layers(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"))
# human_layers$name

# hu_ll <- st_read(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"),
#   layer = "hu_co_demersalfishing_bottomlongline_d")


# -----------------------------------------------------------------------------
# Commercial longline data - raw
# -----------------------------------------------------------------------------
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
library(terra)
library(tidyterra)

strata <- readRDS(here::here("data-raw", "grouping-table.rds")) |>
  select(grouping_code, min_depth_m, max_depth_m, depth_operator, strata_depth_label) |>
  mutate(max_depth_m = as.numeric(max_depth_m))

# st_layers("~/R_DFO/gfdata/scratch/canada_west_coast_DEM_original.gdb")
# dem <- st_read("~/R_DFO/gfdata/scratch/canada_west_coast_DEM_original.gdb", layer = "WEST_COAST_DEM") |>
#   janitor::clean_names()
# dem <- terra::classify(dem0, rcl = cbind(0, Inf, NA))
# dem <- terra::mask(dem0, dem0 <= 0)
# writeRaster(dem, here::here("data-generated/spatial/dem_excluded_land.tif"), overwrite = TRUE)
dem <- terra::rast(here::here("data-generated/spatial/dem_excluded_land.tif"))

hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(survey_series_id %in% hbll_ssids)
# Transform HBLL grid to match DEM CRS
hbll_grid_poly_transformed <- hbll_grid_poly |>
  st_transform(crs = terra::crs(dem)) |>
  left_join(strata, by = "grouping_code")

# Step 1: Filter for HBLL survey and crop DEM
# Save the cropped DEM
if (!file.exists(here::here("data-generated/spatial/dem_hbll_cropped.tif"))) {
  # Crop DEM to HBLL OUT N extent
  dem_cropped <- terra::crop(dem, hbll_grid_poly_transformed)
  meep()
  terra::writeRaster(dem_cropped, here::here("data-generated/spatial/dem_hbll_cropped.tif"), overwrite = TRUE)
  # Step 2: Aggregate DEM to 100m resolution
  dem_100m <- terra::aggregate(dem_cropped, fact = 10, fun = "mean", na.rm = TRUE)
  terra::writeRaster(dem_100m, here::here("data-generated/spatial/dem_hbll_cropped_100m.tif"), overwrite = TRUE)
} else {
  dem_cropped <- terra::rast(here::here("data-generated/spatial/dem_hbll_cropped.tif"))
  dem_100m <- terra::rast(here::here("data-generated/spatial/dem_hbll_cropped_100m.tif"))
}

# get 100 m DEM for full coast
if (!file.exists(here::here("data-generated/spatial/dem_full_coast_100m.tif"))) {
  # Aggregate DEM to 100m resolution
  dem_100m <- terra::aggregate(dem, fact = 10, fun = "mean", na.rm = TRUE)
  terra::writeRaster(dem_100m, here::here("data-generated/spatial/dem_full_coast_100m.tif"), overwrite = TRUE)
} else {
  dem_100m <- terra::rast(here::here("data-generated/spatial/dem_full_coast_100m.tif"))
}

dem_100m <- abs(dem_100m) # go back and apply to excluded land raster

# Step 3: Calculate depth summaries for each HBLL OUT N grid cell
grid_depths <- hbll_grid_poly_transformed |>
  st_drop_geometry() |>
  select(block_id, grouping_code, min_depth_m, max_depth_m, depth_operator)

# Get all depth values from the DEM that overlap with HBLL grid cells
all_extracts <- terra::extract(dem, hbll_grid_poly_transformed, df = TRUE) |>
  rename(dem_depth = west_coast_dem)
meep()
all_extracts$dem_depth <- abs(all_extracts$dem_depth)
saveRDS(all_extracts, here::here("data-generated/spatial/dem_hbll_extracts.rds"))
all_extracts <- readRDS(here::here("data-generated/spatial/dem_hbll_extracts.rds"))

# Join with grid cell metadata
dem_hbll <- all_extracts |>
  left_join(
    hbll_grid_poly_transformed |>
      st_drop_geometry() |>
      mutate(ID = row_number()) |>
      select(ID, survey_series_id, block_id, depth_m, min_depth_m, max_depth_m, depth_operator),
    by = "ID"
  ) |>
  filter(!is.na(dem_depth)) |>
  # Apply depth filtering using specific depth operators
  filter(
    is.na(min_depth_m) | is.na(max_depth_m) |
    case_when(
      str_detect(depth_operator, ">= MIN_DEPTH_M and < MAX_DEPTH_M") ~
        dem_depth >= min_depth_m & dem_depth < max_depth_m,
      str_detect(depth_operator, ">= MIN_DEPTH_M and <= MAX_DEPTH_M") ~
        dem_depth >= min_depth_m & dem_depth <= max_depth_m,
      str_detect(depth_operator, "> MIN_DEPTH_M and < MAX_DEPTH_M") ~
        dem_depth > min_depth_m & dem_depth < max_depth_m,
      str_detect(depth_operator, "> MIN_DEPTH_M and <= MAX_DEPTH_M") ~
        dem_depth > min_depth_m & dem_depth <= max_depth_m,
      str_detect(depth_operator, ">= MIN_DEPTH_M") ~
        dem_depth >= min_depth_m,
      str_detect(depth_operator, "> MIN_DEPTH_M") ~
        dem_depth > min_depth_m,
      TRUE ~ TRUE  # Default case if operator is unclear
    )
  )


left_join(depth_summaries, dem_hbll |> distinct(ID, .keep_all = TRUE)) |>
  select(-depth_operator) %>%
  mutate(across(starts_with("depth_"), ~ round(., 1))) |>
  # Check if the dem values are better than the original strata depth_m values
  # filter(block_id %in% (hbll_out_n |> filter(depth_m < 0) |> pull(block_id))) |>
    # 1. Focus on depth comparisons
  # select(block_id, depth_m, depth_dem_mean, depth_dem_min, depth_dem_max, n_pixels) |>
  # 2. Check constraint compliance
  # select(block_id, min_depth_m, max_depth_m, depth_dem_min, depth_dem_mean, depth_dem_max, n_pixels) |>
  # 3. Assess data quality
  # select(block_id, depth_m, depth_dem_mean, depth_dem_sd, n_pixels) |>
  # Compare og values with new dem values
  mutate(dem_minus_og = round(depth_dem_mean - depth_m, 1)) |>
  arrange(desc(abs(dem_minus_og))) |>
  select(block_id, dem_minus_og,depth_m, depth_dem_mean, depth_dem_sd, n_pixels) |>
  view_dtb()

# Save depth summaries
# saveRDS(depth_summaries, here::here("data-generated/spatial/hbll_out_n_depth_summaries.rds"))
saveRDS(depth_summaries, here::here("data-generated/spatial/hbll_depth_summaries.rds"))


# Step 4: Update historical sp_dat depths for HBLL OUT N using three methods
# Load historical data (assuming this exists from your setup)
sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- sp_dat0 |>
  # filter(stringr::str_detect(survey_abbrev, "HBLL OUT N")) |>
  filter(survey_series_id.x %in% hbll_ssids) |>
  select(survey_abbrev, fishing_event_id, longitude, latitude, longitude_end, latitude_end, depth_m)

# Midpoint depth
extract_midpoint_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster) {
  # Calculate midpoints
  mid_lons <- (start_lon + end_lon) / 2
  mid_lats <- (start_lat + end_lat) / 2

  all_points <- terra::vect(data.frame(x = mid_lons, y = mid_lats),
                            geom = c("x", "y"), crs = "EPSG:4326")

  # Transform and extract
  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  terra::extract(dem_raster, points_trans)[, 2]
}

# Endpoint depth
extract_endpoint_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster) {
  # Combine start and end points
  all_lons <- c(start_lon, end_lon)
  all_lats <- c(start_lat, end_lat)
  all_points <- terra::vect(data.frame(x = all_lons, y = all_lats),
                            geom = c("x", "y"), crs = "EPSG:4326")

  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  depths <- terra::extract(dem_raster, points_trans)[, 2]

  # Average pairs: (start1 + end1) / 2, (start2 + end2) / 2
  n <- length(start_lon)
  (depths[1:n] + depths[(n + 1):(2 * n)]) / 2
}

# Line depth
extract_line_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster, n_points = 5) {
  n_lines <- length(start_lon)

  # Create all line points
  all_lons <- numeric(n_lines * n_points)
  all_lats <- numeric(n_lines * n_points)
  line_ids <- numeric(n_lines * n_points)

  for(i in 1:n_lines) {
    idx <- ((i-1) * n_points + 1):(i * n_points)
    all_lons[idx] <- seq(start_lon[i], end_lon[i], length.out = n_points)
    all_lats[idx] <- seq(start_lat[i], end_lat[i], length.out = n_points)
    line_ids[idx] <- i
  }
  # Single transform and extract
  all_points <- terra::vect(data.frame(x = all_lons, y = all_lats),
                           geom = c("x", "y"), crs = "EPSG:4326")
  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  all_depths <- terra::extract(dem_raster, points_trans)[, 2]

  # Calculate statistics for each line
  line_means <- tapply(all_depths, line_ids, mean, na.rm = TRUE)
  line_ranges <- tapply(all_depths, line_ids, function(x) {
    if(all(is.na(x))) return(NA)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })

  # Get start and end depths (first and last point of each line)
  start_indices <- seq(1, n_lines * n_points, by = n_points)
  end_indices <- seq(n_points, n_lines * n_points, by = n_points)
  start_depths <- all_depths[start_indices]
  end_depths <- all_depths[end_indices]

  list(
    mean = as.numeric(line_means),
    range = as.numeric(line_ranges),
    start = start_depths,
    end = end_depths
  )
}

# Extract line depth statistics
line_depths <- extract_line_depth_vectorized(sp_dat$longitude, sp_dat$latitude, sp_dat$longitude_end, sp_dat$latitude_end, dem_100m)

sp_dat_depths <- sp_dat |>
  mutate(
    depth_midpoint = extract_midpoint_depth_vectorized(longitude, latitude, longitude_end, latitude_end, dem_100m),
    depth_line_mean = line_depths$mean,
    depth_line_range = line_depths$range,
    depth_line_start = line_depths$start,
    depth_line_end = line_depths$end,
    depth_endpoint = extract_endpoint_depth_vectorized(longitude, latitude, longitude_end, latitude_end, dem_100m)
  ) |>
  mutate(across(starts_with("depth_"), ~ round(., 1)))

view_dtb(sp_dat_depths |> mutate(across(is.numeric, ~ round(., 1))))

sp_dat_depths |>
  pivot_longer(cols = starts_with("depth_"), names_to = "depth_method", values_to = "depth") |>
  filter(depth_method != "depth_line_range") |>
  ggplot(aes(x = depth, fill = depth_method)) +
  geom_histogram() +
  scale_fill_viridis_d() +
  # scale_x_continuous(trans = "log10") +
  facet_wrap(~ depth_method, ncol = 1, scales = "free_x") +
  theme_minimal()

sp_dat_depths$depth_line_start |> summary()

hist(sp_dat_depths$depth_line_mean - sp_dat_depths$depth_m)

p1 <- sp_dat_depths |>
  mutate(depth_diff = round(depth_line_mean - depth_m, 1)) |>
  arrange(desc(abs(depth_diff))) |>
  # view_dtb()
  slice(1:30) |>
  XY_to_sf(coords = c("longitude", "latitude"), crs_from = 4326, crs_to = 4326) |>
  ggplot() +
  geom_sf(aes(color = abs(depth_diff))) +
  theme_minimal() +
  geom_sf(data = pacea::bc_coast) +
  get_plot_limits(hbll_grid_poly_transformed)

plotly::ggplotly(p1)
plot(dem_100m )

# Summary of what I see in the tables for just HBLL OUT N:
# - line_mean might make the most sense because it reduces the chances of NA values
# - line_mean has max difference from recorded depth_m of ~116, but overall it's ok,
#   +/- 50 at least still seems to put things into the right strata depth
# - maybe I should make a better summary of diagnostics here if we add it to gfdata

sp_dat_depths |> select(survey_abbrev, fishing_event_id, depth_line_mean)
