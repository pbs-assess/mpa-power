# Extract DEM depths for historical HBLL survey transects
# Saves: data-generated/spatial/hbll_survey_dem_depths.rds
# -----------------------------------------------------------------------------
library(terra)
library(sf)
library(dplyr)

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-utils.R"))

# Load raw DEM (positive depths, land masked to NA) ----
dem <- terra::rast(here::here("data-generated/spatial/dem_with_land_mask_positive_depths.tif"))

# Load survey set data ----
sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- sp_dat0 |>
  filter(survey_series_id.x %in% hbll_ssids) |>
  select(survey_abbrev, fishing_event_id, longitude, latitude, longitude_end, latitude_end, depth_m)

# DEM extraction functions ----
extract_midpoint_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster) {
  mid_lons <- (start_lon + end_lon) / 2
  mid_lats <- (start_lat + end_lat) / 2
  all_points <- terra::vect(data.frame(x = mid_lons, y = mid_lats),
    geom = c("x", "y"), crs = "EPSG:4326")
  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  terra::extract(dem_raster, points_trans)[, 2]
}

extract_endpoint_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster) {
  all_lons <- c(start_lon, end_lon)
  all_lats <- c(start_lat, end_lat)
  all_points <- terra::vect(data.frame(x = all_lons, y = all_lats),
    geom = c("x", "y"), crs = "EPSG:4326")
  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  depths <- terra::extract(dem_raster, points_trans)[, 2]
  n <- length(start_lon)
  (depths[1:n] + depths[(n + 1):(2 * n)]) / 2
}

extract_line_depth_vectorized <- function(
  start_lon, start_lat, end_lon, end_lat, dem_raster, n_points = 5) {
  n_lines <- length(start_lon)
  all_lons <- numeric(n_lines * n_points)
  all_lats <- numeric(n_lines * n_points)
  line_ids <- numeric(n_lines * n_points)

  for (i in seq_len(n_lines)) {
    idx <- ((i - 1) * n_points + 1):(i * n_points)
    all_lons[idx] <- seq(start_lon[i], end_lon[i], length.out = n_points)
    all_lats[idx] <- seq(start_lat[i], end_lat[i], length.out = n_points)
    line_ids[idx] <- i
  }
  all_points <- terra::vect(data.frame(x = all_lons, y = all_lats),
    geom = c("x", "y"), crs = "EPSG:4326")
  points_trans <- terra::project(all_points, terra::crs(dem_raster))
  all_depths <- terra::extract(dem_raster, points_trans)[, 2]

  line_means <- tapply(all_depths, line_ids, mean, na.rm = TRUE)
  line_ranges <- tapply(all_depths, line_ids, function(x) {
    if (all(is.na(x))) return(NA)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  start_indices <- seq(1, n_lines * n_points, by = n_points)
  end_indices   <- seq(n_points, n_lines * n_points, by = n_points)

  list(
    mean  = as.numeric(line_means),
    range = as.numeric(line_ranges),
    start = all_depths[start_indices],
    end   = all_depths[end_indices]
  )
}

# Extract depths ----
line_depths <- extract_line_depth_vectorized(
  sp_dat$longitude, sp_dat$latitude,
  sp_dat$longitude_end, sp_dat$latitude_end,
  dem
)

sp_dat_depths <- sp_dat |>
  mutate(
    depth_midpoint   = extract_midpoint_depth_vectorized(
      longitude, latitude, longitude_end, latitude_end, dem),
    depth_line_mean  = line_depths$mean,
    depth_line_range = line_depths$range,
    depth_line_start = line_depths$start,
    depth_line_end   = line_depths$end,
    depth_endpoint   = extract_endpoint_depth_vectorized(
      longitude, latitude, longitude_end, latitude_end, dem)
  ) |>
  mutate(across(starts_with("depth_"), ~ round(., 1)))
meep()
saveRDS(sp_dat_depths, here::here("data-generated/spatial/hbll_survey_dem_depths.rds"))
