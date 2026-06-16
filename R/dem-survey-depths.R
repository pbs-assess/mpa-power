# Extract DEM depths for historical HBLL survey transects
# Saves: data-generated/spatial/hbll_survey_dem_depths.rds
# -----------------------------------------------------------------------------
library(terra)
library(sf)

source(here::here("R", "00-utils.R"))
set_synopsis_cache()

# Load raw DEM (positive depths, land masked to NA) ----
dem <- terra::rast(here::here("data-generated/spatial/dem_with_land_mask_positive_depths.tif"))

# Load survey set data ----
sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- sp_dat0 |>
  filter(survey_series_id.x %in% hbll_ssids) |>
  select(survey_abbrev, fishing_event_id, longitude, latitude, longitude_end, latitude_end, depth_m) |>
  distinct(survey_abbrev, fishing_event_id, .keep_all = TRUE) |> # there are duplicate fishing events with same data except for start/retrieval time
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N", "HBLL INS S"))

# DEM extraction functions ----
extract_line_depth_vectorized <- function(start_lon, start_lat, end_lon, end_lat, dem_raster) {
  n <- length(start_lon)

  extract_pts <- function(lons, lats) {
    terra::vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326") |>
      terra::project(terra::crs(dem_raster)) |>
      (\(pts) terra::extract(dem_raster, pts)[, 2])()
  }

  lines_vect <- purrr::pmap(
    list(start_lon, start_lat, end_lon, end_lat),
    ~ sf::st_linestring(matrix(c(..1, ..2, ..3, ..4), ncol = 2, byrow = TRUE))
  ) |>
    sf::st_sfc(crs = 4326) |>
    terra::vect() |>
    terra::project(terra::crs(dem_raster))

  extracted   <- terra::extract(dem_raster, lines_vect)
  depth_col   <- names(extracted)[2]
  line_means  <- tapply(extracted[[depth_col]], extracted$ID, mean, na.rm = TRUE)
  line_ranges <- tapply(extracted[[depth_col]], extracted$ID, function(x) {
    if (all(is.na(x))) return(NA)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })

  means_out  <- rep(NA_real_, n)
  ranges_out <- rep(NA_real_, n)
  means_out[as.integer(names(line_means))]   <- as.numeric(line_means)
  ranges_out[as.integer(names(line_ranges))] <- as.numeric(line_ranges)

  list(
    mean     = means_out,
    range    = ranges_out,
    start    = extract_pts(start_lon, start_lat),
    midpoint = extract_pts((start_lon + end_lon) / 2, (start_lat + end_lat) / 2),
    end      = extract_pts(end_lon, end_lat)
  )
}

# Extract depths ----
if (!file.exists(here::here("data-generated", "hbll-dem-survey-depths.rds"))) {
  line_depths <- extract_line_depth_vectorized(
    sp_dat$longitude, sp_dat$latitude,
    sp_dat$longitude_end, sp_dat$latitude_end,
    dem
  )
  sp_dat_depths <- sp_dat |>
    mutate(
      depth_line_mean  = line_depths$mean, # this probably also doesn't make sense without the track data
      depth_line_range = line_depths$range, # this probably also doesn't make sense without the track data
      depth_start      = line_depths$start,
      depth_midpoint   = line_depths$midpoint, # this may not make sense since the midpoint could fall on land; transects are not always a straight line
      depth_end        = line_depths$end
    ) |>
    mutate(across(starts_with("depth_"), ~ round(., 1)))
  meep()
  saveRDS(sp_dat_depths, here::here("data-generated", "hbll-dem-survey-depths.rds"))
} else {
  sp_dat_depths <- readRDS(here::here("data-generated", "hbll-dem-survey-depths.rds"))
}