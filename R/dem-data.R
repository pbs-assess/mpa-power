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

# Yelloweye rockfish data ----
ye_dat <- purrr::map_df(
  list.files(here::here("data-generated", "cleaned-species-data"), pattern = "yelloweye-rockfish",
  full.names = TRUE), function(x) {
    readRDS(x)
})

restricted_df <- readRDS(here::here("data-generated", "hbll-restricted-sf.rds"))
ye_dat0 <- readRDS(file.path(synopsis_cache, paste0("yelloweye-rockfish", ".rds")))$survey_sets
ye_dat <- ye_dat0 |> filter(survey_series_id.x %in% hbll_ssids) |>
  filter(survey_abbrev != "HBLL INS S") # keep dataset consistent with what we have used
  # but in the future it would be good to realise that it will be better to
  # filter out sampling locations that do not overlap with the grid rather than
  # by survey name since that can changed based on extent sampled in a given year
  # Something that could help would be to include block_id in the survey_sets

ye_dat_sf <- ye_dat |>
  XY_to_sf(coords = c("longitude", "latitude"), crs_from = 4326, crs_to = 4326)

# Historical transects ---
ye_line_geoms <- mapply(
  create_linestring,
  ye_dat$longitude, ye_dat$latitude,
  ye_dat$longitude_end, ye_dat$latitude_end,
  SIMPLIFY = FALSE
)

hbll_transects <- st_sfc(ye_line_geoms, crs = 4326) %>%
  st_sf(ye_dat, geometry = .) |>
  st_transform(crs = st_crs(hbll_grid_poly_transformed))

hbll_grid_sample_lu <- st_join(hbll_transects, hbll_grid_poly_transformed,
  join = st_intersects, suffix = c("_sampled", "_grid"))

# hbll_grid_sample_lu |>
#   filter(is.na(block_id)) |>
#   glimpse()


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
if (!file.exists(here::here("data-generated/spatial/hbll_depth_summaries.rds"))) {
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

  saveRDS(depth_summaries, here::here("data-generated/spatial/hbll_depth_summaries.rds"))
} else {
  depth_summaries <- readRDS(here::here("data-generated/spatial/hbll_depth_summaries.rds"))
}

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


depth_summaries_sf <- hbll_grid_poly_transformed |>
  left_join(depth_summaries, by = c("block_id", "survey_series_id"))

tmap_mode("view")
  tm_shape(depth_summaries_sf) +
    tm_polygons(
      fill = "strata_depth_label",
      fill_alpha = 0.1,
      popup.vars = c("block_id", "strata_depth_label", "depth_dem_mean", "depth_dem_sd", "n_pixels", "depth_dem_min", "depth_dem_max")
    )

ssid_labels <- c("HBLL OUT N" = 22, "HBLL OUT S" = 36, "HBLL INS N" = 39, "HBLL INS S" = 40)

depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(dem_minus_og = depth_dem_mean - depth_m) |>
  ggplot() +
  aes(depth_m, depth_dem_mean, colour = factor(survey_series_id, levels = ssid_labels, labels = names(ssid_labels))) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Dark2") +
  # scale_y_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  labs(x = "Grid depth_m", y = "DEM mean depth", colour = "Survey") +
  facet_wrap(~ survey_series_id, ncol = 2)

depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(
    in_stratum = case_when(
      is.na(min_depth_m) | is.na(max_depth_m) ~ NA,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M and < MAX_DEPTH_M")  ~ depth_dem_mean >= min_depth_m & depth_dem_mean < max_depth_m,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M and <= MAX_DEPTH_M") ~ depth_dem_mean >= min_depth_m & depth_dem_mean <= max_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M and < MAX_DEPTH_M")   ~ depth_dem_mean > min_depth_m  & depth_dem_mean < max_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M and <= MAX_DEPTH_M")  ~ depth_dem_mean > min_depth_m  & depth_dem_mean <= max_depth_m,
      stringr::str_detect(depth_operator, ">= MIN_DEPTH_M") ~ depth_dem_mean >= min_depth_m,
      stringr::str_detect(depth_operator, "> MIN_DEPTH_M")  ~ depth_dem_mean > min_depth_m,
      TRUE ~ NA
    )
  ) |>
  ggplot() +
  aes(depth_m, depth_dem_mean, colour = in_stratum) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
    na.value = "grey60",
    labels = c("TRUE" = "In stratum", "FALSE" = "Out of stratum", "NA" = "No bounds")
  ) +
  # scale_x_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  # scale_y_continuous(trans = ggsidekick::fourth_root_power_trans()) +
  labs(x = "Grid depth_m", y = "DEM mean depth", colour = NULL) +
  facet_wrap(~ survey_series_id, ncol = 2)

strata_rects <- distinct(hbll_grid_poly_transformed, survey_abbrev, survey_series_id, strata_depth_label) |>
  tidyr::separate(strata_depth_label, into = c("depth_min", "depth_max"),
    sep = " - ", convert = TRUE, remove = FALSE)

depth_summaries_sf |>
  st_drop_geometry() |>
  mutate(dem_minus_og = depth_dem_mean - depth_m) |>
  ggplot() +
    # aes(x = depth_m) +
    aes(x = depth_dem_mean) +
    geom_histogram() +
    geom_rect(
      data = strata_rects,
      aes(xmin = depth_min, xmax = depth_max, ymin = -Inf, ymax = Inf, fill = strata_depth_label), alpha = 0.3, colour = NA, inherit.aes = FALSE
    ) +
    facet_wrap(~ survey_series_id, ncol = 2)



# DEM depth data from survey sets -----
# use start of transect for now?
# use 100 m aggregate?
test <- terra::extract(dem_100m, ye_dat, df = TRUE) |>

# -------
# Visualise cells and their appropriate depth strata
poly_vect <- terra::vect(hbll_grid_poly_transformed)

# Rasterize the stratum depth bounds onto the DEM grid
min_depth_rast <- terra::rasterize(poly_vect, dem_100m, field = "min_depth_m")
max_depth_rast <- terra::rasterize(poly_vect, dem_100m, field = "max_depth_m")
meep()

# Keep only pixels within their block's valid depth range
dem_filtered <- terra::mask(
  dem_100m,
  dem_100m >= min_depth_rast & dem_100m < max_depth_rast,
  maskvalue = FALSE
)

tmap_mode("view")
tm_shape(dem_filtered) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0, 500),
      values = "Blues"
    )
  ) +
  tm_shape(hbll_grid_poly_transformed) +
  tm_polygons(
      fill = "strata_depth_label",
      fill_alpha = 0.1,
      popup.vars = c("block_id", "strata_depth_label", "depth_dem_mean", "depth_dem_sd", "n_pixels", "min_depth_m", "max_depth_m")
    )
  # tm_borders(col = "black", lwd = 0.5)
# -------

# Join with grid cell metadata
# dem_hbll <-
test <- all_extracts |>
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




library(DT)

# View any data frame in browser
hbll_grid_sample_lu |>
  filter(is.na(block_id)) |>
datatable() |> htmltools::browsable()

tmap_mode("view")
tm_shape(hbll_grid_sample_lu |>
  filter(is.na(block_id))) +
  tm_lines(col = "red") +
  tm_shape(hbll_grid_poly_transformed) +
  tm_polygons(
    fill_alpha = 0.1,
    col = "grey30",
    lwd = 0.4,
    popup.vars = c("block_id", "grouping_code", "strata_depth_label", "depth_m")
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
# saveRDS(depth_summaries, here::here("data-generated/spatial/hbll_depth_summaries.rds"))


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
