#' Fit sdmTMB model to HBLL survey data
#'
#' @param dat Data frame containing survey data
#' @param survey_type Survey type ("HBLL INS" or "HBLL OUT")
#' @param species Species name
#' @param fit_dir Directory to save fitted models
#' @param mesh_cutoff Mesh cutoff distance (default: 10)
#' @param check_cache Check for cached model (default: TRUE)
#' @param formula Model formula (default: catch_count ~ 1)
#' @param ... Additional arguments passed to sdmTMB
#'
#' @return Fitted sdmTMB model object
#'
fit_hbll <- function(dat, survey_type, species, fit_dir, mesh_cutoff = 10,
                        check_cache = TRUE,
                        formula = catch_count ~ 1,
                        family = sdmTMB::nbinom2(link = "log"),
                        spatial = "on",
                        spatiotemporal = "rw",
                        extra_time = NULL,
                        offset = 'offset',
                        time = "year",
                        anisotropy = TRUE,
                        tag = NULL,
                        ...) {

  dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
  fname <- paste(c(species, survey_type, tag), collapse = "-") |>
    gsub("[^a-zA-Z0-9_.-]", "-", x = _)
  rds_file <- file.path(fit_dir, paste0(fname, ".rds"))
  hash_file <- file.path(fit_dir, paste0(fname, ".hash"))

  dat <- dat |> filter(str_detect(survey_abbrev, survey_type))

  model_state <- list(
    dat,
    mesh_cutoff,
    formula,
    family,
    spatial,
    spatiotemporal,
    time,
    offset,
    anisotropy,
    packageVersion("sdmTMB"),
    list(...)
  )
  current_hash <- digest::digest(model_state)

  if (check_cache && file.exists(rds_file) && file.exists(hash_file)) {
    cached_hash <- readLines(hash_file, warn = FALSE)
    if (identical(cached_hash, current_hash)) {
      message("Cache hit. Loading model from: ", rds_file)
      return(readRDS(rds_file))
    }
  }

  message(
    "Cache missing or invalid. Fitting model for: ", fname,
    "\nCache file will be saved to: ", rds_file
  )

  mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)

  if (any(spatiotemporal != "off")) {
    message("Spatiotemporal = ", spatiotemporal, ". Finding missing years.")
    extra_time <- sdmTMB:::find_missing_time(dat$year)
  }

  fit <- sdmTMB::sdmTMB(
    formula = formula,
    data = dat,
    mesh = mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = "on",
    spatiotemporal = "rw",
    extra_time = extra_time,
    offset = 'offset',
    time = "year",
    anisotropy = TRUE,
    ...
  )

  saveRDS(fit, file = rds_file)
  writeLines(current_hash, con = hash_file)

  return(fit)
}

#' Predict from fitted sdmTMB model
#'
#' @param fit Fitted sdmTMB model object
#' @param grid Data frame containing prediction grid
#'
#' @return Data frame of predictions
predict_hbll <- function(fit, grid, re_form = NULL) {
  # Filter grid for survey type
  survey <- unique(fit$data$survey_abbrev)
  pred_grid <- filter(grid, survey_abbrev %in% survey)

  years <- if (all(fit$spatiotemporal != "off")) {
    sort(union(fit$fitted_time, fit$extra_time))
  } else {
    fit$fitted_time
  }

  message("Predicting on ", paste(survey, collapse = ", "), " grid ",
    ifelse(length(years) == length(c(fit$fitted_time, fit$extra_time)), "with", "without"), " extra time")

  nd <- sdmTMB::replicate_df(
    dat = pred_grid,
    time_name = "year",
    time_values = years
  )

  # Make predictions
  pred <- predict(fit, newdata = nd, se_fit = FALSE, re_form = re_form)
}

#' Plot predictions from sdmTMB model
#'
#' @param pred Data frame of predictions
#' @param xlim Optional x-axis limits
#' @param ylim Optional y-axis limits
#' @param rotation Rotation angle in degrees (default: NULL)
#' @param buffer Buffer distance in meters (default: 40000)
#' @param type Type of prediction to plot ("link" or "response")
#'
#' @return ggplot object
plot_hbll_predictions <- function(pred,
                                  xlim = NULL,
                                  ylim = NULL,
                                  rotation = NULL,
                                  buffer = 40000,
                                  type = c("link", "response")) {

  type <- match.arg(type)

  # Calculate pred_value based on type
  pred <- mutate(pred, pred_value = if (type == "link") est else exp(est))

  # Convert to spatial data
  pred_sf <- pred |>
    mutate(lon = 1000 * X, lat = 1000 * Y) |>
    st_as_sf(coords = c("lon", "lat"), crs = 3156)

  bbox <- pred_sf |>
    st_buffer(dist = buffer) |>
    rotate_a(a = rotation) |>
    st_bbox(pred_sf)

  if (is.null(xlim)) xlim <- bbox[c("xmin", "xmax")]
  if (is.null(ylim)) ylim <- bbox[c("ymin", "ymax")]

  # Create plot
  ggplot() +
    geom_sf(data = pacea::bc_coast |> rotate_a(a = rotation), fill = "grey90") +
    geom_sf(data = pred_sf |> rotate_a(a = rotation), aes(colour = pred_value)) +
    viridis::scale_colour_viridis(option = "plasma") +
    theme_light() +
    coord_sf(xlim = xlim, ylim = ylim) +
    labs(colour = if (type == "link") "est" else "exp(est)")
}

#' Convert coordinates to longitude/latitude
#'
#' @param data Data frame containing coordinate columns
#' @param x_col Name of X coordinate column (default: "X")
#' @param y_col Name of Y coordinate column (default: "Y")
#' @param mult Multiplier for coordinates (default: 1000)
#' @param to_sf Convert to sf object (default: TRUE)
#' @param crs_from Source coordinate reference system (default: 3156)
#' @param crs_to Target coordinate reference system (default: 4326)
#'
#' @return Data frame with lon/lat columns, optionally as sf object
#' @export
make_lonlat <- function(data, x_col = "X", y_col = "Y",
                         mult = 1000,
                         to_sf = TRUE,
                         crs_from = 3156, crs_to = 4326) {
  df <- data |>
    dplyr::mutate(
      lon = .data[[x_col]] * mult,
      lat = .data[[y_col]] * mult
    )

  if (to_sf) {
    df <- df |>
      sf::st_as_sf(coords = c("lon", "lat"), crs = crs_from) |>
      sf::st_transform(crs = crs_to)
  }
}

#' Rotate spatial features while maintaining appropriate coordinate reference system
#'
#' @param sf_obj An sf object to rotate
#' @param a Angle in degrees to rotate (default: 90)
#'
#' @return Rotated sf object in oblique Mercator projection
#' @export
rotate_a <- function(sf_obj, a = 90){
  rotated_crs <- paste0("+proj=omerc +lat_0=0 +lonc=-9 +gamma=", -a)
  sf_obj <- sf_obj |> st_transform(rotated_crs)
}

# Fix me: - edit this documentation and double check how the centroids and coords are calculated
  # - add a function to convert from km to m?
  # - double check if the polygons for HBLL account for coastline - I think they do not...
  # this means that the points might fall on land - this is problematic if using
  # to get depths or matching with other oceanographic data.

#' Load survey block data (polygon, point, or coordinate table)
#'
#' Returns the built-in `survey_blocks` dataset as either polygons, centroids, or
#' a data frame with coordinates (X, Y) in kilometres.
#' Uses `sf::st_points_on_surface()` to use points that fall within the polygon,
#' rather than using `sf::st_centroid()` which could lead to points that fall
#' outside of irregularly shaped polygons.
#'
#' @param type Character string specifying the output format. One of:
#'   - `"polygon"` (default): returns an `sf` object with polygon geometries.
#'   - `"centroid"`: returns an `sf` object with the centroid point for each block.
#'   - `"coords"`: returns a `tibble` with columns `X` and `Y` (in kilometres),
#'     representing point-on-surface coordinates extracted from each polygon.
#'
#' @return Either an `sf` object or a `tbl` depending on `type`.
#' @export
#'
#' @examples
#' survey_blocks_data("polygon") |> plot()
#' survey_blocks_data("centroid")
#' survey_blocks_data("coords")
load_survey_blocks <- function(type = c("polygon", "centroid", "coords")) {
  type <- match.arg(type)
  dat <- gfdata::survey_blocks

  if (type == "centroid") {
    return(sf::st_point_on_surface(dat)) # sf points
  }

  if (type == "coords") {
    pts <- sf::st_point_on_surface(dat)
    coords <- sf::st_coordinates(pts) / 1000  # convert metres to km
    df <- sf::st_drop_geometry(pts)
    df$X <- coords[, 1]
    df$Y <- coords[, 2]
    df <- dplyr::as_tibble(df)
    return(df)
  }

  return(dat)  # default: polygon sf
}
