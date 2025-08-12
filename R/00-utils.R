#' Get Model Years
#'
#' Extracts years from an sdmTMB model object.
#'
#' @param fit A sdmTMB model object containing spatiotemporal information and time data.
#' @return A sorted numeric vector of years. For spatiotemporal models, includes all years.
#'   For non-spatiotemporal models, excludes extra time years.
#'
#' @examples
#' years <- get_model_years(fit)
#'
get_model_years <- function(fit) {
  if (any(fit$spatiotemporal != "off")) {
    sort(fit$time_lu$time_from_data)
  } else {
    sort(fit$time_lu$time_from_data[!fit$time_lu$extra_time])
  }
}

#' Get model parameters
#' @param fit Fitted model object
#' @return Combined fixed and random parameters
get_model_pars <- function(fit) {
  bind_rows(tidy(fit), tidy(fit, "ran_pars")) |>
  as.data.frame()
}

#' Get marginal standard deviation of the spatiotemporal random field
#' @param fit Fitted model object
#' @return Marginal standard deviation of the spatiotemporal random field
get_marginal_sigma_E <- function(fit) {
  get_model_pars(fit) |>
    filter(term == "sigma_E") |> pull(estimate)
}

#' Convert sdmTMB XY (km) to sf object
#'
#' @param data Data frame containing coordinate columns
#' @param x_col Name of X coordinate column (default: "X")
#' @param y_col Name of Y coordinate column (default: "Y")
#' @param mult Multiplier for coordinates (default: 1000, unless crs_from is 4326 then 1)
#' @param crs_from Source coordinate reference system (default: 3156)
#' @param crs_to Target coordinate reference system (default: 4326)
#'
#' @return sf object, defaults to WGS84 (EPSG:4326)
#' @export
XY_to_sf <- function(data, x_col = "X", y_col = "Y",
                     mult = 1000,
                     crs_from = 3156, crs_to = 4326) {
  if (crs_from == 4326) mult <- 1

  df <- data |>
    dplyr::mutate(
      x = .data[[x_col]] * mult,
      y = .data[[y_col]] * mult
    )

  df <- df |>
    sf::st_as_sf(coords = c("x", "y"), crs = crs_from) |>
    sf::st_transform(crs = crs_to)
  df
}

#' Get plot limits for sf object
#'
#' @param sf_obj An sf object
#' @param xlim x-axis limits (default: NULL)
#' @param ylim y-axis limits (default: NULL)
#' @param buffer Buffer distance in m (default: 1000)
#'
get_plot_limits <- function(sf_obj, xlim = NULL, ylim = NULL, buffer = 1000, crs_out = 4326) {
  stopifnot(inherits(sf_obj, "sf"))

  if (!is.null(buffer)) {

    sf_obj <- sf_obj |>
      sf::st_transform(crs = 3156) |>
      sf::st_buffer(dist = buffer)
  }

  bbox <- sf::st_bbox(sf_obj)

  if (is.null(xlim)) xlim <- bbox[c("xmin", "xmax")]
  if (is.null(ylim)) ylim <- bbox[c("ymin", "ymax")]

  coord_sf(xlim = xlim, ylim = ylim, crs = st_crs(sf_obj))
}

#' Rotate spatial features while maintaining appropriate coordinate reference system
#'
#' @param sf_obj An sf object to rotate
#' @param a Angle in degrees to rotate (default: 90)
#'
#' @return Rotated sf object in oblique Mercator projection
#' @export
rotate_a <- function(sf_obj, a = 90){
  if (is.null(a)) return(sf_obj)

  rotated_crs <- paste0("+proj=omerc +lat_0=0 +lonc=-9 +gamma=", -a)
  sf_obj <- sf_obj |> st_transform(rotated_crs)
}

#' Beep only for a specific user.
#'
#' @param target_user The username to check for.
#' @param ... Additional arguments to pass to beepr::beep(), e.g., 'sound', 'expr'.
#'
meep <- function(user = "jilliandunic", ...) {
  current_user <- Sys.info()['user']

  if (current_user == user) {
    beepr::beep(...)
  }
}

sp_to_hyphens <- function(sp) {
  sp |>
    gsub(" ", "-", x = _) |>
    gsub("/", "-", x = _)
}

sp_from_hyphens <- function(sp) {
  if (grepl("rougheye-blackspotted", sp)) sp <- gsub("rougheye-blackspotted", "rougheye/blackspotted", sp)

  sp |> gsub("-", " ", x = _)
}
