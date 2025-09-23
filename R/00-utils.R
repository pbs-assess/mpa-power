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
#' @param coords Vector of coordinate column names (default: c("X", "Y"))
#' @param mult Multiplier for coordinates (default: 1000, unless crs_from is 4326 then 1)
#' @param crs_from Source coordinate reference system (default: 32609)
#' @param crs_to Target coordinate reference system (default: 4326)
#'
#' @return sf object, defaults to WGS84 (EPSG:4326)
#' @export
XY_to_sf <- function(data, coords = c("X", "Y"),
                     mult = 1000,
                     crs_from = 32609, crs_to = 4326) {
  if (crs_from == 4326) mult <- 1

  df <- data |>
    dplyr::mutate(
      x = .data[[coords[1]]] * mult,
      y = .data[[coords[2]]] * mult
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

# Caching functions
# Extract parameters from a fitted sdmTMB model for updating
extract_model_params <- function(model) {
  list(
    formula = formula(model),
    data = model$data,
    mesh = model$spde$mesh,
    family = family(model),
    spatial = model$spatial,
    spatiotemporal = model$spatiotemporal,
    time = model$time,
    share_range = model$share_range,
    time_varying = model$time_varying,
    time_varying_type = model$time_varying_type,
    spatial_varying = model$spatial_varying,
    weights = if(!is.null(model$weights)) digest::digest(model$weights, algo = "xxhash64") else NULL,
    offset = model$offset,
    extra_time = model$extra_time,
    reml = model$reml,
    anisotropy = model$anisotropy,
    silent = model$silent,
    bayesian = model$bayesian
  )
}

# Create a stable hash for model parameters
create_model_hash <- function(params, debug = FALSE) {
  # Remove NULLs and exclude parameters that don't affect model specification
  params <- params[!sapply(params, is.null)]
  params$silent <- NULL  # Exclude silent - it's about output, not model specification

  # Convert complex objects to strings for stable hashing
  hash_params <- params

  # Handle family objects for stable hashing (avoid function environments)
  if (!is.null(params$family) && inherits(params$family, "family")) {
    hash_params$family <- list(
      family = params$family$family,
      link = params$family$link
    )
  }

  hash_params$formula <- paste(deparse(params$formula), collapse = "")
  hash_params$data <- digest::digest(as.data.frame(params$data), algo = "xxhash64")

  # Handle mesh objects for stable hashing
  if (!is.null(params$mesh) && inherits(params$mesh, "sdmTMBmesh")) {
    # Use only the stable components of the mesh for hashing
    hash_params$mesh <- list(
      cutoff = params$mesh$cutoff,
      max_edge = params$mesh$max_edge,
      n_vertices = nrow(params$mesh$mesh$loc),
      mesh_bounds = range(params$mesh$mesh$loc)
    )
  }

  # Create hash
  stable_state <- list(
    params = hash_params[order(names(hash_params))],
    sdmTMB_version = as.character(packageVersion("sdmTMB"))
  )

  if (debug) {
    message("Debug: Hash components:")
    str(stable_state, max.level = 2)
  }

  digest::digest(stable_state, algo = "xxhash64")
}

# Simple file-based caching
cache_model <- function(model_name, fit_dir, fit_function, check_cache = TRUE) {
  dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
  rds_file <- file.path(fit_dir, paste0(model_name, ".rds"))

  if (check_cache && file.exists(rds_file)) {
    message("Cache hit. Loading model from: ", rds_file)
    return(readRDS(rds_file))
  }

  message("Cache missing. Fitting model for: ", model_name)
  fit <- fit_function()

  saveRDS(fit, rds_file)
  message("Model saved to: ", rds_file)

  return(fit)
}

# Main caching interface
fit_cached_sdmTMB <- function(fit_dir, check_cache = TRUE, update_from = NULL,
                              model_tag = NULL, debug = FALSE, ...) {

  if (!is.null(update_from)) {
    # For model updates: merge base parameters with new ones
    base_params <- extract_model_params(update_from)
    update_args <- list(...)
    final_params <- modifyList(base_params, update_args)

    # Create hash and fit function
    model_hash <- create_model_hash(final_params, debug)
    fit_function <- local({
      function() {
        do.call(update, c(list(object = update_from), update_args))
      }
    })

  } else {
    # For new models: use parameters as provided
    final_params <- list(...)

    # Create hash and fit function
    model_hash <- create_model_hash(final_params, debug)
    fit_function <- local({
      function() {
        do.call(sdmTMB, final_params)
      }
    })
  }

  # Generate model name
  if (!is.null(model_tag)) {
    model_name <- paste0(model_tag, "-", model_hash)
  } else {
    model_name <- paste0("sdmTMB-", model_hash)
  }

  # Cache and return
  cache_model(
    model_name = model_name,
    fit_dir = fit_dir,
    fit_function = fit_function,
    check_cache = check_cache
  )
}