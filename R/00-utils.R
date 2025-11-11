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

  # Use deparse1() if available (R >= 4.0) for stable formula representation
  # deparse1() always returns a single string, unlike deparse() which can return a vector
  # See: https://stackoverflow.com/questions/70850546/why-does-deparse-return-a-vector-of-length-two-here
  # Falls back to deparse() for older R versions
  if (exists("deparse1")) {
    hash_params$formula <- deparse1(params$formula)
  } else {
    hash_params$formula <- paste(deparse(params$formula), collapse = "")
  }

  # Make data hashing robust across systems
  # Round numeric columns to avoid floating point precision differences
  data_for_hash <- as.data.frame(params$data)
  rownames(data_for_hash) <- NULL

  # Sort columns alphabetically for consistent ordering across systems
  data_for_hash <- data_for_hash[, order(names(data_for_hash)), drop = FALSE]

  # Exclude spatial coordinate columns - spatial info is captured in mesh hash
  # This prevents GDAL/PROJ version differences from affecting the hash
  coord_cols <- c("X", "Y", "lon", "lat", "longitude", "latitude", "x", "y")
  data_for_hash <- data_for_hash[, !names(data_for_hash) %in% coord_cols, drop = FALSE]

  # Round numeric columns to 10 decimal places for stable hashing
  numeric_cols <- sapply(data_for_hash, is.numeric)
  data_for_hash[numeric_cols] <- lapply(data_for_hash[numeric_cols], round, digits = 10)

  # Sort rows to ensure hash is independent of row order
  # This makes mesh generation order-independent for caching
  data_for_hash <- data_for_hash[do.call(order, data_for_hash), , drop = FALSE]

  hash_params$data <- digest::digest(data_for_hash, algo = "xxhash64")

  # Hash numeric vectors (offset, weights) to avoid R version differences in as.character()
  # Different R versions format floats differently, causing hash mismatches
  if (!is.null(params$offset)) {
    hash_params$offset <- digest::digest(round(params$offset, 10), algo = "xxhash64")
  }
  if (!is.null(params$weights)) {
    hash_params$weights <- digest::digest(round(params$weights, 10), algo = "xxhash64")
  }

  # Handle mesh objects for stable hashing
  if (!is.null(params$mesh) && inherits(params$mesh, "sdmTMBmesh")) {
    # For cross-platform stability, only hash coarse spatial extent
    # Exact mesh triangulation (n_vertices) depends on precise X,Y coordinates
    # which vary with GDAL/PROJ versions. What matters is the general spatial domain.
    # Round bounds to nearest 1km for robustness to coordinate transformation differences
    mesh_bounds <- range(params$mesh$loc_xy)
    hash_params$mesh <- list(
      data_bounds = round(mesh_bounds / 1000) * 1000,  # Round to nearest km
      xy_cols = params$mesh$xy_cols,
      manifold = params$mesh$mesh$manifold
    )
  }

  # Create hash using string-based approach for platform independence
  hash_params_ordered <- hash_params[order(names(hash_params))]

  # Convert all components to strings and concatenate
  hash_components <- sapply(names(hash_params_ordered), function(nm) {
    paste0(nm, "=", paste(as.character(hash_params_ordered[[nm]]), collapse = ","))
  })

  hash_string <- paste(hash_components, collapse = "|")

  if (debug) {
    message("\n=== Hash Debug Information ===")
    message("Individual hash components:")
    for (nm in names(hash_params_ordered)) {
      component_str <- paste0(nm, "=", paste(as.character(hash_params_ordered[[nm]]), collapse = ","))
      message(sprintf("  %s: %s", nm, substr(component_str, 1, 100)))  # Truncate long strings
      if (nchar(component_str) > 100) message(sprintf("    ... (truncated, length: %d)", nchar(component_str)))
    }
    message("\nFull hash string:")
    message(substr(hash_string, 1, 200))
    if (nchar(hash_string) > 200) message(sprintf("... (truncated, full length: %d)", nchar(hash_string)))
    message("\nFinal hash: ", digest::digest(hash_string, algo = "xxhash64"))
    message("===============================\n")
  }

  digest::digest(hash_string, algo = "xxhash64")
}

# # Simple file-based caching
# cache_model <- function(model_name, fit_dir, fit_function, check_cache = TRUE) {
#   dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
#   rds_file <- file.path(fit_dir, paste0(model_name, ".rds"))

#   if (check_cache && file.exists(rds_file)) {
#     message("Cache hit. Loading model from: ", rds_file)
#     return(readRDS(rds_file))
#   }

#   message("Cache missing. Fitting model for: ", model_name)
#   fit <- fit_function()
# browser()
#   saveRDS(fit, rds_file)
#   message("Model saved to: ", rds_file)

#   return(fit)
# }

# TODO document
summarise_sanity <- function(fit) {
    sanity_list <- sanity(fit, silent = TRUE) |> unlist()
    sanity_false_names <- names(sanity_list)[!sanity_list]
    sanity_str <- ifelse(length(sanity_false_names) == 0, "ok", paste(sanity_false_names, collapse = "; "))
    gsub("_ok", "", sanity_str)
  }

update_collapsed_rf <- function(fit) {
  rp <- sdmTMB::tidy(fit, "ran_pars")
  rp$collapsed <- rp$estimate < 0.01

  spatial <- if ("sigma_O" %in% rp$term && rp$collapsed[rp$term == "sigma_O"]) {
    "off"
  } else {
    fit$spatial
  }

  spatiotemporal <- if ("sigma_E" %in% rp$term && rp$collapsed[rp$term == "sigma_E"]) {
    "off"
  } else {
    fit$spatiotemporal
  }

  needs_refit <- (spatial != fit$spatial) || (spatiotemporal != fit$spatiotemporal)

  list(
    spatial = spatial,
    spatiotemporal = spatiotemporal,
    needs_refit = needs_refit
  )
}

#' Setup parallel or sequential processing
#'
#' @param use_parallel Logical. Use parallel processing?
#' @param n_workers Number of workers for parallel processing
#'
#' @return A map function (either purrr::map or furrr::future_map)
setup_parallel <- function(use_parallel, n_workers = NULL) {
  if (use_parallel) {
    if (is.null(n_workers)) n_workers <- floor(future::availableCores() / 2)

    # Use multicore on hake server (Unix fork), multisession elsewhere (Windows-safe)
    if (Sys.info()['user'] %in% c("dunic", "anderson")) {
      future::plan(future::multicore, workers = n_workers)
      message("Using ", n_workers, " parallel workers (multicore)")
    } else {
      future::plan(future::multisession, workers = n_workers)
      message("Using ", n_workers, " parallel workers (multisession)")
    }
    return(furrr::future_map)
  } else {
    future::plan(future::sequential)
    message("Using sequential processing")
    return(purrr::map)
  }
}

clean_family_name <- function(fit) {
    if (!is.null(fit$family$delta)) {
      # Delta model - sdmTMB style: "Delta family1(link1), family2(link2)"
      fam1 <- fit$family$family[[1]]
      fam2 <- fit$family$family[[2]]
      link1 <- fit$family$link[[1]]
      link2 <- fit$family$link[[2]]
      paste0("Delta ", fam1, "(", link1, "), ", fam2, "(", link2, ")")
    } else {
      # Regular model - sdmTMB style: "family(link)"
      fam <- fit$family$family
      link <- fit$family$link
      paste0(fam, "(", link, ")")
    }
  }