prep_hbll_data <- function(dat, bait_counts) {
  dat |>
    rename(ssid = "survey_series_id.x") |>
    left_join(bait_counts, by = c("year", "fishing_event_id", "ssid")) |>
    distinct(ssid, fishing_event_id, year, .keep_all = TRUE) |>
    mutate(count_bait_only = replace(count_bait_only, which(count_bait_only == 0), 1),
      prop_bait_hooks = count_bait_only / hook_count,
      hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
      prop_removed = 1 - prop_bait_hooks,
      offset = log(hook_count / hook_adjust_factor),
      depth_mean = mean(log(depth_m), na.rm = TRUE),
      depth_sd = sd(log(depth_m), na.rm = TRUE),
      depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1],
      depth_scaled2 = depth_scaled^2,
      fyear = as.factor(year)
    ) |>
    sdmTMB::add_utm_columns()
}

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
                        spatiotemporal = "iid",
                        use_extra_time = FALSE,
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
    use_extra_time,
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

  if (any(spatiotemporal != "off") && use_extra_time) {
    message("Spatiotemporal = ", spatiotemporal, ". Finding missing years.")
    extra_time <- sdmTMB:::find_missing_time(dat$year)
  } else {
    message("Spatiotemporal = ", spatiotemporal, ". No extra time used.")
  }

  fit <- sdmTMB::sdmTMB(
    formula = formula,
    data = dat,
    mesh = mesh,
    family = sdmTMB::nbinom2(link = "log"),
    spatial = spatial,
    spatiotemporal = spatiotemporal,
    extra_time = extra_time,
    offset = offset,
    time = time,
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
                                  crs = 4326,
                                  buffer = 40000,
                                  type = c("link", "response")) {

  if (buffer <= 0) buffer <- 1

  type <- match.arg(type)

  # Calculate pred_value based on type
  pred <- mutate(pred, pred_value = if (type == "link") est else exp(est))

  # Convert to spatial data
  pred_sf <- pred |>
    mutate(lon = 1000 * X, lat = 1000 * Y) |>
    st_as_sf(coords = c("lon", "lat"), crs = 3156)

  if (is.null(rotation)) {
    bbox <- pred_sf |>
      st_buffer(dist = buffer) |>
      st_bbox(pred_sf) |>
      st_transform(crs = crs)

    pred_sf <- st_transform(pred_sf, crs = crs)
  } else {
    pred_sf <- pred_sf |> rotate_a(a = rotation)

    bbox <- pred_sf |>
      st_buffer(dist = buffer) |>
      st_bbox(pred_sf)
  }

  if (is.null(xlim)) xlim <- bbox[c("xmin", "xmax")]
  if (is.null(ylim)) ylim <- bbox[c("ymin", "ymax")]

  # Create plot
  ggplot() +
    geom_sf(data = pacea::bc_coast |> rotate_a(a = rotation), fill = "grey90") +
    geom_sf(data = pred_sf, aes(colour = pred_value)) +
    viridis::scale_colour_viridis(option = "plasma") +
    theme_light() +
    coord_sf(xlim = xlim, ylim = ylim, crs = st_crs(pred_sf)) +
    labs(colour = if (type == "link") "est" else "exp(est)")
}
