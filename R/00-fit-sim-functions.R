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

#' Draw a single sample from the posterior (but not posterior??) distribution of random effects
#' #question - what do I call the not posterior?
#'
#' @param object A fitted sdmTMB model object
#' @param use_names Logical. Return named vector. (default: TRUE)
#'
#' @return A numeric vector containing parameter values with random effects sampled from
#'   their posterior distribution and fixed effects at their MLE values.
#'
#' @details
#' TMB `MC()` to draw a single sample from the posterior distribution of random
#' effects while keeping fixed effects at their maximum likelihood estimates.
#'
#' See: https://github.com/pbs-assess/sdmTMB/blob/228363611f891462b6cb9b50fd19afb5eab5d5e0/R/residuals.R#L491-L501
#'
one_sample_posterior <- function(object, use_names = TRUE) {
  tmp <- object$tmb_obj$env$MC(n = 1L, keep = TRUE, antithetic = FALSE)
  re_samp <- as.vector(attr(tmp, "samples"))
  lp <- object$tmb_obj$env$last.par.best
  p <- numeric(length(lp))
  fe <- object$tmb_obj$env$lfixed()
  re <- object$tmb_obj$env$lrandom()
  p[re] <- re_samp
  p[fe] <- lp[fe]
  if (use_names) names(p) <- names(lp)
  p
}

#' Simulate data from fitted sdmTMB model with MPA recovery trends
#'
#' @param fit Fitted sdmTMB model object
#' @param restricted_df Data frame containing spatial grid with restricted area indicators
#' @param sim_dir Directory to save simulated data (default: "data-generated/sim-dat")
#' @param check_cache Check for cached simulation (default: TRUE)
#' @param year_covariate Vector of time values for simulation (default: seq(0, 20, 2))
#' @param mpa_trend Log-scale trend in restricted areas (default: log(1.05) for 5% increase/year)
#' @param seed Random seed for reproducibility
#' @param formula Simulation formula (default: ~ 1 + restricted * year_covariate)
#' @param family Distribution family (default: nbinom2(link = "log"))
#' @param fixed_spatial_re Use fixed spatial random effects or use spatial sd (default: TRUE)
#' @param fixed_spatiotemporal_re Use fixed spatiotemporal random effects or use spatiotemporal sd (default: FALSE)
#' @param tag Optional tag for file naming
#' @param ... Additional arguments passed to sdmTMB::sdmTMB_simulate
#'
#' @return Data frame of simulated data
#'
simulate_hbll <- function(fit,
                          restricted_df,
                          sim_dir = "data-generated/sim-dat",
                          check_cache = TRUE,
                          year_covariate = seq(from = 0, to = 20, by = 2),
                          mpa_trend = log(1.05), # 5% increase per year
                          seed = NULL,
                          formula = ~ 1 + restricted * year_covariate,
                          family = nbinom2(link = "log"),
                          fixed_spatial_re = TRUE,
                          fixed_spatiotemporal_re = FALSE,
                          tag = NULL,
                          ...) {

  # Create directory for simulated data
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

  # Generate filename based on fit and parameters
  survey_type <- unique(fit$data$survey_abbrev)
  species <- unique(fit$data$species_common_name)

  fname <- paste(c(species, survey_type, "sim", tag), collapse = "-") |>
    gsub("[^a-zA-Z0-9_.-]", "-", x = _)
  rds_file <- file.path(sim_dir, paste0(fname, ".rds"))
  hash_file <- file.path(sim_dir, paste0(fname, ".hash"))

  # Create sim state for hashing (similar to fit_hbll)
  sim_state <- list(
    fit$data,  # Original data used for fitting
    restricted_df,
    year_covariate,
    mpa_trend,
    seed,
    formula,
    family,
    fixed_spatial_re,
    fixed_spatiotemporal_re,
    fit$spde$mesh,  # Mesh from fitted model
    packageVersion("sdmTMB"),
    list(...)
  )
  current_hash <- digest::digest(sim_state)

  # Check cache
  if (check_cache && file.exists(rds_file) && file.exists(hash_file)) {
    cached_hash <- readLines(hash_file, warn = FALSE)
    if (identical(cached_hash, current_hash)) {
      message("Cache hit. Loading simulation from: ", rds_file)
      return(readRDS(rds_file))
    }
  }

  message(
    "Cache missing or invalid. Running simulation for: ", fname,
    "\nCache file will be saved to: ", rds_file
  )

  # Get the model parameters
  b <- get_model_pars(fit)

  # Fixed random effects (get single draw from rf distributions)
  osp <- one_sample_posterior(fit)
  omega_s <- if (fixed_spatial_re) {
    osp[grepl("omega_s", names(osp))] |> matrix()
  } else {
    NULL
  }

  # Random effect SDs
  omega_s_sd <- b$estimate[b$term == "sigma_O"]
  epsilon_st_sd <- b$estimate[b$term == "sigma_E"]
  if (fit$spatiotemporal == "off") epsilon_st_sd <- 0

  # Prepare input data for simulation
  input_dat <- restricted_df |>
    filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
    select(X, Y, restricted) |>
    sdmTMB::replicate_df(
      time_name = "year_covariate",
      time_values = year_covariate
    ) |>
    mutate(
      restricted = restricted,
      year = as.numeric(year_covariate)
    )

  input_mesh <- make_mesh(input_dat, xy_cols = c("X", "Y"), mesh = fit$spde$mesh)

  # Prepare fixed random effects list
  fixed_re <- list(
    omega_s = omega_s,
    epsilon_st = NULL,
    zeta_s = NULL
  )

  # Calculate intercept from year effects (if using year as factor)
  intercept_value <- if (any(grepl("year", b$term))) {
    mean(b[grep("year", b$term), "estimate"])
  } else {
    b$estimate[b$term == "(Intercept)"]
  }

  # Coefficient vector: (Intercept), restrictedTRUE, year_covariate, restrictedTRUE:year_covariate
  B <- c(intercept_value, 0, 0, mpa_trend)

  # Simulate data
  sim_dat <- sdmTMB::sdmTMB_simulate(
    formula = formula,
    data = input_dat,
    mesh = input_mesh,
    family = family,
    time = "year",
    # rho = rho,
    sigma_E = if (fixed_spatiotemporal_re) epsilon_st_sd else 0,
    phi = b$estimate[b$term == "phi"],
    range = b$estimate[b$term == "range"],
    fixed_re = fixed_re,
    B = B,
    seed = seed,
    ...
  ) |>
    as_tibble()

  # Save to cache
  saveRDS(sim_dat, file = rds_file)
  writeLines(current_hash, con = hash_file)

  return(sim_dat)
}

# Example usage:
# sim_dat <- simulate_mpa_data(
#   fit = fit_IN,
#   restricted_df = restricted_df,
#   year_covariate = seq(0, 20, 2),
#   mpa_trend = log(1.05),  # 5% increase per year
#   seed = 714
# )


#' Sample simulated data according to a sampling plan
#'
#' This function samples simulated data based on a specified sampling effort plan.
#' It groups the data by specified variables and samples the appropriate number
#' of observations from each group according to the sampling effort specifications.
#'
#' @param sim_dat A data frame containing the simulated data to be sampled.
#'   Must contain columns that match the grouping variables and join keys
#'   for the sampling effort data frame.
#' @param sampling_effort A data frame containing sampling effort specifications.
#'   Must have a column `n_samps` indicating the number of samples to take
#'   from each group, and columns that can be joined with `sim_dat` (typically
#'   survey identifiers like `survey_abbrev`).
#' @param grouping_vars A character vector of column names to group by for sampling.
#'   Common grouping variables include `c("survey_abbrev", "year")` to sample
#'   separately for each survey and year combination. If `NULL`, no grouping
#'   is applied and sampling is done across the entire dataset.
#'
#' @return A data frame containing the sampled observations. The structure
#'   matches the input `sim_dat` but with fewer rows based on the sampling
#'   effort specifications.
#'
#' @details
#' The function works by:
#' 1. Joining the simulated data with the sampling effort specifications
#' 2. Grouping the data by the specified grouping variables
#' 3. For each group, sampling the number of observations specified in `n_samps`
#' 4. Combining all sampled groups back into a single data frame
#'
#' This is particularly useful for simulating survey sampling scenarios where
#' different areas or time periods may have different sampling intensities.
#'
sample_by_plan <- function(sim_dat, sampling_effort, grouping_vars = NULL) {
  group_list <- sim_dat |>
    left_join(sampling_effort) |>
    group_by(!!!syms(grouping_vars)) |>
    group_split()

  sampled_list <- map(group_list, function(g) {
    n_samps <- unique(g$n_samps)
    slice_sample(g, n = n_samps, replace = FALSE)
  })
  bind_rows(sampled_list)
}