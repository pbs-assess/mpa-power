prep_hbll_data <- function(dat, bait_counts) {
  dat |>
    rename(ssid = "survey_series_id.x") |>
    left_join(bait_counts, by = c("year", "fishing_event_id", "ssid")) |>
    distinct(ssid, fishing_event_id, year, .keep_all = TRUE) |>
    mutate(
      present = ifelse(catch_count > 0, 1, 0),
      count_bait_only = replace(count_bait_only, which(count_bait_only == 0), 1),
      prop_bait_hooks = count_bait_only / hook_count,
      hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
      prop_removed = 1 - prop_bait_hooks,
      offset = log(hook_count / hook_adjust_factor),
      log_depth = log(depth_m),
      fyear = as.factor(year)
    ) |>
    sdmTMB::add_utm_columns()
}

#' Fit sdmTMB model to HBLL survey data
# TODO: document
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

  message(
    "Predicting on ", paste(survey, collapse = ", "), " grid ",
    ifelse(length(years) == length(c(fit$fitted_time, fit$extra_time)), "with", "without"), " extra time"
  )

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

#' Simulate AR(1) deviations from trend
#'
#' Generates an AR(1) time series with specified correlation and marginal variance.
#' The implementation follows TMB's AR1_t parameterization where innovations are
#' scaled by sqrt(1-rho^2) to achieve the target marginal variance regardless of
#' temporal correlation strength.
#'
#' See: https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html
#'
#' \deqn{x_0 \sim N(0, \sigma_E^2)}
#' \deqn{x_1 = \rho x_0 + \sqrt{1-\rho^2} \varepsilon_1, \quad \varepsilon_1 \sim N(0, \sigma_E^2)}
#' \deqn{x_i = \rho x_{i-1} + \sqrt{1-\rho^2} \varepsilon_i, \quad \varepsilon_i \sim N(0, \sigma_E^2)}
#'
#' where \eqn{\rho} is the AR(1) correlation parameter and \eqn{\sigma_E} is the
#' marginal standard deviation of the process.
#'
#' @param rho AR(1) correlation parameter (-1 < \rho < 1)
#' @param sigma Marginal standard deviation of the AR(1) process (sigma_V in sdmTMB documentation - #question correct?)
#' @param years Vector of years to simulate deviations for
#'
sim_ar1_deviations <- function(rho, sigma, years) {
  n_years <- length(years)

  # Get annual innovations/epsilon/deviations:
  # innovations before applying sqrt(1-rho^2) scaling to get stationary variance; epsilon in TMB docs
  epsilon <- rnorm(n_years) * sigma

  # AR1 scaling factor:
  # innovation scaling factor (NOT an sd); sigma in TMB documentation
  # shrinks innovations such that the marginal sd = sigma_V (sigma_V in sdmTMB documentation)
  ar1_scale <- sqrt(1 - rho^2)

  # Get annual deviations:
  year_devs <- numeric(length = n_years)

  # First year from stationary distribution (uses marginal SD)
  year_devs[1] <- rnorm(1) * sigma

  # Subsequent years from AR1 process with marginal SD scaled to appropriate innovation SD
  for (i in seq(2, n_years)) {
    year_devs[i] <- rho * year_devs[i - 1] + ar1_scale * epsilon[i]
  }

  return(year_devs)
}

#' Simulate data from fitted sdmTMB model with MPA recovery trends
#'
#' @param fit Fitted sdmTMB model object
#' @param restricted_df Data frame containing spatial grid with restricted area indicators
#' @param sim_dir Directory to save simulated data (default: "data-generated/sim-dat")
#' @param check_cache Check for cached simulation (default: TRUE)
#' @param save_sim Save simulated data to cache (default: TRUE)
#' @param year_covariate Vector of time values for simulation (default: seq(0, 20, 2))
#' @param mpa_trend Log-scale trend in restricted areas (default: log(1.05) for 5% increase/year)
#' @param seed Random seed for reproducibility
#' @param formula Simulation formula (default: ~ 1 + restricted * year_covariate)
#' @param family Distribution family (default: nbinom2(link = "log"))
#' @param fixed_spatial_re Use fixed spatial random effects or use spatial sd (default: TRUE)
#' @param fixed_spatiotemporal_re Use fixed spatiotemporal random effects or use spatiotemporal sd (default: FALSE)
#' @param rho_V AR(1) correlation parameter. If NULL, no temporal AR(1) process is used (default: NULL)
#' @param sigma_V Marginal standard deviation for AR(1) process. Required if rho_V is provided.
#' @param tag Optional tag for file naming
#' @param ... Additional arguments passed to sdmTMB::sdmTMB_simulate
#'
#' @return Data frame of simulated data
#'
simulate_hbll <- function(fit,
                          restricted_df,
                          sim_dir = "data-generated/sim-dat",
                          check_cache = TRUE,
                          save_sim = TRUE,
                          year_covariate = 1,
                          mpa_trend = log(1.05), # 5% increase per year
                          seed = NULL,
                          formula = ~ 1 + restricted * year_covariate,
                          family = nbinom2(link = "log"),
                          fixed_spatial_re = TRUE,
                          fixed_spatiotemporal_re = FALSE,
                          rho_V = NULL, # NULL = no AR1 deviations
                          sigma_V = NULL,
                          phi = NULL,
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
    rho_V, # NULL if not using AR1 deviations on mean
    sigma_V, # NULL if not provided
    fit$spde$mesh,  # Mesh from fitted model
    phi,
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
  omega_s_sd <- 0
  epsilon_st_sd <- 0
  if (fit$spatial != "off") omega_s_sd <- b$estimate[b$term == "sigma_O"]
  if (fit$spatiotemporal != "off") epsilon_st_sd <- b$estimate[b$term == "sigma_E"]

  # Prepare input data for simulation
  input_dat <- restricted_df |>
    filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
    # select(X, Y, restricted) |>
    sdmTMB::replicate_df(
      time_name = "year_covariate",
      time_values = year_covariate
    ) |>
    mutate(
      restricted = restricted,
      year = as.numeric(year_covariate)
    )

# TODO: make this more general
  if (any(grepl("log_depth", formula))) {
    input_dat <- input_dat |>
      filter(!(is.na(log_depth) | is.infinite(log_depth)))
  }

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
  # Build coefficient vector
  X <- model.matrix(object = formula, data = input_dat)
  n_coef <- ncol(X)
  coef_names <- colnames(X)
  B <- numeric(n_coef)
  # Coefficients
  B[grep("(Intercept)", coef_names)] <- intercept_value # Only if intercept in formula
  B[grep("restrictedTRUE", coef_names)] <- 0
  B[grep("year_covariate$", coef_names)] <- 0  # Main effect (not interaction)
  B[grep("restricted:year_covariate", coef_names)] <- mpa_trend

  if (!is.null(rho_V)) {
    if (is.null(sigma_V)) {
      stop("sigma_V must be provided if rho_V is provided")
    } else {
      message("Simulating with AR1 process: rho = ", round(rho_V, 2), ", sigma_V = ", round(sigma_V, 2))
    }

    # generate AR1 deviations
    message("- Generating AR1 deviations for ", length(year_covariate))
    year_devs <- sim_ar1_deviations(rho = rho_V,
      sigma = sigma_V,
      years = year_covariate
    )
    year_indices <- grep("as.factor\\(year\\)", coef_names)

    if (length(year_indices) > 0) {
      message("- Using ", length(year_indices), " year factors")
        if ("(Intercept)" %in% coef_names) {
          stop("Don't use both intercept and factor years. Use either '~ 0 + as.factor(year) + ...' or '~ 1 + restricted * year_covariate'")
        }
      # Year coefficients (for each factor level)
      last_year_intercept <- b |>
        filter(grepl("^fyear", term)) |>
        mutate(year = as.numeric(gsub("fyear", "", term))) |>
        filter(year == max(year))
      B[year_indices] <- year_devs + last_year_intercept$estimate # add AR1 deviations starting from last year intercept
    # test <- b |> filter(grepl("fyear", term)) |> pull(estimate)
    # plot(c(test, B[year_indices]))
    }
  }

  log_mean_offset <- log(mean(exp(fit$data$offset), na.rm = TRUE))
  offset <- rep(log_mean_offset, nrow(input_dat))

  message("- MPA trend: ", round(mpa_trend, 2) * 100, "%")

  # Simulate data
  sim_dat <- sdmTMB::sdmTMB_simulate(
    formula = formula,
    data = input_dat,
    mesh = input_mesh,
    family = family,
    time = "year",
    # rho = ar1_rho, # affects AR1 deviations of the GMRF
    sigma_E = if (fixed_spatiotemporal_re) epsilon_st_sd else 0,
    phi = if (is.null(phi)) b$estimate[b$term == "phi"] else phi,
    range = b$estimate[b$term == "range"],
    fixed_re = fixed_re,
    B = B,
    offset = offset,
    seed = seed,
    ...
  ) |>
    as_tibble()

  sim_dat$offset <- log_mean_offset

  # # Add simulation parameters as attributes for tracking
  # attr(sim_dat, "simulation_params") <- list(
  #   survey_abbrev = unique(fit$data$survey_abbrev),
  #   B = B,
  #   B_names = coef_names,
  #   phi = if (is.null(phi)) b$estimate[b$term == "phi"] else phi,
  #   range = b$estimate[b$term == "range"],
  #   mpa_trend = mpa_trend,
  #   rho_V = rho_V,
  #   sigma_V = sigma_V,
  #   seed = seed
  # )

  # Save to cache if save_sim is TRUE
  if (save_sim) {
    message("Cache file will be saved to: ", rds_file)
    saveRDS(sim_dat, file = rds_file)
    writeLines(current_hash, con = hash_file)
  }

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
sample_by_plan <- function(sim_dat,
  sampling_effort,
  grouping_vars = NULL
  ) {
  group_list <- sim_dat |>
    left_join(sampling_effort) |>
    group_by(!!!syms(grouping_vars)) |>
    group_split()

  sampled_list <- purrr::map(group_list, function(g) {
    n_samps <- unique(g$n_samps)
    slice_sample(g, n = n_samps, replace = FALSE)
  })
  bind_rows(sampled_list)
}

#' Simple function to plot sampling plans
#'
#' @param sampled_data sampled `sim_dat`
#' @param plan_name plot tite
#'
#' @return ggplot object
#'
plot_sampling_plan <- function(sampled_data, plan_name) {

  # Create summary for text labels
  samp_summary <- sampled_data |>
    group_by(year, restricted) |>
    summarise(n = n(), .groups = "drop") |>
    filter(restricted == 1)

  # Prepare plot data
  plot_dat <- sampled_data |>
    XY_to_sf()

  # Create the plot
  ggplot(data = plot_dat) +
    geom_sf(data = mpa_shape_simplified, fill = "#0072B2",
      colour = NA, alpha = 0.3) +
    # scale_fill_manual(name = "MPA status",
    #   values = c("not identified as concern" = "#0072B2",
    #     "currently restricted" = "#D55E00", "identified as concern" = "#F0E442",
    #     "not applicable" = "#999999"), na.value = "grey90") +
    # ggnewscale::new_scale_fill() +
    geom_sf(data = pacea::bc_coast, fill = "grey94", colour = "grey90") +
    ggnewscale::new_scale_fill() +
    geom_sf(aes(colour = eta, shape = factor(restricted)), size = 1.2) +
    scale_shape_manual(name = "Restricted", values = c(`0` = 21, `1` = 19)) +
    scale_colour_viridis_c(name = "eta", option = "A", end = 0.8) +
    facet_wrap(~ year, nrow = 5) +
    plot_limits_combined +
    theme(legend.position = "bottom",
      # legend.position.inside = c(0.87, 0.1),
      legend.box = "horizontal") +
    geom_sf_text(data = samp_summary |> mutate(X = -126, Y = 54) |>
      XY_to_sf(crs_from = 4326, crs_to = 4326),
      aes(x = X, y = Y, label = paste0("n = ", n)), size = 3) +
    ggtitle(plan_name)
}
