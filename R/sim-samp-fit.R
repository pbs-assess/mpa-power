source('R/00-utils.R')

library(sdmTMB)
library(dplyr)

# TODO PUT THIS SOMEWHERE BETTER
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds")) |>
  as_tibble()
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL"))

allocation_lu <- left_join(hbll_grid, hbll_allocations) |>
  select(-depth_m, -active_block, -area)
saveRDS(allocation_lu, file.path("data-generated", "allocation-lu.rds"))


# 1. SIMULATE ----------------------------------------------------------------
simulate_hbll <- function(fit,
                          restricted_df,
                          sim_dir = "data-generated/sim-dat",
                          check_cache = TRUE,
                          save_sim = TRUE,
                          year_covariate = 1,
                          mpa_trend = log(1.05), # 5% increase per year
                          seed = NULL,
                          formula = ~ 1 + restricted * year_covariate,
                          family = betabinomial(link = "cloglog"),
                          use_fixed_spatial_field = TRUE,
                          sigma_E = NULL, # NULL = use fitted sigma_E, 0 = no spatiotemporal variation
                          rho_V = NULL, # NULL = no AR1 deviations
                          sigma_V = NULL,
                          phi = NULL,
                          range_val = NULL,
                          tag = NULL,
                          B = NULL,
                          effect_modifer_fn = NULL,
                          mpa_trend_gmrf_sd = 0,
                          ...) {
  # Create directory for simulated data
  # dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(dir.exists(sim_dir))
  survey_type <- unique(fit$data$survey_abbrev)
  species <- unique(fit$data$species_common_name)


  # Set up simulation input parameters  ---------------------------------------
  # Get the model parameters
  b <- get_model_pars(fit)

  # Calculate sigma_V default if NULL (SD of year effects)
  if (is.null(sigma_V) && !is.null(rho_V)) {
    sigma_V <- sd(b$estimate[grepl("fyear", b$term)])
    message("Using default sigma_V = ", round(sigma_V, 3), " (SD of year effects)")
  }

  # Fixed random effects (get single draw from rf distributions)
  osp <- one_sample_posterior(fit)
  omega_s <- if (use_fixed_spatial_field) {
    osp[grepl("omega_s", names(osp))] |> matrix()
  } else {
    NULL
  }

  # Random effect SDs
  omega_s_sd <- 0
  epsilon_st_sd <- 0
  if (fit$spatial != "off") omega_s_sd <- b$estimate[b$term == "sigma_O"]
  if (fit$spatiotemporal != "off") epsilon_st_sd <- b$estimate[b$term == "sigma_E"]

  # Set sigma_E: use fitted value if NULL, otherwise use provided value
  sigma_E_sim <- if (is.null(sigma_E)) {
    epsilon_st_sd
  } else {
    sigma_E
  }
  message("Using sigma_E = ", round(sigma_E_sim, 3))

  # Prepare input data for simulation
  input_dat <- restricted_df |>
    filter(survey_abbrev %in% unique(fit$data$survey_abbrev)) |>
    sdmTMB::replicate_df(
      time_name = "year_covariate",
      time_values = year_covariate
    ) |>
    mutate(
      year = as.numeric(year_covariate),
      fyear = as.factor(year)
    )

  # # TODO: make this more general
  # if (any(grepl("log_depth", formula))) {
  #   input_dat <- input_dat |>
  #     filter(!(is.na(log_depth) | is.infinite(log_depth)))
  # }

  input_mesh <- make_mesh(input_dat, xy_cols = c("X", "Y"), mesh = fit$spde$mesh)

  # Prepare fixed random effects list
  fixed_re <- list(
    omega_s = omega_s,
    epsilon_st = NULL,
    zeta_s = NULL
  )

  # Build coefficient vector
  X <- model.matrix(object = formula, data = input_dat)
  n_coef <- ncol(X)
  coef_names <- colnames(X)

  # Calculate intercept from year effects (if conditioning model uses year as factor)
  intercept_value <- if (any(grepl("year", b$term))) {
    # mean(b[grep("year", b$term), "estimate"]) # use mean of year effects
    b[grepl("fyear", b$term), ]$estimate[nrow(b[grepl("fyear", b$term), ])] # use last year
  } else {
    b$estimate[b$term == "(Intercept)"]
  }

  if (is.null(B)) {
    B <- numeric(n_coef)
    # Coefficients - @TODO: generalise this...
    B[grep("(Intercept)", coef_names)] <- intercept_value
    # If simulating with year as factor, option to use random draws of year effects
    B[grep("fyear", coef_names)] <- sample(b[grepl("fyear", b$term), "estimate"], #
      size = length(B[grep("fyear", coef_names)]), replace = TRUE)
    B[grep("restrictedTRUE", coef_names)] <- 0
    B[grep("year_covariate$", coef_names)] <- 0 # Main effect (not interaction)
    B[grep("restricted:year_covariate", coef_names)] <- mpa_trend
    # B[grep("poly(log_depth, 2)1", coef_names)] <- b$estimate[b$term == "poly(log_depth, 2)1"]
    # B[grep("poly(log_depth, 2)2", coef_names)] <- b$estimate[b$term == "poly(log_depth, 2)2"]
  }

  # Generate offsets or weights using draws from original data
  if (!is.null(offset) && family(fit)$family != "betabinomial") {
    if (!is.null(seed)) set.seed(seed)
    offset <- sample(fit$data$offset, size = nrow(input_dat), replace = TRUE)
  }

  if (family(fit)$family == "betabinomial") {
    if (!is.null(seed)) set.seed(seed)
    weights <- sample(round(fit$data$adjusted_hook_count), size = nrow(input_dat), replace = TRUE)
  } else {
    weights <- NULL
  }

  phi <- ifelse(is.null(phi), b$estimate[b$term == "phi"], phi)
  range_val <- ifelse(is.null(range_val), b$estimate[b$term == "range"], range_val)

  message("Simulating data for ", species, " ", survey_type)
  message("- Formula: ", formula)
  message("- MPA trend: ", round(mpa_trend, 2) * 100, "%")
  message("- Coef names: ", paste(coef_names, collapse = ", "))
  message("- B: ", paste(round(B, 2), collapse = ", "))
  message("- Fitted family: ", family$family)
  message("- Parameter values: phi = ", round(phi, 1))

  # Prepare cache parameters ---------------------------------------------------
  # Extract stable parameters - avoid raw fit object instability
  final_params <- list(
    fit_params = extract_model_params(fit),  # Stable fit components
    data = input_dat,
    weights = weights,
    year_covariate = year_covariate,
    mpa_trend = mpa_trend,
    seed = seed,
    formula = formula,
    family = family,
    use_fixed_spatial_field = use_fixed_spatial_field,
    sigma_E = sigma_E_sim,
    rho_V = rho_V,
    sigma_V = sigma_V,
    time_varying = if (!is.null(rho_V)) ~ 1 else NULL,
    time_varying_type = "ar1",
    phi = phi,
    range = range_val,
    B = B,
    mpa_trend_gmrf_sd = mpa_trend_gmrf_sd,
    ...
  )

  # Create stable hash and check cache
  sim_hash <- create_model_hash(final_params)

  fname <- paste(c(species, survey_type, "sim", tag, substr(sim_hash, 1, 8)), collapse = "-") |>
    gsub("[^a-zA-Z0-9_.-]", "-", x = _)
  rds_file <- file.path(sim_dir, paste0(fname, ".rds"))

  if (check_cache && file.exists(rds_file)) {
    message("Cache hit. Loading simulation from: ", rds_file)
    return(readRDS(rds_file))
  }

  message("Cache missing. Running simulation for: ", fname)

  tvc_intercept <- if (!is.null(rho_V)) ~ 1 else NULL

  # Simulate data --------------------------------------------------------------
  sim_dat <- sdmTMB::simulate_new(
    formula = formula,
    data = input_dat,
    mesh = input_mesh,
    family = family,
    time = "year",
    time_varying = tvc_intercept,
    time_varying_type = "ar1",
    sigma_V = sigma_V,
    rho_time = rho_V,
    # rho = ar1_rho, # affects AR1 deviations of the GMRF
    sigma_E = sigma_E_sim,
    phi = phi,
    range = range_val,
    fixed_re = fixed_re,
    B = B,
    # offset = offset,
    weights = weights,
    seed = seed,
    ...
  ) |>
    as_tibble() |>
    left_join(input_dat |> distinct(X, Y, restricted))

  if (mpa_trend_gmrf_sd > 0) {
    # grab a GMRF
    svc <- sdmTMB::simulate_new(
      formula = ~ 1,
      data = input_dat,
      family = gaussian(),
      phi = 0.0001,
      range = range_val,
      sigma_O = mpa_trend_gmrf_sd,
      mesh = input_mesh,
      seed = seed * 18273,
      B = 0.01
    )
    # ggplot(svc, aes(X, Y, colour = mu)) + geom_point()
    # hist(svc$mu); abline(v = mpa_trend, col = "red"); abline(v = mpa_trend * 2, col = "blue")
    # mean(svc$mu)
    # sd(svc$mu)
    sim_dat$svc_log_lambda <- rep(svc$mu, length.out = length(sim_dat$year))
    sim_dat <- mutate(sim_dat, eta = ifelse(restricted == 1, eta + year * svc_log_lambda, eta))

    set.seed(seed * 291818)
    rbetabinom <- function(n, size, prob, phi) {
      # phi > 0; larger phi = less overdispersion (phi -> Inf approaches binomial)
      p <- rbeta(n, shape1 = prob * phi, shape2 = (1 - prob) * phi)
      rbinom(n, size = size, prob = p)
    }
    inv_cloglog <- function(x) 1 - exp(-exp(x))
    sim_dat$observed <- rbetabinom(
      length(sim_dat$observed),
      size = weights,
      prob = inv_cloglog(sim_dat$eta),
      phi = phi
    )
  }

  sim_dat$hook_count <- weights
  sim_dat$survey_abbrev <- survey_type

  # Apply modifier function to posthoc control year effects within-outside MPA network:
  # if (!is.null(effect_modifer_fn)) sim_dat <- effect_modifer_fn(sim_dat)
  # df <- data.frame(year = seq_len(25L), depletion_multiplier = seq(1, 1.3, length.out = 25L))
  # sim_dat <- left_join(sim_dat, df)
  # sim_dat$eta <- sim_dat$eta + log(sim_dat$B_recovery_val)
  # rbetabinom <- function(n, size, prob, phi) {
  #   # phi > 0; larger phi = less overdispersion (phi -> Inf approaches binomial)
  #   p <- rbeta(n, shape1 = prob * phi, shape2 = (1 - prob) * phi)
  #   rbinom(n, size = size, prob = p)
  # }
  # inv_cloglog <- function(x) 1 - exp(-exp(x))
  # sim_dat$observed <- rbetabinom(length(sim_dat$observed), size = weights, prob = inv_cloglog(sim_dat$eta), phi = phi)

  attr(sim_dat, "simulation_params") <- list(
    sim_hash = substr(sim_hash, 1, 8),
    species = species,
    survey_abbrev = survey_type,
    mpa_trend = mpa_trend,
    formula = deparse1(formula),
    time = "year",
    time_varying = deparse1(tvc_intercept),
    time_varying_type = "ar1",
    sigma_V = sigma_V,
    rho_time = rho_V,
    sigma_E = sigma_E_sim,
    phi = phi,
    range = range_val,
    B = B,
    seed = seed,
    mpa_trend_gmrf_sd = mpa_trend_gmrf_sd
  )

  # Save to cache
  if (save_sim) {
    saveRDS(sim_dat, rds_file)
    message("Simulation saved to: ", rds_file)
  }

  return(sim_dat)
}

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


# Simulate data for yelloweye --------------------------------------------------
# Make dir with tag
sim_dir <- here::here("data-generated", "test-sim-dat-svc-rates")
sim_dir <- here::here("data-generated", "test-sim-dat-no-svc-rates")
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

# Load objects for sim
ar1_parameters <- readRDS(file.path("data-generated", "ar1-parameters.rds"))
recovery_rates <- readRDS(file.path("data-generated", "recovery-rates-lambda.rds"))
restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  sf::st_drop_geometry() |>
  select(survey_abbrev, block_id, X, Y, restricted)

sp <- "yelloweye rockfish"
# Load conditioning models
fit_files <- list.files(here::here("data-generated", "fits"),
                        pattern = paste0("^", sp_to_hyphens(sp), ".*-betabinomial-on-iid-"),
                        full.names = TRUE)
reps <- 1:20
reps <- 1:5

sim_check_grid <- tidyr::crossing(rep = reps, fit_file = fit_files)

#future::plan(future::multicore, workers = 8)
future::plan(future::multisession, workers = 8)
map_fn <- furrr::future_pmap_dfr
# map_fn <- purrr::pmap_dfr # sequential
#
test_sims <- map_fn(sim_check_grid, function(rep, fit_file) {
  fit <- readRDS(fit_file)
  survey_abbrev <- unique(fit$data$survey_abbrev)[[1]]

  ar1_row <- readRDS(file.path("data-generated", "ar1-parameters.rds")) |>
    dplyr::filter(species == sp, survey_abbrev == !!survey_abbrev) |>
    dplyr::slice(1)

  lambda_val <- readRDS(file.path("data-generated", "recovery-rates-lambda.rds")) |>
    dplyr::filter(species == sp, case == "50% rate") |>
    dplyr::slice(1) |>
    dplyr::pull(lambda)

  out_file <- file.path(sim_dir, paste0(sp_to_hyphens(sp), "-", survey_abbrev, "-rep", sprintf("%03d", rep), ".rds"))

  if (!file.exists(out_file)) {
  test <- simulate_hbll(
    fit = fit,
    restricted_df = restricted_df,
    sim_dir = sim_dir,
    check_cache = FALSE,
    save_sim = FALSE,
    formula = ~ 1,
    seed = 999L + rep,
    year_covariate = 1:25,
    mpa_trend = log(lambda_val),
    use_fixed_spatial_field = TRUE,
    rho_V = ar1_row$rho_V,
    sigma_V = ar1_row$sigma_V,
    mpa_trend_gmrf_sd = 0 #0.01
  ) |>
    dplyr::mutate(replicate = rep)

  attr(test, "simulation_params")$replicate <- rep

    saveRDS(test, out_file)
  }
  # return(invisible(NULL))
  # readRDS(out_file)
}, .options = furrr::furrr_options(seed = TRUE))

future::plan(future::sequential)

catch_prop <- test_sims$observed / test_sims$hook_count
checks <- c(
  `NaN` = sum(is.nan(test_sims$observed)),
  `Inf` = sum(is.infinite(test_sims$observed)),
  all_zero = all(test_sims$observed == 0, na.rm = TRUE),
  all_NAs = all(is.na(test_sims$observed)),
  bad_range = any(catch_prop < 0 | catch_prop > 1, na.rm = TRUE),
  missing_columns = any(!colnames(test_sims) %in%
    c("survey_abbrev", "replicate",
      "year_covariate", "mpa_trend", "restricted",
      "observed", "hook_count"))
)

if (any(checks > 0)) {
  stop("Simulation check failed: ", paste(names(checks)[checks > 0], collapse = ", "))
}
message("✓ Simulation check passed")

# 2. SAMPLE --------------------------------------------------------------------
run_sampling <- function(sim_dat, plan = "status_quo") {

  filter_hbll_survey_years <- function(sim_dat) {
    sim_dat |>
      mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd")) |>
      filter(
        (survey_abbrev %in% c("HBLL INS N", "HBLL OUT N") & odd_even == "odd") |
          (survey_abbrev == "HBLL OUT S" & odd_even == "even")
      )
  }

  sample_by_plan <- function(sim_dat, sampling_effort, grouping_vars = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    group_list <- sim_dat |>
      left_join(sampling_effort, by = join_by(year, X, Y, restricted,
                                                survey_abbrev, block_id, grouping_code, survey_series_id,
                                                pfma, strata_depth, allocation)) |>
      tidyr::drop_na(n_samps) |>
      group_by(!!!syms(grouping_vars)) |>
      group_split()

    sampled_list <- purrr::map(group_list, function(g) {
      n_samps <- unique(g$n_samps)
      slice_sample(g, n = n_samps, replace = FALSE)
    })
    bind_rows(sampled_list)
  }

  # Verify sim_dat contains single replicate
  rep <- unique(sim_dat$replicate)
  if (length(rep) != 1) {
    stop("sim_dat must contain exactly one replicate, found: ", length(rep))
  }

  # Base sampling effort
  sample_effort_base <- sim_dat |>
    mutate(n_samps = allocation) |>
    select(survey_series_id, survey_abbrev, year, X, Y, block_id, grouping_code,
            pfma, strata_depth, restricted, allocation, n_samps) |>
    filter_hbll_survey_years()

  # Define sampling effort based on plan
  sampling_effort <- switch(plan,
    status_quo = sample_effort_base,

    mpas_4_years = sample_effort_base |>
      group_by(survey_abbrev) |>
      mutate(first_year = min(year)) |>
      filter(restricted == 0 | (year - first_year) %% 4 == 0) |>
      select(-first_year) |>
      ungroup(),

    historical = sim_dat |>
      filter(historical_location == 1) |>
      mutate(n_samps = allocation) |>
      select(survey_series_id, survey_abbrev, year, X, Y, block_id, grouping_code,
              pfma, strata_depth, restricted, allocation, n_samps) |>
      filter_hbll_survey_years(),

    stop("Invalid plan: '", plan, "'. Choose 'status_quo', 'mpas_4_years', or 'historical'")
  )

  plan_label <- switch(plan,
    status_quo = "status quo",
    mpas_4_years = "MPAs every 4 years",
    historical = "historical locations only"
  )

  # Sample and return (using same seed for all plans)
  sample_by_plan(
    sim_dat = sim_dat,
    sampling_effort = sampling_effort,
    grouping_vars = c("survey_abbrev", "year", "grouping_code"),
    seed = 7000 + rep
  ) |>
    mutate(plan = plan_label, replicate = rep)
}


sample_dir <- here::here("data-generated", "test-sampled-data-no-svc-rates")
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

sim_files <- list.files(sim_dir, pattern = ".*rds", full.names = TRUE)

sim_data <- purrr::map(sim_files, function(f) {
  allocation_lu <- readRDS(file.path("data-generated", "allocation-lu.rds"))
  rep <- as.integer(stringr::str_extract(basename(f), "(?<=rep)\\d+"))
  sim_dat <- readRDS(f) |> left_join(allocation_lu, by = join_by(X, Y, survey_abbrev))

  historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
    tidyr::drop_na(block_id) |>
    mutate(historical_location = 1) |>
    select(survey_abbrev, block_id, historical_location)

  sim_dat <- sim_dat |>
    mutate(historical_location = ifelse(block_id %in% historical_locations$block_id, 1, 0))
})

samp_data <- purrr::map_dfr(sim_data, function(d) {
  species <- attr(d, "simulation_params")$species |> sp_to_hyphens()
  survey_abbrev <- attr(d, "simulation_params")$survey_abbrev
  replicate <- attr(d, "simulation_params")$replicate
  plan <- "status_quo"

  out_file <- file.path(sample_dir, paste0(species, "-", survey_abbrev, "-", plan, "-rep", sprintf("%03d", replicate), ".rds"))
  if (!file.exists(out_file)) {
    samps <- run_sampling(d, plan = plan) |>
      select(-pfma, -strata_depth, -grouping_code, -allocation)

    attr(samps, "simulation_params") <- attr(d, "simulation_params")
    attr(samps, "simulation_params")$plan <- plan
    str(attr(samps, "simulation_params"))
    saveRDS(samps, out_file)
  }
  return(invisible(NULL))
})



# 3. FIT MONITORING MODEL ------------------------------------------------------

combine_hist_sim_data <- function(sim_data, hist_data, eval_year) {
  hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))
  sim_data_prep <- sim_data |>
    left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
    mutate(
      catch_count = observed,
      historical = FALSE
    )

  # Get unique surveys present in simulated data
  sim_surveys <- unique(sim_data_prep$survey_abbrev)

  # Filter historical data to only include those surveys
  hist_data_filtered <- hist_data |>
    filter(survey_abbrev %in% sim_surveys)

  combined_data <- bind_rows(hist_data_filtered, sim_data_prep) |>
    select(replicate, survey_abbrev, X, Y, block_id, restricted, historical,
           year, year_covariate, last_sampled_year,
           catch_count, hook_count, offset,
           plan) |># sim_mpa_trend, sim_ar1_scenario, sim_time_scenario) |>
    mutate(
      replicate = ifelse(historical, 0, replicate),
      catch_prop = catch_count / hook_count,
      fyear_value = ifelse(historical, year, last_sampled_year),
      fyear = as.factor(fyear_value)
    ) |>
    filter(year <= eval_year)

  if (nrow(combined_data) == 0) {
    stop("No data remaining after filtering to eval_year = ", eval_year)
  }

  return(combined_data)
}

#' Fit sdmTMB model to sampled data
fit_simulation <- function(dat,
                           formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
                           spatial = "on",
                           spatiotemporal = "iid",
                           family = betabinomial(link = "cloglog"),
                           cutoff = 10,
                           control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
                           silent = TRUE) {

  survey_type <- unique(dat$survey_abbrev)

  weights <- cdat$hook_count
  offset <- NULL

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  params <- list(
    formula = formula,
    data = dat,
    mesh = mesh,
    family = family,
    spatial = spatial,
    spatiotemporal = spatiotemporal,
    time = "year",
    weights = weights,
    offset = offset,
    silent = silent,
    control = control
  )

  fit <- local({
    tryCatch({
      do.call(sdmTMB, params)
    }, error = function(e) {
      list(error = TRUE, message = e$message)
    })
  })
}

# Monitoring model setup -------------------------------------------------------
# FIXME: move this into extract_trend estimate?
#' Extract random effect parameters from fitted model
extract_re_pars <- function(fit) {
  if (!is.null(fit$error) && fit$error) {
    return(list(
      sigma_O = NA_real_,
      sigma_E = NA_real_,
      range = NA_real_
    ))
  }

  pars <- tidy(fit, effects = "ran_pars")

  get_par <- function(term_name) {
    val <- pars$estimate[pars$term == term_name]
    if (length(val) == 0) NA_real_ else val[1]
  }

  list(
    sigma_O = get_par("sigma_O"),
    sigma_E = get_par("sigma_E"),
    range = get_par("range")
  )
}

#' Extract MPA trend estimate from fitted model
extract_trend_estimate <- function(fit, trend_param = "restricted:year_covariate") {
  if (!is.null(fit$error) && fit$error) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      converged = FALSE,
      sanity = NA_character_,
      error_msg = fit$message
    ))
  }

  coefs <- tidy(fit, conf.int = TRUE)
  trend_row <- coefs |> filter(term == trend_param)

  if (nrow(trend_row) == 0) {
    stop("Parameter ", trend_param, " not found in model")
  }

  return(list(
    estimate = trend_row$estimate,
    se = trend_row$std.error,
    ci_lower = trend_row$conf.low,
    ci_upper = trend_row$conf.high,
    converged = TRUE,
    sanity = summarise_sanity(fit),
    error_msg = NA_character_
  ))
}

species <- "yelloweye rockfish"
sample_dir <- here::here("data-generated", "test-sampled-data-no-svc-rates")
hist_dir <- here::here("data-generated", "cleaned-species-data")
fit_dir <- here::here("data-generated", "test-fits-no-svc-rates")
results_dir <- here::here("data-generated", "test-results-no-svc-rates")
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
samp_files <- list.files(file.path(sample_dir), full.names = TRUE)
hist_files <- list.files(hist_dir, pattern = paste0(sp_to_hyphens(species), ".*"), full.names = TRUE)

eval_years <- c(2030, 2034, 2038, 2042, 2046)

# FIXME: year_covariate went missing in simulated data?????

# FIXME: sampling overwrites sampling plans?



# Redefine replicates and eval_years for clarity
USE_PARALLEL <- FALSE
replicates <- 3
eval_years <- c(2030)

USE_PARALLEL <- TRUE
N_WORKERS <- 8

replicates <- 1:8
eval_years <- c(2030, 2034, 2038, 2042, 2046)

# Set up parallel backend
if (USE_PARALLEL) {
  future::plan(future::multisession, workers = N_WORKERS)
  message("Using ", N_WORKERS, " parallel workers (multisession)")
} else {
  future::plan(future::sequential)
  message("Using sequential processing")
}

#TODO: add plan tag to sampled filename: /*-plan-status_quo-.*
# Process each replicate in parallel (each rep does all eval years)
results_parallel <- furrr::future_map_dfr(
# results_parallel <- purrr::map_dfr(
  replicates,
  function(rep_num) {
    # Get files for this replicate
    pattern <- paste0("-HBLL.*rep", sprintf("%03d", rep_num), ".*")
    files <- list.files(file.path(sample_dir), full.names = TRUE, pattern = pattern)
    # Fit all eval years for this replicate
    rep_results <- purrr::map_dfr(eval_years, function(eval_yr) {
      message("Rep ", rep_num, ", eval year ", eval_yr)
      # Combine data
      cdat <- purrr::map_dfr(files, function(file) {
        sim_data <- readRDS(file)
        sim_survey <- unique(sim_data$survey_abbrev)
        hist_data <- readRDS(hist_files[grepl(sp_to_hyphens(sim_survey), hist_files)])

        comb_dat <-combine_hist_sim_data(
          sim_data = sim_data,
          hist_data = hist_data,
          eval_year = eval_yr
        )
        # FIXME: attributes for tracking
        # sim_params <- attr(sim_data, "simulation_params")
        # att_name <- paste0("simulation_params.", sim_survey)
        # attr(comb_dat, att_name) <- sim_params
        comb_dat
      }) |>
      mutate(survey_combination = "HBLL combined") |>
      mutate(year_covariate = ifelse(is.na(year_covariate), year, year_covariate))

#       filter(cdat, is.na(year_covariate)) |>
#         distinct(survey_abbrev, year, year_covariate, last_sampled_year, .keep_all = TRUE) |>
#         glimpse()
# # FIXME - why are there some empty year_covariate rows????

#       filter(cdat, !is.na(year_covariate)) |>
#       filter(year != year_covariate) |>
#       distinct(survey_abbrev, year, year_covariate) |>
#       print(n = Inf)

      # sample_plan <- attr(cdat, "simulation_params")$plan # FIXME
      sample_plan <- "status_quo" # FIXME need to be able to pull outside of loop
      fit_path <- file.path(fit_dir, paste0(species, "-", sample_plan, "-rep", sprintf("%03d", rep_num), "_", "eval", eval_yr, "_fit.rds"))
      if (!file.exists(fit_path)) {
        # Fit model
        message('Fit model')
        fit <- fit_simulation(
          dat = cdat,
          formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
          spatial = "on",
          spatiotemporal = "iid",
          # spatial = "on",
          # spatiotemporal = "iid",
          family = betabinomial(link = "cloglog"),
          cutoff = 10,
          control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
          silent = FALSE
        )
      #
        attr(fit, "simulation_params") <- attr(cdat, "simulation_params")

       saveRDS(fit, file.path(fit_dir, paste0(species, "-", sample_plan, "-rep", sprintf("%03d", rep_num), "_", "eval", eval_yr, "_fit.rds")))
      } else {
        fit <- readRDS(fit_path)
      }
     # Extract estimates
      # sink();browser();
     trend_est <- extract_trend_estimate(fit, trend_param = "restricted:year_covariate")
     re_pars <- extract_re_pars(fit)

  message('Extract estimates')
     tibble(
        replicate = rep_num,
        eval_year = eval_yr,
        estimate = trend_est$estimate,
        se = trend_est$se,
        ci_lower = trend_est$ci_lower,
        ci_upper = trend_est$ci_upper,
        converged = trend_est$converged,
        sanity = trend_est$sanity,
        error_msg = trend_est$error_msg,
        sigma_O = re_pars$sigma_O,
        sigma_E = re_pars$sigma_E,
        range = re_pars$range
      )
    })
    # Save per-replicate file
    sample_plan <- "status_quo" # FIXME want to set this once at top of run
    saveRDS(
      rep_results,
      file.path(results_dir, paste0(sp_to_hyphens(species), "-", sample_plan, "-rep", sprintf("%03d", rep_num), "_results.rds"))
    )

    return(rep_results)
  },
  .options = furrr::furrr_options(
    seed = TRUE,
    globals = c("fit_simulation", "combine_hist_sim_data", "extract_trend_estimate",
                "extract_re_pars", "sp_to_hyphens", "eval_years", "sample_dir",
                "hist_files", "species", "fit_dir", "results_dir", "summarise_sanity"),
    packages = c("dplyr", "sdmTMB", "purrr")
  )
)

meep()

# 4. POWER VISUALISATION -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs())

trend_colours <- c(
  "low"    = "#56B4E9",  # sky blue
  "moderate" = "#E69F00",  # golden orange
  "high"   = "#D55E00"   # vermillion
)

summarise_power <- function(power_df,
  by = c("species", "survey_abbrev", "sim_mpa_trend", "sim_ar1_scenario",
         "sampling_plan", "eval_year")) {
  power_df |>
    group_by(!!!syms(by)) |>
    summarise(
      # mpa_effect_label = first(mpa_effect_label),
      n_reps = n(),
      n_converged = sum(converged),
      convergence_rate = n_converged / n_reps,
      n_significant = sum(significant & converged),
      power = n_significant / n_converged,
      power_allreps = n_significant / n_reps,
      type_s_error = sum(!sign_correct & significant & converged) / n_significant,
      type_m_error = mean(ratio_to_true[significant & converged], na.rm = TRUE),
      mean_estimate = mean(estimate[converged], na.rm = TRUE),
      true_effect = first(true_effect),
      mean_bias = mean(estimate[converged] - true_effect),
      .groups = "drop"
    ) |>
    mutate(
      mpa_effect_pct = round(100 * true_effect, 2)
    )
}

get_cumul_power <- function(power_df, combo) {
  power_df0_n <- power_df0 |>
    add_count(!!!syms("eval_year"), name = "combo_n_reps")
  samples <- purrr::map_dfr(1:max(power_df$n_reps), \(x) {
    power_df0_n |>
      filter(combo_n_reps >= x) |>
      group_by(!!!syms(combo)) |>
      slice_sample(n = x, replace = FALSE) |>
      summarise_power(by = combo) |>
      mutate(n_samps = x)
  })
  # samples |> mutate(species = factor(species, levels = spp_levels))
}

species <- "yelloweye rockfish"
# fit_dir <- here::here("data-generated", "test-fits-svc-rates")
results_dir <- here::here("data-generated", "test-results-svc-rates")
results_dir <- here::here("data-generated", "test-results-no-svc-rates")
eval_years <- c(2030, 2034, 2038, 2042, 2046)

fit_files <- list.files(fit_dir)
test_fit <- readRDS(file.path(fit_dir, fit_files[2]))
# sdmTMB::sanity(test_fit)
results_files <- list.files(results_dir)
test_res <- readRDS(file.path(results_dir, results_files[2]))

all_fitted_results <- purrr::map_dfr(results_files, function(f) {
  readRDS(file.path(results_dir, f))
}) |>
  mutate(mpa_trend = log(lambda_val)) # FIXME UPSTREAM; then remove this

# FIXME: mpa_trend, sampling_plan are missing from all_fitted_results
# combo <- c("species", "survey_abbrev",
#            "sim_mpa_trend", "sim_ar1_scenario",
#            "sampling_plan", "eval_year"
# )
combo <- c("eval_year")

# Replicate-level metrics
power_df0 <- all_fitted_results |>
  mutate(
    true_effect = mpa_trend,
    significant = !(ci_lower < 0 & ci_upper > 0), # Significance: CI doesn't include 0
    sign_correct = estimate * true_effect > 0, # more robust than assuming positive effect (e.g., lingcod)
    ratio_to_true = abs(estimate) / abs(true_effect)
  )

power_df <- power_df0 |>
  summarise_power(by = combo)


# Bias check on estimate -------------------------------------------------------
b1 <- power_df0 |>
  filter(replicate <= max(replicate)) |>
ggplot(data = _) +
  geom_point(aes(x = replicate, y = estimate)) +
  geom_hline(aes(yintercept = log(lambda_val)), colour = "red") +
  facet_grid(. ~ eval_year) +
  ggtitle("Bias check")
b1

# Cumulative power plot --------------------------------------------------------
cpower <- get_cumul_power(power_df, combo)
c1 <- cpower |>
ggplot(data = ) +
  aes(x = n_samps, y = power, colour = mpa_effect_pct) +
  geom_point() +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  facet_grid(. ~ eval_year) +
  ggtitle("Cumulative power")
c1
# Main power plot --------------------------------------------------------------
p1 <- power_df |>
ggplot(data = ) +
  aes(x = eval_year, y = power, colour = mpa_effect_pct) +
  geom_point() +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  ggtitle("Power")
p1

(p1 / b1 / c1) + plot_annotation(title = "SVC rates")
