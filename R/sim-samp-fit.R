source('R/00-utils.R')

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs())

# TODO PUT THIS SOMEWHERE BETTER
run_ssf <- function(species, mpa_rate, tag) {
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
    message("- MPA trend: ", round(mpa_trend, 4) * 100, "%")
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

    # message("Cache missing. Running simulation for: ", fname)

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

  ################################################################################
  # FULL RUN SETUP CONFIG --------------------------------------------------------
  lambda_values <- exp(log(c(1.05, 1.10, 1.25, 1.50)) / 25)
  lambda_lu <- data.frame(
    percent_increase = c("5%", "10%", "25%", "50%"),
    lambda = lambda_values
  )

  sp <- species
  mpa_rate <- mpa_rate
  # sp <- "yelloweye rockfish"
  # mpa_rate <- "10%"
  lambda_val <- lambda_lu |>
    filter(percent_increase == mpa_rate) |>
    pull(lambda)

  # Type of simulation model
  # tag <- "no-svc-rates"
  # mpa_trend_gmrf_sd <- 0

  # tag <- "svc-rates-sd-0.01"
  # mpa_trend_gmrf_sd <- 0.01

  # FIXME: easy to make mistake - this should become a check
  if (tag == "svc-rates-sd-0.01") mpa_trend_gmrf_sd <- 0.01
  if (tag == "no-svc-rates") mpa_trend_gmrf_sd <- 0

  sampling_plan <- "status_quo" # eventually include as input
  eval_years <- c(2030, 2034, 2038, 2042, 2046) # eventually include as input
  reps <- 1:25 # eventually include as input

message("Running ssf for ", sp, " ", mpa_rate, " ", tag)
message("  - Recovery rate: ", round(lambda_val, 4))
message("  - Repetitions: ", length(reps))
message("  - mpa_trend_gmrf_sd: ", mpa_trend_gmrf_sd)
message("  - sampling_plan: ", sampling_plan)
message("  - eval_years: ", paste(eval_years, collapse = ", "))


  hist_dir <- here::here("data-generated", "cleaned-species-data")

  # 1. SIMULATE
  sim_dir <- here::here("data-generated", paste0("1-sim-dat-", tag))
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

  # 2. SAMPLE
  sample_dir <- here::here("data-generated", paste0("2-sampled-data-", tag))
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

  # 3. FIT
  fit_dir <- here::here("data-generated", paste0("3-fits-", tag))
  dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)

  # 4. RESULTS
  results_dir <- here::here("data-generated", paste0("4-results-", tag))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  ################################################################################


  # Simulate data for yelloweye --------------------------------------------------
  # Load objects for sim
  ar1_parameters <- readRDS(file.path("data-generated", "ar1-parameters.rds"))
  restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
    sf::st_drop_geometry() |>
    select(survey_abbrev, block_id, X, Y, restricted)

  message("Simulating data for ", sp)
  message("  - Recovery rate: ", round(lambda_val, 4))
  message("  - Repetitions: ", length(reps))

  # Load conditioning models
  fit_files <- list.files(here::here("data-generated", "fits"),
                          pattern = paste0("^", sp_to_hyphens(sp), ".*-betabinomial-on-iid-"),
                          full.names = TRUE)

  sim_task_grid <- tidyr::crossing(rep = reps, fit_file = fit_files, lambda_val = lambda_val)

  future::plan(future::multisession, workers = 8)
  map_fn <- furrr::future_pmap
  # future::plan(future::sequential)
  # map_fn <- purrr::pmap # sequential
  #
  reps <- 1:20
  sims <- map_fn(sim_task_grid, function(rep, fit_file, lambda_val) {
    fit <- readRDS(fit_file)
    survey_abbrev <- unique(fit$data$survey_abbrev)[[1]]
    message("Simulating data for ", sp, "\n  - Survey: ", survey_abbrev,
      "\n  - Rep: ", rep, "\n  - Lambda: ", round(lambda_val, 4))

    ar1_row <- readRDS(file.path("data-generated", "ar1-parameters.rds")) |>
      dplyr::filter(species == sp, survey_abbrev == !!survey_abbrev)

    out_file <- file.path(sim_dir, paste0(sp_to_hyphens(sp), "-", survey_abbrev, "-", round(lambda_val, 4), "-rep", sprintf("%03d", rep), ".rds"))

    if (!file.exists(out_file)) {
      message("Cache missing. Running simulation for: ", basename(out_file))
        year_covariate <- 1:25
      sim_dat <- simulate_hbll(
        fit = fit,
        restricted_df = restricted_df,
        sim_dir = sim_dir,
        check_cache = FALSE,
        save_sim = FALSE,
        formula = ~ 1,
        seed = 999L + rep,
        year_covariate = year_covariate,
        mpa_trend = log(lambda_val),
        use_fixed_spatial_field = TRUE,
        rho_V = ar1_row$rho_V,
        sigma_V = ar1_row$sigma_V,
        mpa_trend_gmrf_sd = mpa_trend_gmrf_sd
      ) |>
        dplyr::mutate(replicate = rep)

      attr(sim_dat, "simulation_params")$replicate <- rep

      saveRDS(sim_dat, out_file)
    } else {
      message("Cache hit. Loading simulation from: ", basename(out_file))
      # return(invisible(NULL))
      readRDS(out_file)
    }
  }, .options = furrr::furrr_options(seed = TRUE))
  # })
  meep()

  future::plan(future::sequential)

  # SIMULATION CHECKS --------
  test_sim <- purrr::map_dfr(c("HBLL INS N", "HBLL OUT N", "HBLL OUT S"), function(survey) {
    readRDS(file.path(
      sim_dir, paste0(sp_to_hyphens(sp),
        "-", survey,
        "-", round(lambda_val, 4),
        "-rep", sprintf("%03d", 1), ".rds")
      )
    )
  })

  ggplot(test_sim) +
    aes(X, Y, colour = observed, shape = factor(restricted)) +
    geom_point() +
    scale_colour_viridis_c(trans = ggsidekick::fourth_root_power_trans()) +
    facet_wrap(~ year)

  catch_prop <- test_sim$observed / test_sim$hook_count
  checks <- c(
    `NaN` = sum(is.nan(test_sim$observed)),
    `Inf` = sum(is.infinite(test_sim$observed)),
    all_zero = all(test_sim$observed == 0, na.rm = TRUE),
    all_NAs = all(is.na(test_sim$observed)),
    bad_range = any(catch_prop < 0 | catch_prop > 1, na.rm = TRUE),
    missing_columns = any(!colnames(test_sim) %in%
      c("survey_abbrev", "replicate",
        "year_covariate", "mpa_trend", "restricted",
        "observed", "hook_count"))
  )

  if (any(checks > 0)) {
    warning("Simulation check failed: ", paste(names(checks)[checks > 0], collapse = ", "))
  } else {
    message("✓ Simulation check passed")
  }

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

  # Sample data ----
  message("Using sim data from: ", basename(sim_dir))
  message("Out dir: ", basename(sample_dir))
  message("Sampling data for ", sp)
  message("  - Plan: ", sampling_plan)
  message("  - Lambda: ", round(lambda_val, 4))
  message("  - Repetitions: ", length(reps))

  # TODO: DON'T LOAD EVERYTHING INTO MEMORY
  sim_files <- list.files(sim_dir,
    pattern = paste0(sp_to_hyphens(sp),
    "-(HBLL INS N|HBLL OUT N|HBLL OUT S)-",
    round(lambda_val, 4), "-rep\\d{3}.*rds"),
    full.names = TRUE)

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

  samp_data <- purrr::walk(sim_data, function(d) {
    species <- attr(d, "simulation_params")$species |> sp_to_hyphens()
    survey_abbrev <- attr(d, "simulation_params")$survey_abbrev
    replicate <- attr(d, "simulation_params")$replicate
    sampling_plan <- "status_quo" # FIXME - hard coded sampling plan

    out_file <- file.path(sample_dir, paste0(species,
      "-", survey_abbrev,
      "-", round(lambda_val, 4),
      "-", sampling_plan,
      "-rep", sprintf("%03d", replicate),
      ".rds")
    )

    if (!file.exists(out_file)) {
      samps <- run_sampling(d, plan = sampling_plan) |>
        select(-pfma, -strata_depth, -grouping_code, -allocation)
      attr(samps, "simulation_params") <- attr(d, "simulation_params")
      attr(samps, "simulation_params")$sampling_plan <- sampling_plan
      str(attr(samps, "simulation_params"))
      saveRDS(samps, out_file)
    } else {
      samps <- readRDS(out_file)
    }
  })
  # attr(samp_data[[2]], "simulation_params")

  # SAMPLING CHECKS --------
  samp_data <- purrr::map_dfr(c("HBLL INS N", "HBLL OUT N", "HBLL OUT S"), function(survey) {
    readRDS(file.path(sample_dir, paste0(sp_to_hyphens(sp),
      "-", survey,
      "-", round(lambda_val, 4),
      "-", sampling_plan,
      "-rep", sprintf("%03d", 1), ".rds")))
  })

  ggplot(samp_data) +
    aes(X, Y, colour = observed, shape = factor(restricted)) +
    geom_point() +
    scale_colour_viridis_c(trans = ggsidekick::fourth_root_power_trans()) +
    scale_shape_manual(values = c(21, 19)) +
    facet_wrap(~ year)


  # 3. FIT MONITORING MODEL ------------------------------------------------------

  # combine_hist_sim_data <- function(sim_data, hist_data, eval_year) {
  #   hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))
  #   sim_data_prep <- sim_data |>
  #     left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
  #     mutate(
  #       catch_count = observed,
  #       historical = FALSE
  #     )

  #   # Get unique surveys present in simulated data
  #   sim_surveys <- unique(sim_data_prep$survey_abbrev)

  #   # Filter historical data to only include those surveys
  #   hist_data_filtered <- hist_data |>
  #     filter(survey_abbrev %in% sim_surveys)

  #   combined_data <- bind_rows(hist_data_filtered, sim_data_prep) |>
  #     select(replicate, survey_abbrev, X, Y, block_id, restricted, historical,
  #            year, year_covariate, last_sampled_year,
  #            catch_count, hook_count, offset,
  #            plan) |># sim_mpa_trend, sim_ar1_scenario, sim_time_scenario) |>
  #     mutate(
  #       replicate = ifelse(historical, 0, replicate),
  #       catch_prop = catch_count / hook_count,
  #       fyear_value = ifelse(historical, year, last_sampled_year),
  #       fyear = as.factor(fyear_value)
  #     ) |>
  #     filter(year <= eval_year)

  #   if (nrow(combined_data) == 0) {
  #     stop("No data remaining after filtering to eval_year = ", eval_year)
  #   }

  #   return(combined_data)
  # }


  # ---
  # species <- "yelloweye rockfish"
  # plan <- "status_quo"
  # sample_dir <- here::here("data-generated", "test-sampled-data-no-svc-rates")
  # rep_num <- 1
  # samp_files <- list.files(sample_dir, pattern = paste0(sp_to_hyphens(species), ".*", plan, ".*", "rep", sprintf("%03d", rep_num), ".*rds"), full.names = TRUE)
  # test <- purrr::map_dfr(samp_files, function(f) {
  #   hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))
  #   readRDS(f) |>
  #     mutate(
  #       catch_count = observed,
  #       historical = FALSE
  #     ) |>
  #     left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
  #     mutate(
  #       year_post_imp = year, # post implementation year
  #       year = last_sampled_year + year, # this is calendar year
  #       fyear = as.factor(last_sampled_year), # this connects with historical data
  #     )
  # }) |>
  #   mutate(survey_domain = gsub("HBLL ", "", paste(unique(survey_abbrev), collapse = "/")))

  #TODO: add mpa_trend, ar1_scenario, tag
  prep_full_timeseries <- function(species, sampling_plan, rep_num, sample_dir,
    survey_domain_name = NULL) {

    # Find files matching criteria
    samp_files <- list.files(
      sample_dir,
      pattern = paste0(sp_to_hyphens(species), ".*", sampling_plan, ".*", "rep", sprintf("%03d", rep_num), ".*rds"),
      full.names = TRUE
    )
    n_files <- length(samp_files)

    if (length(samp_files) > 3) {
      warning("Found ", n_files, " files, using first 3 only")
      samp_files <- samp_files[1:3]
    } else {
      message("Found ", length(samp_files), " files")
    }

    # Load lookup table once
    hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))

    # Process files
    sim_dat <- purrr::map_dfr(samp_files, function(f) {
      readRDS(f) |>
        mutate(
          catch_count = observed,
          catch_prop = catch_count / hook_count,
          historical = FALSE
        ) |>
        left_join(hbll_last_sampled_year, by = "survey_abbrev")
    })

    hist_survey_abbrevs <- unique(sim_dat$survey_abbrev)
    hist_dat <- purrr::map_dfr(hist_survey_abbrevs, function(survey_abbrev) {
      readRDS(file.path(hist_dir, paste0(sp_to_hyphens(species), "-", sp_to_hyphens(survey_abbrev), ".rds"))) |>
        mutate(
          historical = TRUE,
          replicate = rep_num
        ) |>
        select(ssid, survey_abbrev, year, X, Y, fishing_event_id, restricted, historical,
          catch_count, catch_prop, hook_count, adjusted_hook_count, replicate, last_sampled_year)
    })

    combined <- bind_rows(sim_dat, hist_dat) |>
      mutate(
          year_post_imp = ifelse(historical, 0, year), # post implementation year
          year = ifelse(historical, year, last_sampled_year + year), # calendar year
          fyear = factor(ifelse(historical, year, last_sampled_year)) # match with end of historical data
      )
    # Create survey domain label
    if (is.null(survey_domain_name)) {
      survey_domain_name <- gsub("HBLL ", "", paste(unique(combined$survey_abbrev), collapse = "/"))
    }

    combined |>
      mutate(survey_domain = survey_domain_name)
  }

  # DATA STRUCTURE CHECKING ------------------------------------------------------
  # This is what goes into the model fitting
  check_data_prep <- prep_full_timeseries(sp, sampling_plan, rep_num = 1, sample_dir)

  # Check that years are correct
  p1 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = year, y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "Year", y = "Catch proportion", colour = "Historical")

  # Check that fyear is correct
  p2 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = as.numeric(fyear), y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "fyear", y = "Catch proportion", colour = "Historical")

  # Check that year_post_imp is correct
  p3 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = year_post_imp, y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "Post-implementation year", y = "Catch proportion", colour = "Historical")
  (p1 / p2 / p3) + plot_annotation(title = "Data alignment / simulation + historical comparison") +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  # ---


  # MONITORING MODEL SETUP -------------------------------------------------------
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

    weights <- dat$hook_count
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
    fit
  }

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


  # samp_files <- list.files(file.path(sample_dir), full.names = TRUE)
  # hist_files <- list.files(hist_dir, pattern = paste0(sp_to_hyphens(species), ".*"), full.names = TRUE)

  # Redefine replicates and eval_years for clarity
  # USE_PARALLEL <- FALSE
  # replicates <- 1
  # eval_years <- c(2034)

  USE_PARALLEL <- TRUE
  if (Sys.info()['user'] %in% c("dunic", "anderson")) {
    USE_PARALLEL <- TRUE
    N_WORKERS <- 80
  }
  if (Sys.info()['user'] %in% c("jillian", "jilliandunic")) {
    USE_PARALLEL <- TRUE
    N_WORKERS <- ifelse(Sys.info()['user'] == "jillian", 10, 8)
  }
  # replicates <- 1:20
  eval_years <- c(2030, 2034, 2038, 2042, 2046)

  SAVE_FITS <- TRUE
  SAVE_FITS <- FALSE

  # eval_years <- c(2030, 2034, 2038, 2042, 2046)

  # Set up parallel backend
  if (USE_PARALLEL) {
    future::plan(future::multisession, workers = N_WORKERS)
    message("Using ", N_WORKERS, " parallel workers (multisession)")
  } else {
    future::plan(future::sequential)
    message("Using sequential processing")
  }


  # Process each replicate in parallel (each rep does all eval years)
  results_parallel <- furrr::future_walk(
  # results_parallel <- purrr::walk(
    replicates,
    function(rep_num) {
      # Fit all eval years for this replicate
      rep_results <- purrr::map_dfr(eval_years, function(eval_yr) {
        message("Fitting model for ", sp,
          "\n - lambda ", round(lambda_val, 4),
          "\n - sampling plan ", sampling_plan,
          "\n - replicate ", rep_num,
          "\n - eval year ", eval_yr)
        message("Rep ", rep_num, ", eval year ", eval_yr)

        # Create output filename
        out_fname <- paste0(
          sp_to_hyphens(sp),
          "-", round(lambda_val, 4),
          "-", sampling_plan,
          "-rep", sprintf("%03d", rep_num),
          "_", "eval", eval_yr, ".rds")

        # Prepare full timeseries
        cdat <- prep_full_timeseries(sp, sampling_plan, rep_num, sample_dir) |>
          filter(year <= eval_yr)

        fit_path <- file.path(fit_dir, out_fname)

        if (!file.exists(fit_path) || SAVE_FITS) {
          # Fit model
          message('Fitting model for eval year ', eval_yr, '; replicate ', rep_num)
          fit <- fit_simulation(
            dat = cdat,
            formula = catch_prop ~ 0 + fyear + restricted:year_post_imp,
            spatial = "on",
            spatiotemporal = "iid",
            family = betabinomial(link = "cloglog"),
            cutoff = 10,
            control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
            silent = FALSE
          )

          saveRDS(fit, fit_path)
        } else {
          fit <- readRDS(fit_path)
        }

        message('Extracting estimates')
        trend_est <- extract_trend_estimate(fit, trend_param = "restricted:year_post_imp")
        re_pars <- extract_re_pars(fit)
        result_row <- tibble(
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
        saveRDS(result_row, file.path(results_dir, out_fname))

        return(result_row)
      })  # End of purrr::map_dfr

    },  # End of function(rep_num)
    .options = furrr::furrr_options(
      seed = TRUE,
      globals = c("fit_simulation", "prep_full_timeseries", "extract_trend_estimate",
                  "extract_re_pars", "sp_to_hyphens", "eval_years", "sample_dir",
                  "hist_files", "species", "fit_dir", "results_dir", "summarise_sanity",
                  "sp", "lambda_val", "sampling_plan", "SAVE_FITS", "hist_dir"),
      packages = c("dplyr", "sdmTMB", "purrr")
    )
  )  # End of furrr::future_map_dfr
  # })
} # end of run_ssf

# sp_list <- c(
#   "yelloweye rockfish", "lingcod", "quillback rockfish",
#   "canary rockfish", "silvergray rockfish",
#   "north pacific spiny dogfish", "pacific halibut"
# )
# mpa_rate_list <- c("5%", "10%", "25%", "50%")

sp_list <- c("yelloweye rockfish")
# mpa_rate_list <- c("10%")
mpa_rate_list <- c("25%")

replicates <- 1:25 # FIXME not used in function right now

stop('stop here', call. = FALSE)

tag_list <- c("svc-rates-sd-0.01", "no-svc-rates")
for (sp in sp_list) {
  for (mpa_rate in mpa_rate_list) {
    for (tag in tag_list) {
      run_ssf(sp, mpa_rate, tag = tag)
    }
  }
}

check_data_prep <- prep_full_timeseries("yelloweye rockfish",
  sampling_plan = "status_quo", rep_num = 1,
  sample_dir = here::here("data-generated", "3-sampled-data-svc-rates-sd-0.01"))

# Check that years are correct
p1 <- ggplot(data = check_data_prep) +
  geom_point(aes(x = year, y = catch_prop, colour = factor(historical)), shape = 21) +
  scale_colour_manual(values = c("orange", "dodgerblue")) +
  facet_wrap(~ survey_abbrev, scales = "free_y") +
  labs(x = "Year", y = "Catch proportion", colour = "Historical")

# Check that fyear is correct
p2 <- ggplot(data = check_data_prep) +
  geom_point(aes(x = as.numeric(fyear), y = catch_prop, colour = factor(historical)), shape = 21) +
  scale_colour_manual(values = c("orange", "dodgerblue")) +
  facet_wrap(~ survey_abbrev, scales = "free_y") +
  labs(x = "fyear", y = "Catch proportion", colour = "Historical")

# Check that year_post_imp is correct
p3 <- ggplot(data = check_data_prep) +
  geom_point(aes(x = year_post_imp, y = catch_prop, colour = factor(historical)), shape = 21) +
  scale_colour_manual(values = c("orange", "dodgerblue")) +
  facet_wrap(~ survey_abbrev, scales = "free_y") +
  labs(x = "Post-implementation year", y = "Catch proportion", colour = "Historical")
(p1 / p2 / p3) + plot_annotation(title = "Data alignment / simulation + historical comparison") +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
# ---

meep()

result_check <- readRDS(file.path(results_dir,
  paste0(sp_to_hyphens(sp),
    "-", round(lambda_val, 4),
    "-", sampling_plan,
    "-rep", sprintf("%03d", 2),
    "_", "eval", 2046, ".rds"
  )
))
glimpse(result_check)

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

# fit_dir <- here::here("data-generated", "test-fits-svc-rates")
# results_dir <- here::here("data-generated", "test-results-svc-rates")
# results_dir <- here::here("data-generated", "test-results-no-svc-rates")
# eval_years <- c(2030, 2034, 2038, 2042, 2046)

# fit_files <- list.files(fit_dir)
# test_fit <- readRDS(file.path(fit_dir, fit_files[2]))
# # sdmTMB::sanity(test_fit)
# results_files <- list.files(results_dir)
# test_res <- readRDS(file.path(results_dir, results_files[2]))

# lambda_val <- readRDS(file.path("data-generated", "recovery-rates-lambda.rds")) |>
#   filter(species == .env$species)
#   filter(case == "50% rate") |>
#   pull(lambda)

message("Reading results from ", basename(results_dir),
  "\n - lambda: ", round(lambda_val, 4),
  "\n - sampling plan: ", sampling_plan,
  "\n - replicate: ", min(replicates), " to ", max(replicates),
  "\n - eval year: ", paste(eval_years, collapse = ", "))

results_files <- list.files(results_dir,
  pattern = paste0(sp_to_hyphens(sp),
  "-", round(lambda_val, 4),
  "-", sampling_plan,
   "-rep[0-9]{3}_eval[0-9]{4}\\.rds$"),
  full.names = TRUE)

all_fitted_results <- purrr::map_dfr(results_files, function(f) {
  readRDS(f)
}) |>
  mutate(mpa_trend = log(lambda_val),
         plan = sampling_plan,
         tag = tag,
         species = sp
         )

# combo <- c("species", "mpa_trend", "tag", "plan", "eval_year")

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

(p1 / b1 / c1) + plot_annotation(title = tag)


