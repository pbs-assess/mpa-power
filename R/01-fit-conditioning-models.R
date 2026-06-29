# =============================================================================
# Fit conditioning models
# =============================================================================

# Load setup and functions
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# devtools::load_all("~/R_DFO/sdmTMB") # need betabinomial branch

library(tidyr)
library(patchwork)
library(digest)
library(purrr)

# =============================================================================
# Configuration
# =============================================================================
USE_PARALLEL <- TRUE
N_WORKERS <- 8

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 40 #NULL
}

if (Sys.info()['user'] == "jilliandunic") {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 8
}

# Setup directories
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cleaned_data_dir, recursive = TRUE, showWarnings = FALSE)

shared_mesh_species <- "yelloweye rockfish"
shared_mesh_file <- file.path(
  mesh_dir,
  paste0("HBLL-combined-mesh.rds")
)

if (!file.exists(shared_mesh_file)) {
  stop(
    "Shared HBLL mesh file not found: ", shared_mesh_file, ". ",
    "Build it first with R/build-shared-hbll-mesh.R",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
hbll_allocations <- readRDS(hbll_allocations_file)
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))
survey_set_depths <- readRDS(here::here("data-generated", "hbll-dem-survey-depths.rds")) |>
  group_by(survey_abbrev, fishing_event_id) |>
  mutate(dem_depth = mean(c(depth_start, depth_end), na.rm = TRUE)) |>
  ungroup() |>
  select(survey_abbrev, fishing_event_id, dem_depth)

# Fitting parameters
check_cache <- TRUE
silent <- FALSE

# get_unscaled_rho <- function(rho_time) qlogis((rho_time + 1) / 2)
# get_unscaled_rho(0.2)

# Fit models for a single species
fit_species <- function(sp_name, check_cache = FALSE, silent = FALSE,
                        save_cleaned_data = TRUE,
                        # INS_prior = NULL,
                        # INS_control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
                        .options = furrr::furrr_options(seed = TRUE)) {

  # get_unscaled_rho <- function(rho_time) qlogis((rho_time + 1) / 2)

  historical_locations <- readRDS(historical_locations_file) |>
    st_drop_geometry() |>
    select(X, Y, uid, fishing_event_id, restricted)

  sp <- sp_to_hyphens(sp_name)
  message(paste0("Fitting conditioning models for ", sp))
  sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

  sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
    filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
    prep_hbll_data(bait_counts = bait_counts, restricted_df = historical_locations) |>
    distinct(survey_abbrev, fishing_event_id, .keep_all = TRUE) |> # one duplicate with different deployment/retrieval time
    left_join(survey_set_depths, by = c("survey_abbrev", "fishing_event_id")) |>
    mutate(
      obs_id = factor(row_number()),
      catch_prop = catch_count / hook_count,
      log_depth = log(dem_depth), # use DEM depth to match prediction grid DEM depth
      last_sampled_year = max(year),
      year_covariate = 0,
      historical = TRUE
    ) |>
    drop_na(log_depth)

  d_ON <- filter(sp_dat, survey_abbrev == "HBLL OUT N")
  d_OS <- filter(sp_dat, survey_abbrev == "HBLL OUT S")
  d_IN <- filter(sp_dat, survey_abbrev == "HBLL INS N")

  # extra_time_ON <- setdiff(seq(min(d_ON$year), max(d_ON$year)), unique(d_ON$year))
  # extra_time_OS <- setdiff(seq(min(d_OS$year), max(d_OS$year)), unique(d_OS$year))
  # extra_time_IN <- setdiff(seq(min(d_IN$year), max(d_IN$year)), unique(d_IN$year))

  # Validate each dataset
  val_ON <- validate_hbll_data(d_ON, sp_name, "HBLL OUT N")
  val_OS <- validate_hbll_data(d_OS, sp_name, "HBLL OUT S")
  val_IN <- validate_hbll_data(d_IN, sp_name, "HBLL INS N")

  # Report validation results
  for (val in list(val_ON, val_OS, val_IN)) {
    if (!val$passed) {
      warning(sprintf("%s - %s: Data validation failed: %s",
                     val$species, val$survey,
                     paste(val$issues, collapse = "; ")),
              call. = FALSE)
    } else {
      years_with_pos <- val$n_years - val$years_with_all_zeroes
      message(sprintf("  %s: %d / %d (%.1f%%) pos sets; %d years with pos sets",
                     val$survey, val$n_positive, val$n_obs,
                     val$prop_pos * 100, years_with_pos))
    }
  }

  # Cache vertices to make sure these are consistent across platforms
  f <- file.path(mesh_dir, paste0(sp, "-HBLL-OUT-N.rds"))
  if (!file.exists(f)) {
    mesh_ON <- local(make_mesh(d_ON, xy_cols = c("X", "Y"), cutoff = 10))
    saveRDS(mesh_ON, file.path(mesh_dir, paste0(sp, "-HBLL-OUT-N.rds")))
  } else {
    message("Reading cached mesh")
    mesh_ON <- readRDS(f)
  }

  f <- file.path(mesh_dir, paste0(sp, "-HBLL-OUT-S.rds"))
  if (!file.exists(f)) {
    mesh_OS <- local(make_mesh(d_OS, xy_cols = c("X", "Y"), cutoff = 10))
    saveRDS(mesh_OS, file.path(mesh_dir, paste0(sp, "-HBLL-OUT-S.rds")))
  } else {
    message("Reading cached mesh")
    mesh_OS <- readRDS(f)
  }

  f <- file.path(mesh_dir, paste0(sp, "-HBLL-INS-N.rds"))
  if (!file.exists(f)) {
    mesh_IN <- local(make_mesh(d_IN, xy_cols = c("X", "Y"), cutoff = 10))
    saveRDS(mesh_IN, file.path(mesh_dir, paste0(sp, "-HBLL-INS-N.rds")))
  } else {
    message("Reading cached mesh")
    mesh_IN <- readRDS(f)
  }

  # Save cleaned datasets
  if (save_cleaned_data) {
    saveRDS(d_ON, file.path(cleaned_data_dir, paste0(sp, "-HBLL-OUT-N.rds")))
    saveRDS(d_OS, file.path(cleaned_data_dir, paste0(sp, "-HBLL-OUT-S.rds")))
    saveRDS(d_IN, file.path(cleaned_data_dir, paste0(sp, "-HBLL-INS-N.rds")))
    message("  Saved cleaned data for ", sp)
  }

  # Beta binomial ----------------------------------------------------------------
  sprf <- "on"
  strf <- "iid"
  # conditioning_formula <- catch_prop ~ 0 + fyear + restricted + log_depth + I(log_depth^2)
  conditioning_formula <- catch_prop ~ 0 + fyear + restricted + log_depth + I(log_depth^2)

  fit_ON <- if (val_ON$passed) {
    fit_cached_sdmTMB(
      model_tag = paste0(sp, "-HBLL-OUT-N-betabinomial-restricted-depth-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_ON,
      formula = conditioning_formula,
      mesh = mesh_ON,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
      # time_varying = ~ 1,
      # time_varying_type = "ar1",
      # extra_time = extra_time_ON,
      anisotropy = FALSE,
      weights = d_ON$adjusted_hook_count,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
      check_cache = check_cache,
      silent = silent
    )
  } else {
    message("  Skipping HBLL OUT N fit due to data quality issues")
    NULL
  }

  fit_OS <- if (val_OS$passed) {
    fit_cached_sdmTMB(
      model_tag = paste0(sp, "-HBLL-OUT-S-betabinomial-restricted-depth-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_OS,
      formula = conditioning_formula,
      mesh = mesh_OS,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
      # time_varying = ~ 1,
      # time_varying_type = "ar1",
      # extra_time = extra_time_OS,
      anisotropy = FALSE,
      weights = d_OS$adjusted_hook_count,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
      check_cache = check_cache,
      silent = silent
    )
  } else {
    message("  Skipping HBLL OUT S fit due to data quality issues")
    NULL
  }

  fit_IN <- if (val_IN$passed) {
    fit_cached_sdmTMB(
      model_tag = paste0(sp, "-HBLL-INS-N-betabinomial-restricted-depth-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_IN,
      formula = conditioning_formula,
      mesh = mesh_IN,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
      # time_varying = ~ 1,
      # time_varying_type = "ar1",
      # extra_time = extra_time_IN,
      anisotropy = FALSE,
      weights = d_IN$adjusted_hook_count,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
      # control = INS_control,
      # priors = INS_prior,
      check_cache = check_cache,
      silent = silent
    )
  } else {
    message("  Skipping HBLL INS N fit due to data quality issues")
    NULL
  }

  message(paste0("Completed: ", sp_name))
  return(invisible(list(
    fit_ON = fit_ON,
    fit_OS = fit_OS,
    fit_IN = fit_IN,
    validation_ON = val_ON,
    validation_OS = val_OS,
    validation_IN = val_IN
  )))
}

assert_fits_have_omega_s <- function(all_fits_flat) {
  fitted_models <- all_fits_flat |>
    purrr::keep(~!is.null(.x) && inherits(.x, "sdmTMB"))

  omega_check <- purrr::map_dfr(fitted_models, function(fit) {
    osp <- one_sample_posterior(fit)
    omega_draw <- osp[grepl("omega_s", names(osp))]
    tibble(
      species = unique(fit$data$species_common_name),
      survey_abbrev = unique(fit$data$survey_abbrev),
      n_omega_s = length(omega_draw),
      has_omega_s = length(omega_draw) > 0
    )
  })

  missing_omega <- omega_check |>
    filter(!has_omega_s)

  if (nrow(missing_omega) > 0) {
    missing_labels <- paste(
      missing_omega$species,
      "(",
      missing_omega$survey_abbrev,
      ")"
    )

    stop(
      "Conditioning fit validation failed: omega_s was empty for ",
      nrow(missing_omega),
      " fit(s): ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  message("Validated omega_s posterior draws for ", nrow(omega_check), " conditioning fit(s)")
}

# -----------------------------------------------------------------------------
# Run fits
# -----------------------------------------------------------------------------
# Single species (for testing)
# Dogfish OUT S had one year estimate not meet the gradient threshold but pretty close (gradien = 0.00381)
# test_fit <- fit_species("north pacific spiny dogfish", save_cleaned_data = FALSE, check_cache = TRUE)
# tibble::tibble(param = names(test_fit$fit_OS$model$par), gradient = test_fit$fit_OS$gradients) |>
#   dplyr::arrange(dplyr::desc(abs(gradient)))

# test_fit[[3]] |> sanity()
# test_fit <- fit_species("yelloweye rockfish", save_cleaned_data = FALSE)

# # look at new species fits:
# test0 <- fit_species("yelloweye rockfish", save_cleaned_data = F, check_cache = TRUE)
# # test <- fit_species("silvergray rockfish", save_cleaned_data = F)
# test_fit[[1]] |> sanity()
# test_fit[[2]] |> sanity()
# test_fit[[3]] |> sanity()

# test_fit[[3]]

# Attempted the below; but it ended up being rather fickle, and would likely make
# the final estimation model also very fickle, so including depth, but going back
# to the two step approach.
# - Rather than fitting independent year effects and fitting a post hoc ARIMA
#   model to estimate AR1 parameters, we directly model the AR1 on the
#   time-varying intercept in the conditioning model. This makes the conditioning
#   model consistent with the simulation and evaluation models (yes?), and
#   jointly estimates these parameters at the same time as the others.
# - Including both the time_varying year intercept and IID spatiotemporal
#   random fields creates competition to explain the year-to-year variation and
#   increases the chance of sigma_V collapsing. A weakly informative prior on
#   the sigma_V we can keep it from collapsing.
# - We used sigma_V = gamma_cv(0.2, 0.5) prior for yelloweye and lingcod INS N.
#   Parameterisation: gamma_cv(mean, CV); CV = sd / mean, which seems reasonable
#   given the original post hoc arima model estimates:
# og_ar_estimates <- readRDS(here::here("data-generated", "ar1-parameters-resdoc-nodepth.rds"))
# og_ar_estimates
# og_ar_estimates |> filter(species %in% c("yelloweye rockfish", "lingcod"), survey_abbrev == "HBLL INS N")
# - In the case of yelloweye HBLL INS N, the rho_time was unidentifiable even
#   with a sigma_V prior, meaning there was little/no temporal variation for the
#   model to estimate autocorrelation from. The post-hoc rho for YE was -0.08
#   (near-zero), so we fixed rho_time_unscaled at 0.
# Note: rho_time_unscaled is on the unconstrained scale, but rho is constrained
# to -1 and 1.
#   See: sdmTMB/src/utils.h:280 (minus_one_to_one), sdmTMB/src/sdmTMB.cpp:710
#   transformation is rho = 2 * plogis(rho_time_unscaled) - 1
#   inverse: rho_time_unscaled = qlogis((rho + 1) / 2)
# # The transformaation is rho_time = 2 * plogis(rho_time_unscaled) - 1
# qlogis((0 + 1) / 2)    # rho = 0  → rho_unscaled = 0
# # example of how we convert if we were fixing at a value other than 0
# rho = 0.32
# qlogis((rho + 1) / 2)

# ye <- fit_species("yelloweye rockfish")
# ye <- fit_species("yelloweye rockfish", save_cleaned_data = TRUE,
#   INS_control = sdmTMBcontrol(collapse_spatial_variance = TRUE,
#     map = list(rho_time_unscaled = factor(matrix(NA, nrow = 1, ncol = 1))),
#     start = list(rho_time_unscaled = matrix(0, nrow = 1, ncol = 1)),
#   ),
#   INS_prior = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5))
# )
# qb <- fit_species("quillback rockfish", save_cleaned_data = TRUE)
# dog <- fit_species("north pacific spiny dogfish", save_cleaned_data = TRUE)
# lng <- fit_species("lingcod", save_cleaned_data = TRUE,
#   INS_prior = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5)))
# hal <- fit_species("pacific halibut", save_cleaned_data = TRUE)
# can <- fit_species("canary rockfish", save_cleaned_data = TRUE)
# sil <- fit_species("silvergray rockfish", save_cleaned_data = TRUE)

# Species list (needed for the assert omega)
sp_list <- c(
  "yelloweye rockfish",
  "north pacific spiny dogfish",
  "lingcod",
  "quillback rockfish",
  "pacific halibut",
  "canary rockfish",
  "silvergray rockfish"
)

# readRDS(file.path(overlay_dir, "hbll-spp-encounter-rate.rds")) |>
#   select(species_common_name, pos_sets) |>
#   slice(1:20)

# # Setup parallel processing

map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)
# Fit all species
if (USE_PARALLEL) {
  all_fits <- map_fn(sp_list, fit_species,
                     check_cache = check_cache,
                     silent = silent,
                     save_cleaned_data = TRUE,
                     .options = furrr::furrr_options(seed = TRUE))
} else {
  all_fits <- map_fn(sp_list, fit_species,
                     check_cache = check_cache,
                     silent = silent,
                     save_cleaned_data = TRUE)
}

# Reset to sequential
future::plan(future::sequential)

# # Get process error parameters from conditioning models
ar1_estimates <- all_fits |>
  purrr::flatten() |>
  purrr::keep(~!is.null(.x) && inherits(.x, "sdmTMB")) |>
  purrr::map_dfr(function(fit) {
    species <- unique(fit$data$species_common_name)
    survey_abbrev <- unique(fit$data$survey_abbrev)

    year_ests <- fit |>
      get_model_pars() |>
      filter(stringr::str_detect(term, "fyear")) |>
      mutate(year = as.numeric(stringr::str_extract(term, "\\d+"))) |>
      select(year, estimate)
    afit <- arima(year_ests$estimate, order = c(1, 0, 0))
    tibble(
           species = species,
           survey_abbrev = survey_abbrev,
           rho_V = afit$coef[["ar1"]],
           sigma_V = sqrt(afit$sigma2))
  })

saveRDS(ar1_estimates, ar1_parameters_file)
message("Saved AR1 parameters for ", nrow(ar1_estimates), " species × survey combinations")

# Summary of fitting outcomes
all_fits_flat <- all_fits |> purrr::flatten()
n_fits_attempted <- length(sp_list) * 3  # 3 surveys per species
n_fits_succeeded <- sum(purrr::map_lgl(all_fits_flat, ~!is.null(.x) && inherits(.x, "sdmTMB")))
n_fits_failed <- n_fits_attempted - n_fits_succeeded

message("\n=== Fitting Summary ===")
message(sprintf("Attempted: %d species × survey combinations", n_fits_attempted))
message(sprintf("Succeeded: %d", n_fits_succeeded))
message(sprintf("Failed: %d", n_fits_failed))

if (n_fits_failed > 0) {
  message("\nFailed combinations:")
  for (i in seq_along(all_fits)) {
    sp <- sp_list[i]
    fits <- all_fits[[i]]
    for (survey in c("fit_ON", "fit_OS", "fit_IN")) {
      if (is.null(fits[[survey]])) {
        val_name <- gsub(survey, pattern = "fit", replacement = "validation")
        survey_name <- fits[[val_name]]$survey
        message(sprintf("  - %s: %s", sp, survey_name))
      }
    }
  }
}

assert_fits_have_omega_s(all_fits_flat)

# Grab some values for comparison later:
fit_characteristics <- purrr::map_dfr(all_fits_flat, \(x) {
  if (!is.null(x) && inherits(x, "sdmTMB")) {
    p <- tidy(x, "ran_pars")
    sp <- unique(x$data$species_common_name)
    surv <- unique(x$data$survey_abbrev)
    encountered_count_per_year <- sum(x$data$catch_count) / length(unique(x$data$year))
    temp <- group_by(x$data, fishing_event_id) |>
      summarise(encountered = sum(catch_count, na.rm = TRUE) > 1) |> ungroup()
    encountered_rate_avg <- mean(temp$encountered)
    pw <- p |> select(1:2) |> tidyr::pivot_wider(names_from  = term, values_from = estimate)

    temp <- group_by(x$data, fishing_event_id, restricted) |>
      summarise(encountered = sum(catch_count, na.rm = TRUE) > 1) |>
      filter(restricted == 1)
    encountered_rate_restricted <- mean(temp$encountered)

    temp <- filter(x$data, restricted == 1)
    encountered_count_per_year_restricted <- sum(temp$catch_count) / length(unique(temp$year))

    out <- tibble(
      species = sp,
      survey = surv,
      encountered_count_per_year = encountered_count_per_year,
      encountered_rate_avg = encountered_rate_avg,
      encountered_rate_restricted = encountered_rate_restricted,
      encountered_count_per_year_restricted = encountered_count_per_year_restricted
    )
    bind_cols(out, pw)
  }
})
saveRDS(fit_characteristics, fit_characteristics_file)
