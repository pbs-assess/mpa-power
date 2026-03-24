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
fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
cleaned_data_dir <- here::here("data-generated", "cleaned-species-data")
dir.create(cleaned_data_dir, recursive = TRUE, showWarnings = FALSE)
mesh_dir <- here::here("data-generated", "mesh-cache")
dir.create(mesh_dir, recursive = TRUE, showWarnings = FALSE)
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

shared_hbll_mesh <- readRDS(shared_mesh_file)
message("Using shared HBLL mesh: ", shared_mesh_file)

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))

# Fitting parameters
check_cache <- TRUE
silent <- TRUE

# Fit models for a single species
fit_species <- function(sp_name, check_cache = TRUE, silent = FALSE,
                        save_cleaned_data = TRUE,
                        .options = furrr::furrr_options(seed = TRUE)) {

  historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
    st_drop_geometry() |>
    select(X, Y, uid, restricted)

  sp <- sp_to_hyphens(sp_name)
  message(paste0("Fitting conditioning models for ", sp))
  sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

  sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
    filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
    prep_hbll_data(bait_counts = bait_counts, restricted_df = historical_locations) |>
    mutate(
      obs_id = factor(row_number()),
      catch_prop = catch_count / hook_count,
      log_depth = log(depth_m),
      last_sampled_year = max(year),
      year_covariate = 0,
      historical = TRUE
    )

  d_ON <- filter(sp_dat, survey_abbrev == "HBLL OUT N")
  d_OS <- filter(sp_dat, survey_abbrev == "HBLL OUT S")
  d_IN <- filter(sp_dat, survey_abbrev == "HBLL INS N")

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

  mesh_ON <- shared_hbll_mesh
  mesh_OS <- shared_hbll_mesh
  mesh_IN <- shared_hbll_mesh

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
  conditioning_formula <- catch_prop ~ 0 + fyear + restricted
  fit_ON <- if (val_ON$passed) {
    fit_cached_sdmTMB(
      model_tag = paste0(sp, "-HBLL-OUT-N-betabinomial-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_ON,
      formula = conditioning_formula,
      mesh = mesh_ON,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
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
      model_tag = paste0(sp, "-HBLL-OUT-S-betabinomial-restricted-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_OS,
      formula = conditioning_formula,
      mesh = mesh_OS,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
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
      model_tag = paste0(sp, "-HBLL-INS-N-betabinomial-restricted-", sprf, "-", strf),
      fit_dir = fit_dir,
      data = d_IN,
      formula = conditioning_formula,
      mesh = mesh_IN,
      family = betabinomial(link = "cloglog"),
      spatial = sprf,
      spatiotemporal = strf,
      time = "year",
      anisotropy = FALSE,
      weights = d_IN$adjusted_hook_count,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
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
# test <- fit_species("yelloweye rockfish", save_cleaned_data = FALSE)

# look at new species fits:
# test0 <- fit_species("yelloweye rockfish", save_cleaned_data = F)
# test <- fit_species("silvergray rockfish", save_cleaned_data = F)
# test[[1]] |> sanity()
# test[[2]] |> sanity()
# test[[3]] |> sanity()

# hbll_grid <- gfdata::load_survey_blocks(type = "XY")
# hbll_grid_sf <- gfdata::load_survey_blocks(type = "polygon") |>
#   select(survey_abbrev, block_id, geometry) |>
#   filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"))
# pON <- predict_hbll(test[[1]], grid = hbll_grid)
# pOS <- predict_hbll(test[[2]], grid = hbll_grid)

# dat <- left_join(
#   hbll_grid_sf,
#   bind_rows(
#     pON |> filter(year == 2023) |> select(survey_abbrev, block_id, year, est),
#     pOS |> filter(year == 2024) |> select(survey_abbrev, block_id, year, est)
#     )
# )

# display_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa-500m.rds"))
# ggplot(data = dat |> rotate_sf()) +
#   geom_sf(aes(fill = exp(est) * 100, colour = exp(est) * 100)) +
#   geom_sf(data = display_mpa |> rotate_sf(),
#     fill = "grey50", colour = "grey50", alpha = 0.5) +
#   viridis::scale_fill_viridis(trans = ggsidekick::fourth_root_power_trans()) +
#   viridis::scale_colour_viridis(trans = ggsidekick::fourth_root_power_trans()) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_sf(), buffer = 0)

# Species list
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

# Get process error parameters from conditioning models
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

saveRDS(ar1_estimates, here::here("data-generated", "ar1-parameters.rds"))
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
saveRDS(fit_characteristics, here::here("data-generated", "fit_characteristics.rds"))
