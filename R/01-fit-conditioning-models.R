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
USE_PARALLEL <- FALSE
N_WORKERS <- NULL

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

  restricted_df <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
    st_drop_geometry() |>
    select(X, Y, uid, restricted)

  sp <- sp_to_hyphens(sp_name)
  message(paste0("Fitting conditioning models for ", sp))
  sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

  sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
    filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
    prep_hbll_data(bait_counts = bait_counts, restricted_df = restricted_df) |>
    mutate(
      obs_id = factor(row_number()),
      catch_prop = catch_count / hook_count,
      log_depth = log(depth_m)
    )

  d_ON <- filter(sp_dat, survey_abbrev == "HBLL OUT N")
  d_OS <- filter(sp_dat, survey_abbrev == "HBLL OUT S")
  d_IN <- filter(sp_dat, survey_abbrev == "HBLL INS N")

  mesh_ON <- local(make_mesh(d_ON, xy_cols = c("X", "Y"), cutoff = 10))
  mesh_OS <- local(make_mesh(d_OS, xy_cols = c("X", "Y"), cutoff = 10))
  mesh_IN <- local(make_mesh(d_IN, xy_cols = c("X", "Y"), cutoff = 10))

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
  fit_ON <- fit_cached_sdmTMB(
    model_tag = paste0(sp, "-HBLL-OUT-N-betabinomial-", sprf, "-", strf),
    fit_dir = fit_dir,
    data = d_ON,
    formula = catch_prop ~ 0 + fyear,
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

  fit_OS <- fit_cached_sdmTMB(
    model_tag = paste0(sp, "-HBLL-OUT-S-betabinomial-", sprf, "-", strf),
    fit_dir = fit_dir,
    data = d_OS,
    formula = catch_prop ~ 0 + fyear,
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

  fit_IN <- fit_cached_sdmTMB(
    model_tag = paste0(sp, "-HBLL-INS-N-betabinomial-", sprf, "-", strf),
    fit_dir = fit_dir,
    data = d_IN,
    formula = catch_prop ~ 0 + fyear,
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

  message(paste0("Completed: ", sp_name))
  return(invisible(list(
    fit_ON = fit_ON,
    fit_OS = fit_OS,
    fit_IN = fit_IN)))
}

# -----------------------------------------------------------------------------
# Run fits
# -----------------------------------------------------------------------------

# Single species (for testing)
# test <- fit_species("yelloweye rockfish", save_cleaned_data = FALSE)

# look at new species fits:
# test0 <- fit_species("yelloweye rockfish", save_cleaned_data = F)
# test <- fit_species("pacific cod", save_cleaned_data = F)
# test[[1]] |> sanity()
# test[[2]] |> sanity()
# test[[3]] |> sanity()

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
