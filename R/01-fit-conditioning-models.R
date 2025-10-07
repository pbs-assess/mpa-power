# =============================================================================
# Fit conditioning models
# =============================================================================

# Load setup and functions
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

devtools::load_all("~/R_DFO/sdmTMB") # need betabinomial branch

library(tidyr)
library(patchwork)
library(digest)

# TODO: remove? no figures needed in this script?
# fig_dir <- here::here("draft-figures", "diagnostics")
# dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Setup directories
fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

# TODO: Move to sim script
# sim_cache <- here::here("data-generated", "sim-cache")
# dir.create(sim_cache, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
# Grids
hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") #|>
  # depth not currently used in sim
  # mutate(depth_mean = mean(log(depth_m), na.rm = TRUE),
  #        depth_sd = sd(log(depth_m), na.rm = TRUE),
  #        depth_scaled = (log(depth_m) - depth_mean[1]) / depth_sd[1])

hbll_grid_poly <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S")

hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
simple_mpa <- readRDS(here::here("data-generated", "spatial", "simple-mpa.rds"))
# comm_ll_activity_status <- readRDS(here::here("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))
# mpa_shape_simplified <- comm_ll_activity_status |> st_simplify(dTolerance = 100)

# HBLL species data
# sp <- sp_to_hyphens("north pacific spiny dogfish")
# sp <- sp_to_hyphens("lingcod")
# sp <- sp_to_hyphens("quillback rockfish")
sp <- sp_to_hyphens("yelloweye rockfish")

# sp_list <- c("yelloweye rockfish", "north pacific spiny dogfish", "lingcod", "quillback rockfish")
# for (sp in sp_to_hyphens(sp_list)) {
message(paste0("Fitting conditioning models for ", sp))
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) |>
  mutate(
    obs_id = factor(row_number()),
    catch_prop = catch_count / hook_count,
    log_depth = log(depth_m)
  )

historical <- sp_dat |>
  mutate(x = X * 1000, y = Y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 3156) |>
  st_join(simple_mpa |> st_transform(crs = 3156), join = st_within) |>
  mutate(restricted = ifelse(is.na(uid), 0, 1)) |>
  st_join(hbll_grid_poly |> select(block_id, grouping_code) |> st_transform(crs = 3156), join = st_within) |>
  st_drop_geometry() |>
  select(ssid, survey_abbrev, year, species_common_name, fishing_event_id, latitude, longitude, X, Y,
    block_id, fe_grouping_code = grouping_code.x, grouping_code = grouping_code.y,
    depth_m, offset, hook_count,
    catch_count, restricted) |>
  mutate(after = 0) |>
  left_join(hbll_allocations, by = c("survey_abbrev", "grouping_code", "ssid" = "survey_series_id"))

# Prepare data and meshes
d_IN <- sp_dat |> filter(survey_abbrev == "HBLL INS N")
d_IN$weights <- d_IN$hook_count / mean(d_IN$hook_count)
mesh_IN <- local(make_mesh(d_IN, xy_cols = c("X", "Y"), cutoff = 10))
d_OS <- sp_dat |> filter(survey_abbrev == "HBLL OUT S")
d_OS$weights <- d_OS$hook_count / mean(d_OS$hook_count)
mesh_OS <- local(make_mesh(d_OS, xy_cols = c("X", "Y"), cutoff = 10))
d_ON <- sp_dat |> filter(survey_abbrev == "HBLL OUT N")
d_ON$weights <- d_ON$hook_count / mean(d_ON$hook_count)
mesh_ON <- local(make_mesh(d_ON, xy_cols = c("X", "Y"), cutoff = 10))

check_cache <- TRUE
silent <- FALSE

# TODO: deal with random fields that collapse to 0

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
  weights = d_ON$hook_count,
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
  weights = d_OS$hook_count,
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
  weights = d_IN$hook_count,
  check_cache = check_cache,
  silent = silent
)

# TODO: turn of random fields that fail - currently this just turns off spatiotemporal
# TODO - add sanity check as tag in filename so that we can reload the appropriate cache?
# and then spatial field sequentially if sanity checks keep failing.
# Function to check sanity and refit if needed
refit_if_failed <- function(fit, survey_name, sp, fit_dir) {
  # Use model-specific values, not global ones
  current_sprf <- fit$spatial
  current_strf <- fit$spatiotemporal

  sanity_check <- all(unlist(sdmTMB::sanity(fit)))
  if (!sanity_check) {
    current_strf <- "off"
    message(paste0("Sanity check failed for ", survey_name, ", refitting with spatiotemporal = 'off'"))
    fit <- fit_cached_sdmTMB(
      model_tag = paste0(sp, "-", survey_name, "-betabinomial-", current_sprf, "-", current_strf),
      fit_dir = fit_dir,
      update_from = fit,
      spatiotemporal = current_strf
    )

    sanity_check2 <- all(unlist(sdmTMB::sanity(fit)))
    if (!sanity_check2) {
      current_sprf <- "off"
      message(paste0("Second sanity check failed for", survey_name, ", refitting with spatial = 'off'"))
      fit <- fit_cached_sdmTMB(
        model_tag = paste0(sp, "-", survey_name, "-betabinomial-", current_sprf, "-", current_strf),
        fit_dir = fit_dir,
        update_from = fit,
        spatial = current_sprf
      )
    }
  }

  return(fit)
}

# Apply sanity checks and refitting to all models
fit_ON <- refit_if_failed(fit_ON, "HBLL-OUT-N", sp, fit_dir)
fit_OS <- refit_if_failed(fit_OS, "HBLL-OUT-S", sp, fit_dir)
fit_IN <- refit_if_failed(fit_IN, "HBLL-INS-N", sp, fit_dir)

meep()