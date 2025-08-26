# =============================================================================
# Fit conditioning models
# =============================================================================

# Load setup and functions
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# Setup directories
fit_dir <- here::here("data-generated", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

# From: https://github.com/pbs-assess/gfmpa/blob/9429210b9da5b5044f3afddcf6eaa9cffaec4d40/analysis/sim.Rmd
# - Simulate recovery in restricted areas to assess whether current survey effort is sufficient to detect population recovery

# - Approach:
#     1. fit geostatistical models to observed data
#     2. use parameters from that model to simulate new data with observations at the actual historically observed locations
#     3. when simulating, simulate recovery at some rate within closed areas and a stationary abundance/density outside of closed areas

# - Dimensions that will likely affect the answer:
#     1. species (therefore estimated spatial and spatiotemporal SD, spatial correlation range, observation error)
#     2. rate of 'recovery'
#     3. number of years observed
#     4. whether one assesses all restricted areas together or individually
#     5. **level of fishing (or activity of concern) that occurred before management actions taken**
#     6. size of restricted area?

# Notes:
# - to start no depth because I think there is a lot of uncertainty in the grid depth
#   - TODO: add depth to prediction grid (part of gfdata updates I have going)

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

# HBLL species data
sp <- sp_to_hyphens("yelloweye rockfish")
bait_counts <- readRDS(file.path(synopsis_cache, "bait-counts.rds"))
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets
comm_ll_activity_status <- readRDS(here::here("data-generated", "spatial", "comm-ll-draft-activity-status.rds"))

mpa_shape_simplified <- comm_ll_activity_status |> st_simplify(dTolerance = 100)

sp_dat <- filter(sp_dat0, stringr::str_detect(survey_abbrev, "HBLL")) |>
  filter(survey_abbrev != "HBLL INS S") |> # may as well remove this up here
  prep_hbll_data(bait_counts = bait_counts) |>
  mutate(log_depth = log(depth_m))

sp_dat |> pull(depth_m) |> summary()

combined <- st_intersection(
    st_as_sfc(st_bbox(hbll_grid_poly |> st_transform(crs = st_crs(mpa_shape_simplified)))),
    st_as_sfc(st_bbox(mpa_shape_simplified))
  ) |>
  st_as_sf()

plot_limits_combined <- get_plot_limits(combined, buffer = 1000)

# Fit conditioning models
# ------------------------------------------------------------
{
  .formula <- catch_count ~ 0 + fyear
  .tag <- NULL
  .family <- nbinom2(link = "log")
  .check_cache <- TRUE
}
# {
# .formula <- catch_count ~ 0 + fyear + poly(log_depth, 2) # poly() creates orthogonal polynomials (uncorrelated terms) vs raw polynomials
# .check_cache <- TRUE
# .tag <- "depth"
# .family <- nbinom2(link = "log")
# .spatiotemporal <- "off"
# # .family <- tweedie(link = "log")
# # .tag <- "depth-tweedie"
# }

# {
#   .formula <- catch_count ~ 0 + fyear
#   .tag <- "nb2mixed"
# }

fit_OS <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT S",
  formula = .formula,
  family = .family,
  species = sp,
  spatiotemporal = "iid",
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir,
  tag = .tag,
  check_cache = .check_cache
)
fit_ON <- fit_hbll(dat = sp_dat,
  survey_type = "HBLL OUT N",
  formula = .formula,
  family = .family,
  species = sp,
  spatiotemporal = "iid",
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir,
  tag = .tag,
  check_cache = .check_cache
)
fit_IN <- fit_hbll(dat = sp_dat, # didn't converge with spatiotemporal = "iid"
  survey_type = "HBLL INS N",
  formula = .formula,
  family = .family,
  species = sp,
  spatiotemporal = "off",
  use_extra_time = FALSE,
  time = "year",
  fit_dir = fit_dir,
  tag = .tag,
  check_cache = .check_cache
)
meep()

stop()
# TODO: add sanity checks
# TODO: evaluate and compare conditioning models: see - https://github.com/mis-assess/shrimp_surveydesign_csas/blob/794abdf0d4657dff5ed3316fe876b58afab0dd83/Reproducible_Examples/coastwide-density.R#L157
# fit <- fit_OS
fit <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S.rds"))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
dev.set(2)
plot(r_nb2, title = "no depth")



model_name <- "HBLL-OUT-N-depth"
model_name <- "HBLL-OUT-S-depth"
model_name <- "HBLL-INS-N-depth"
fit <- readRDS(here::here(fit_dir, paste0("yelloweye-rockfish-", model_name, ".rds")))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
# dev.set(3)
plot(r_nb2, title = model_name)

fit <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S-depth.rds"))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
dev.set(3)
plot(r_nb2, title = "depth")




fit <- readRDS(here::here(fit_dir, "yelloweye-rockfish-HBLL-OUT-S-nb2mixed.rds"))
s_nb2 <- simulate(fit, nsim = 500, type = "mle-mvn")
r_nb2 <- dharma_residuals(s_nb2, fit, return_DHARMa = TRUE)
dev.set(3)
plot(r_nb2, title = "depth")

# DHARMa::testResiduals(r_nb2)
# sp_r <- DHARMa::recalculateResiduals(s_nb2, group = fit$data$fyear)
# DHARMa::testSpatialAutocorrelation(r_nb2,
#   x = fit$data$X,
#   y = fit$data$Y)
# DHARMa::testZeroInflation(r_nb2)
