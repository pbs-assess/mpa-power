# =============================================================================
# Fit sampled data
# =============================================================================

source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))

library(purrr)
library(progressr)

# =============================================================================
# Configuration
# =============================================================================

cleaned_data_dir <- here::here("data-generated", "cleaned-species-data")
sample_dir <- here::here("data-generated", "sampled-data")
results_dir <- here::here("data-generated", "power-results")
hist_path <- here::here("data-generated", "historical-data-processed")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- TRUE
N_WORKERS <- 8 #NULL
N_REPLICATES <- 50

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 70
  N_REPLICATES <- 100
}

if (Sys.info()['user'] %in% c("jillian", "jilliandunic")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- ifelse(Sys.info()['user'] == "jillian", 10, 8)
  N_REPLICATES <- 100
}

FORMULA <- catch_prop ~ 0 + fyear + restricted:year_covariate
EVALUATION_YEARS <- c(2030, 2034, 2038, 2042, 2046)

hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))

# # Testing
# sample_summary <- readRDS(file.path(sample_dir,  "sampling-summary.rds"))

# # species <- "lingcod"
# species <- "yelloweye rockfish"
# ar1_scenarios <- c("fitted_AR1")

# # time_scenarios <- c("twentyfive_years")
# plans <- c(
#   "status quo",
#   "MPAs at 5 year intervals"#,
#   # "status quo + 20% effort"
# )

# # f <- list.files(file.path(sample_dir, sp_to_hyphens(species)))

# sp_files <- filter(sample_summary,
#   species %in% .env$species,
#   ar1_scenario %in% .env$ar1_scenarios,
#   # time_scenario %in% time_scenarios,
#   plan %in% .env$plans
# ) |>
#   pull(file)

# # Create a cache environment (can reuse for multiple calls)
# hist_cache <- new.env(parent = emptyenv())
# # Load historical data
# hist_data <- purrr::map_dfr(c("HBLL OUT N"), function(survey_abbrev) {
#     get_hist_data(
#     species = species,
#     survey_abbrev = survey_abbrev,  # or whatever survey you're testing
#     hist_path = hist_path,
#     cache_env = hist_cache
#   )
# })
# # Combine with your simulated data (sim_dat0 from line 59)

# test_f <- sp_files[grepl("mpas-at-5-year-intervals", sp_files)]
# sim_dat0 <- readRDS(file.path(sample_dir, test_f[1]))
# # sim_dat0 <- bind_rows(
# #   readRDS(file.path(sample_dir, "yelloweye-rockfish/HBLL-OUT-N_mpa1.011_no_AR1_twenty-five_years_mpas-at-5-year-intervals.rds")),
# #   readRDS(file.path(sample_dir, "yelloweye-rockfish/HBLL-OUT-S_mpa1.011_no_AR1_twenty-five_years_mpas-at-5-year-intervals.rds"))
# # )
# sim_dat <- combine_hist_sim_data(sim_dat0, hist_data, 2047) |>
#   filter(replicate %in% 0:1)

# ggplot(data = sim_dat |> filter(year >= last_sampled_year)) +
#   geom_point(aes(x = X, y = Y, colour = factor(restricted), shape = factor(historical))) +
#   # geom_text(data = tibble(year = EVALUATION_YEARS +, X = 500, Y = 5900),
#   #   aes(x = X, y = Y, label = year), size = 5) +
#   scale_shape_manual(values = c(19, 21)) +
#   facet_wrap(~ year)

# test <- fit_simulation(
#   dat = sim_dat,
#   formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
#   spatial = "on",
#   spatiotemporal = "iid",
#   cutoff = 20,
#   control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
#   silent = FALSE
#   )
# meep()
# sanity(test)
# test
# #
# # combine surveys
# f <- sample_summary |>
#   filter(species == "yelloweye rockfish",
#   mpa_trend == 1.011, ar1_scenario == "no_AR1", time_scenario == "twenty-five_years",
#   plan == "MPAs at 5 year intervals") |>
# pull(file)
# sampled_data <- purrr::map_dfr(f, \(x) readRDS(file.path(sample_dir, x)))
# hist_data <- purrr::map_dfr(c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"),
#   \(x) get_hist_data(species, x, hist_path, hist_cache))
# combined_data <- combine_hist_sim_data(sampled_data, hist_data, 2047) |>
#   filter(replicate %in% 0:1) |>
#   mutate(survey_abbrev = "HBLL")
# # hbll_grid_subregion_lu <- readRDS(file.path("data-generated", "spatial","hbll-grid-subregion-lu.rds")) |>
# #   select(survey_abbrev, block_id, subregion, subregion_name)
# # test <- left_join(combined_data, hbll_grid_subregion_lu, by = c("survey_abbrev", "block_id"))

# test <- fit_simulation(
#   dat = combined_data,
#   formula = catch_prop ~ 0 + fyear + restricted:year_covariate,
#   spatial = "on",
#   spatiotemporal = "iid",
#   cutoff = 20,
#   control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
#   silent = FALSE
#   )
# meep()
# sanity(test)
# test

source(here::here("R", "00-fit-power-analysis-functions.R"))

# =============================================================================
# Defensive check: fit one combined HBLL model
# =============================================================================
hist_cache <- new.env(parent = emptyenv())
ye_hist <- purrr::map_dfr(c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"),
  ~get_hist_data("yelloweye rockfish", .x, hist_path, hist_cache))

ye_files <- list.files(file.path(sample_dir, "yelloweye-rockfish"), pattern = "mpa1.011.*status-quo", full.names = TRUE)
ye_samp <- purrr::map_dfr(ye_files, readRDS) |> filter(replicate == 1)
hbll_last_sampled_year
ye_combined <- combine_hist_sim_data(ye_samp, ye_hist, 2034) |> mutate(survey_abbrev = "HBLL")

ye_fit <- fit_simulation(ye_combined, formula = FORMULA, cutoff = 20, silent = FALSE)
meep()
sanity(ye_fit)
tidy(ye_fit, conf.int = TRUE) |> filter(term == "restricted:year_covariate")
get_model_pars(ye_fit)

stopifnot("Model converged" = all(sanity(ye_fit)))
stopifnot("Has historical data" = any(ye_combined$historical))
stopifnot("Has MPA sites" = sum(ye_combined$restricted) > 0)
stopifnot("No fit errors" = is.na(ye_trend$error_msg))
message("✓ Model checks passed")

# =============================================================================
# Main Workflow
# =============================================================================

message("\n=== Power Analysis: Model Fitting ===")
tictoc::tic("Starting power analysis")
future::plan(future::sequential)

setup_parallel(USE_PARALLEL, N_WORKERS)

sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))
# N_REPLICATES <- max(sampling_summary$n_replicates)
# N_REPLICATES
N_REPLICATES <- 50
results_dir <- here::here("data-generated", "power-results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

task_grid <- create_task_grid(sampling_summary, sample_dir) |>
  add_combined_survey_tasks(sampling_summary, sample_dir) |>
  dplyr::filter(
    #species %in% c("yelloweye rockfish", "lingcod", "pacific halibut"),
    survey_abbrev == "HBLL" # this combines all three surveys
    ) |>
  dplyr::filter(ar1_scenario == "fitted_AR1") |>
  dplyr::arrange(species)

cat("Task grid:\n")
print(task_grid)

message("Parameter combinations: ", nrow(task_grid))
message("Replicates per combination: ", N_REPLICATES)
message("Total models to fit: ", nrow(task_grid) * N_REPLICATES)
message("Parallel workers: ", if (is.null(N_WORKERS)) "auto" else N_WORKERS)

message("\n=== Executing Parallel Fitting ===")
message("Evaluation years: ", paste(EVALUATION_YEARS, collapse = ", "))
summary_stats <- execute_parallel_fitting_flat(
  task_grid = task_grid,
  results_dir = results_dir,
  hist_path = hist_path,
  sample_dir = sample_dir,
  sampling_summary = sampling_summary,
  n_reps_to_fit = N_REPLICATES,
  evaluation_years = EVALUATION_YEARS,
  .formula = FORMULA
)
# meep()
message("\n=== Creating Summary Catalog ===")

catalog <- create_summary_catalog(results_dir)

message("\n=== Combining All Results ===")
all_results <- combine_all_results(results_dir)

message("\n=== Fitting Complete ===")
message("Combos processed: ", nrow(summary_stats))
message("Total new fits: ", sum(summary_stats$n_new))
message("Total errors: ", sum(summary_stats$n_errors))
message("\nResults saved to: ", results_dir)

future::plan(future::sequential)
tictoc::toc()

# test <- readRDS(file.path(results_dir, "all-fitted-results.rds"))
# test <- readRDS(file.path(results_dir, "yelloweye-rockfish", "HBLL_mpa1.011_no_AR1_twenty-five_years_mpas-at-5-year-intervals_results.rds"))
# glimpse(test)

# distinct(test, survey_abbrev)

