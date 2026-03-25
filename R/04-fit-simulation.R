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

hist_path <- here::here("data-generated", "cleaned-species-data")
sample_dir <- here::here("data-generated", "sampled-data")
results_dir <- here::here("data-generated", "power-results")
# hist_path <- here::here("data-generated", "historical-data-processed")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- TRUE
N_WORKERS <- 8L #NULL
N_REPLICATES <- 25

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 78
  N_REPLICATES <- 10
}

if (Sys.info()['user'] %in% c("jillian", "jilliandunic")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- ifelse(Sys.info()['user'] == "jillian", 10, 8)
  N_REPLICATES <- 100
}

FORMULA <- catch_prop ~ 0 + fyear + restricted + year_covariate +
  # restricted:future_step +
  restricted:future_year_covariate

# FORMULA <- catch_prop ~ 0 + fyear + restricted + year_covariate +
#   restricted:future_step + restricted:future_year_covariate

TREND_PARAM <- "restricted:future_year_covariate"
EVALUATION_YEARS <- c(2030, 2034, 2038, 2042, 2046)
# EVALUATION_YEARS <- c(2030, 2038, 2046)

# EVALUATION_YEARS <- c(2046)
# EVALUATION_YEARS <- c(2038)

hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))

# =============================================================================
# Testing/Debugging options
# =============================================================================
RUN_DEFENSIVE_CHECKS <- FALSE  # Run defensive checks through execute_parallel_fitting
SAVE_TEST_FITS <- TRUE  # Save test fit objects for inspection
TEST_FITS_DIR <- here::here("data-generated", "test-fits")
dir.create(TEST_FITS_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
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
# Defensive checks: fit sample models through execute_parallel_fitting
# =============================================================================
if (RUN_DEFENSIVE_CHECKS) {
  message("\n=== Running Defensive Checks ===")

  # Load sampling summary early
  sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))

  # Create test task grid for specific species
  test_species <- c("yelloweye rockfish", "silvergray rockfish")

  # Filter sampling summary FIRST
  test_sampling_summary <- sampling_summary |>
    filter(
      species %in% test_species,
      ar1_scenario == "fitted_AR1",
      time_scenario == "twenty-five_years",
      plan == "status quo",
      mpa_trend %in% c(1.021, 1.018),
      replicate == 1  # Just test first replicate
    )

  # Then create task grid from filtered summary
  test_task_grid <- test_sampling_summary |>
    create_task_grid(sample_dir = sample_dir) |>
    add_combined_survey_tasks(test_sampling_summary, sample_dir) |>
    filter(survey_abbrev == "HBLL")

  message("Test task grid: ", nrow(test_task_grid), " combinations")
  print(test_task_grid)

  # Run through actual parallel fitting pathway
  test_results <- execute_parallel_fitting(
    task_grid = test_task_grid,
    results_dir = TEST_FITS_DIR,  # Save to test directory
    hist_path = hist_path,
    sample_dir = sample_dir,
    sampling_summary = test_sampling_summary,  # Use filtered version
    replicate_filter = 1,  # Just first replicate
    evaluation_years_filter = 2046,  # Just one year for testing
    save_fits = SAVE_TEST_FITS,  # Save fit objects
    fits_dir = TEST_FITS_DIR,  # Where to save fits
    evaluation_years = EVALUATION_YEARS,
    .formula = FORMULA,
    .trend_param = TREND_PARAM
  )

  message("\n=== Defensive Check Results ===")
  print(test_results)

  # Load and inspect saved results
  for (sp in test_species) {
    result_file <- list.files(
      file.path(TEST_FITS_DIR, sp_to_hyphens(sp)),
      pattern = "HBLL.*results\\.rds$",
      full.names = TRUE
    )

    if (length(result_file) > 0) {
      results <- readRDS(result_file[1])
      message("\n", sp, " results:")
      print(results)

      if (all(!is.na(results$sanity) & results$sanity == "ok")) {
        message("  ✓ Model converged")
      } else {
        message("  ✗ Model convergence issues")
      }
    }

    # Show saved fit files
    fit_files <- list.files(
      file.path(TEST_FITS_DIR, sp_to_hyphens(sp)),
      pattern = "_fit\\.rds$",
      full.names = TRUE
    )
    if (length(fit_files) > 0) {
      message("  Fit objects saved: ", length(fit_files))
      message("  ", paste(basename(fit_files), collapse = "\n    "))
    }
  }

  message("\n✓ Defensive checks complete")
  message("Results saved to: ", TEST_FITS_DIR)

  # Stop here if just testing
  stop("Defensive checks complete. Set RUN_DEFENSIVE_CHECKS <- FALSE to run full analysis.")
}
meep()

# test_files <- list.files(
#   file.path(TEST_FITS_DIR, sp_to_hyphens("silvergray rockfish")),
#   pattern = "HBLL.*rep.*\\.rds$",
#   full.names = TRUE
# )

# test_fits <- readRDS(test_files)

# test_fits$fit

# test_fits$fit$data |>
#   ggplot() +
#   aes(x = X, y = Y, colour = factor(restricted), shape = factor(historical)) +
#   geom_point() +
#   facet_wrap(~ year)

# =============================================================================
# Main Workflow
# =============================================================================

message("\n=== Power Analysis: Model Fitting ===")
tictoc::tic("Starting power analysis")

setup_parallel(USE_PARALLEL, N_WORKERS)

sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))


### SETTINGS
# =============================================================================
# Task grid filtering (set to NULL to use all available)
# =============================================================================
# FILTER_SPECIES <- "pacific halibut" #"silvergray rockfish"         # NULL = all species
FILTER_SPECIES <- c(
  "yelloweye rockfish",
   "north pacific spiny dogfish",
   "lingcod",
  "quillback rockfish",
  # "pacific halibut",
   "canary rockfish",
   "silvergray rockfish"
)
FILTER_SURVEY <- NULL          # NULL = all surveys
FILTER_MPA_TREND <- NULL #1.018       # NULL = all MPA trends
# FILTER_MPA_TREND <- 1.088 #1.018       # NULL = all MPA trends
# FILTER_MPA_TREND <- 1.009 #1.018       # NULL = all MPA trends
FILTER_AR1_SCENARIO <- "fitted_AR1"    # NULL = all AR1 scenarios
# FILTER_TIME_SCENARIO <- "twenty-five_years"   # NULL = all time scenarios
FILTER_TIME_SCENARIO <- NULL   # NULL = all time scenarios
# FILTER_PLAN <- "status quo" #c("status quo", "MPAs every 4 years")  # NULL = all plans
FILTER_PLAN <- c("historical survey-year bootstrap") #c("status quo", "MPAs every 4 years")  # NULL = all plans
FILTER_TIME_SCENARIO <- "twenty-five_years"   # NULL = all time scenarios
FILTER_REPLICATES <- 1:20    # NULL = all available replicates, e.g., 1:50
FILTER_EVALUATION_YEARS <- NULL #c(2046)  # NULL = all evaluation years

results_dir <- here::here("data-generated", "power-results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Apply task grid filters
message("\n=== Applying Task Grid Filters ===")
sampling_summary_filtered <- apply_filters_to_sampling_summary(
  sampling_summary,
  species_filter = FILTER_SPECIES,
  survey_filter = FILTER_SURVEY,
  mpa_trend_filter = FILTER_MPA_TREND,
  ar1_filter = FILTER_AR1_SCENARIO,
  time_filter = FILTER_TIME_SCENARIO,
  plan_filter = FILTER_PLAN
)

if (!is.null(FILTER_REPLICATES)) {
  sampling_summary_filtered <- sampling_summary_filtered |>
    filter(replicate %in% FILTER_REPLICATES)
  message("Filtering to replicates: ", paste(FILTER_REPLICATES, collapse = ", "))
}

message("Rows after filtering: ", nrow(sampling_summary_filtered), " (from ", nrow(sampling_summary), ")")

# Create task grid from filtered summary
task_grid <- create_task_grid(sampling_summary_filtered, sample_dir) |>
  add_combined_survey_tasks(sampling_summary_filtered, sample_dir) |>
  # Additional hardcoded filters if needed
  dplyr::filter(
    survey_abbrev == "HBLL" # this combines all three surveys
  ) |>
  dplyr::filter(ar1_scenario == "fitted_AR1") |>
  dplyr::arrange(species)

cat("\nTask grid:\n")
print(task_grid)

# Calculate expected number of fits
eval_years_to_use <- if (!is.null(FILTER_EVALUATION_YEARS)) {
  intersect(EVALUATION_YEARS, FILTER_EVALUATION_YEARS)
} else {
  EVALUATION_YEARS
}

# Determine replicates that will be fit
n_reps_msg <- if (!is.null(FILTER_REPLICATES)) {
  paste0(length(FILTER_REPLICATES), " (", paste(FILTER_REPLICATES, collapse = ", "), ")")
} else {
  paste0(max(task_grid$n_replicates), " (all available)")
}

message("\n=== Fitting Configuration ===")
message("Parameter combinations: ", nrow(task_grid))
message("Replicates per combination: ", n_reps_msg)
message("Evaluation years: ", paste(eval_years_to_use, collapse = ", "))
message("Parallel workers: ", if (is.null(N_WORKERS)) "auto" else N_WORKERS)

message("\n=== Executing Parallel Fitting ===")
summary_stats <- execute_parallel_fitting(
  task_grid = task_grid,
  results_dir = results_dir,
  hist_path = hist_path,
  sample_dir = sample_dir,
  sampling_summary = sampling_summary_filtered,
  replicate_filter = FILTER_REPLICATES,
  evaluation_years_filter = FILTER_EVALUATION_YEARS,
  evaluation_years = EVALUATION_YEARS,
  .formula = FORMULA,
  .trend_param = TREND_PARAM
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
meep()
# test <- readRDS(file.path(results_dir, "all-fitted-results.rds"))
# test <- readRDS(file.path(results_dir, "yelloweye-rockfish", "HBLL_mpa1.011_no_AR1_twenty-five_years_mpas-at-5-year-intervals_results.rds"))
# glimpse(test)

# distinct(test, survey_abbrev)
