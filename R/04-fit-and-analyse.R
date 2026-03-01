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
N_REPLICATES <- 100

if (Sys.info()['user'] %in% c("dunic", "anderson")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- 40
  N_REPLICATES <- 50
}

if (Sys.info()['user'] %in% c("jillian", "jilliandunic")) {
  USE_PARALLEL <- TRUE
  N_WORKERS <- ifelse(Sys.info()['user'] == "jillian", 10, 8)
  N_REPLICATES <- 100
}

FORMULA <- catch_prop ~ 0 + fyear + restricted:year_covariate
EVALUATION_YEARS <- c(2030, 2034, 2038, 2042, 2046)

# hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))

# # Testing
# sample_summary <- readRDS(file.path(sample_dir,  "sampling-summary.rds"))

# # species <- "lingcod"
# species <- "yelloweye rockfish"
# ar1_scenarios <- c("no_AR1", "moderate_AR1")
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

# =============================================================================
# Helper Functions
# =============================================================================

#' Create standardized error row
create_error_row <- function(combo, replicate, eval_year, error_message) {
  tibble(
    species = combo$species,
    survey_abbrev = combo$survey_abbrev,
    sim_mpa_trend = combo$mpa_trend,
    sim_ar1_scenario = combo$ar1_scenario,
    sim_time_scenario = combo$time_scenario,
    sampling_plan = combo$plan,
    replicate = replicate,
    eval_year = eval_year,
    estimate = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    converged = FALSE,
    sanity = NA_character_,
    error_msg = error_message,
    fit_formula = NA_character_,
    fit_family = NA_character_,
    fit_spatial = NA_character_,
    fit_spatiotemporal = NA_character_
  )
}

#' Generate result filename matching sampling data convention
generate_result_filename <- function(species, survey_abbrev, mpa_trend,
                                    ar1_scenario, time_scenario, plan) {
  plan_slug <- gsub("[^a-zA-Z0-9]", "-", plan) |>
    gsub("-+", "-", x = _) |>
    tolower()

  paste0(
    survey_abbrev, "_",
    "mpa", round(mpa_trend, digits = 3), "_",
    ar1_scenario, "_",
    time_scenario, "_",
    plan_slug,
    "_results.rds"
  ) |>
    gsub(" ", "-", x = _)
}

#' Lazy load historical data with per-worker caching
get_hist_data <- function(species, survey_abbrev, hist_path, cache_env) {
  key <- paste(species, survey_abbrev, sep = "___")
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }

  sp_hyp <- sp_to_hyphens(species)
  survey_hyp <- sp_to_hyphens(survey_abbrev)
  fname <- file.path(hist_path,
                     paste0(sp_hyp, "-", survey_hyp, "-processed.rds"))
  if (!file.exists(fname)) {
    stop("Historical data not found: ", fname)
  }

  hist_data <- readRDS(fname)
  assign(key, hist_data, envir = cache_env)
  hist_data
}

#' Combine historical and simulated data
combine_hist_sim_data <- function(sim_data, hist_data, eval_year, last_sampled_year) {
  hbll_last_sampled_year <- readRDS(file.path("data-generated", "hbll-last-sampled-year.rds"))

  sim_data_prep <- sim_data |>
    left_join(hbll_last_sampled_year, by = "survey_abbrev") |>
    mutate(
      catch_count = observed,
      historical = FALSE
    )

  combined_data <- bind_rows(hist_data, sim_data_prep) |>
    select(replicate, survey_abbrev, X, Y, block_id, restricted, historical,
           year, year_covariate, last_sampled_year,
           catch_count, hook_count, offset,
           plan, sim_mpa_trend, sim_ar1_scenario, sim_time_scenario) |>
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

  if (grepl("HBLL", survey_type)) {
    weights <- dat$hook_count
    offset <- NULL
  } else {
    weights <- NULL
    offset <- dat$offset
  }

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

#' Process one parameter combo (all replicates)
fit_parameter_combo <- function(combo,
                               .formula = catch_prop ~ 0 + fyear + restricted + restricted:year_covariate,
                               sample_file, results_dir,
                               hist_path, n_replicates = 50,
                               evaluation_years = EVALUATION_YEARS,
                               hist_cache_env = new.env(parent = emptyenv())) {

  result_file <- generate_result_filename(
    combo$species, combo$survey_abbrev, combo$mpa_trend,
    combo$ar1_scenario, combo$time_scenario, combo$plan
  )
  result_path <- file.path(results_dir, sp_to_hyphens(combo$species), result_file)

  if (file.exists(result_path)) {
    existing_results <- readRDS(result_path)
    completed_combos <- existing_results |>
      filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
      distinct(replicate, eval_year) |>
      mutate(combo_id = paste(replicate, eval_year, sep = "_"))

    expected_combos <- expand.grid(
      replicate = 1:n_replicates,
      eval_year = evaluation_years
    ) |>
      mutate(combo_id = paste(replicate, eval_year, sep = "_"))

    combos_to_run <- expected_combos |>
      filter(!combo_id %in% completed_combos$combo_id)

    message("  Resume: ", nrow(completed_combos), " done, ",
            nrow(combos_to_run), " remaining")
  } else {
    existing_results <- NULL
    combos_to_run <- expand.grid(
      replicate = 1:n_replicates,
      eval_year = evaluation_years
    )
    message("  Starting: ", n_replicates, " replicates × ",
            length(evaluation_years), " eval years = ",
            nrow(combos_to_run), " fits")
  }

  if (nrow(combos_to_run) == 0) {
    message("  Complete: skipping")
    return(tibble(n_new = 0, n_total = nrow(existing_results), n_errors = 0))
  }

  sampled_data_all <- readRDS(sample_file)
  hist_data <- get_hist_data(combo$species, combo$survey_abbrev,
                             hist_path, hist_cache_env)

  combos_by_rep <- combos_to_run |>
    group_by(replicate) |>
    summarise(eval_years = list(sort(eval_year)), .groups = "drop")

  new_results <- purrr::map_dfr(seq_len(nrow(combos_by_rep)), function(i) {
    rep <- combos_by_rep$replicate[i]
    eval_years_to_fit <- combos_by_rep$eval_years[[i]]

    sampled_data_rep <- sampled_data_all |> filter(replicate == rep)

    if (nrow(sampled_data_rep) == 0) {
      return(purrr::map_dfr(eval_years_to_fit, function(ey) {
        create_error_row(combo, rep, ey, "Missing replicate in sampled data")
      }))
    }

    purrr::map_dfr(eval_years_to_fit, function(eval_year) {
      tryCatch({
        combined_data <- combine_hist_sim_data(sampled_data_rep, hist_data, eval_year)

        fit <- fit_simulation(
          dat = combined_data,
          formula = .formula,
          spatial = "on",
          spatiotemporal = "iid",
          family = betabinomial(link = "cloglog"),
          cutoff = 20,
          control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
          silent = FALSE
        )

        trend_results <- extract_trend_estimate(fit, "restricted:year_covariate")

        tibble(
          species = combo$species,
          survey_abbrev = combo$survey_abbrev,
          sim_mpa_trend = combo$mpa_trend,
          sim_ar1_scenario = combo$ar1_scenario,
          sim_time_scenario = combo$time_scenario,
          sampling_plan = combo$plan,
          replicate = rep,
          eval_year = eval_year,
          estimate = trend_results$estimate,
          se = trend_results$se,
          ci_lower = trend_results$ci_lower,
          ci_upper = trend_results$ci_upper,
          converged = trend_results$converged,
          sanity = trend_results$sanity,
          error_msg = trend_results$error_msg,
          fit_formula = if (!is.null(fit$error)) NA_character_ else deparse1(formula(fit)),
          fit_family = if (!is.null(fit$error)) NA_character_ else clean_family_name(fit),
          fit_spatial = if (!is.null(fit$error)) NA_character_ else fit$spatial,
          fit_spatiotemporal = if (!is.null(fit$error)) NA_character_ else fit$spatiotemporal
        )

      }, error = function(e) {
        create_error_row(combo, rep, eval_year, paste("Worker error:", e$message))
      })
    })
  })

  if (!is.null(existing_results)) {
    combined_results <- bind_rows(new_results, existing_results) |>
      distinct(replicate, eval_year, .keep_all = TRUE) |>
      arrange(replicate, eval_year)
  } else {
    combined_results <- new_results |>
      arrange(replicate, eval_year)
  }

  dir.create(dirname(result_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(combined_results, result_path)

  tibble(
    n_new = nrow(new_results),
    n_total = nrow(combined_results),
    n_errors = sum(!is.na(new_results$error_msg))
  )
}

#' Create combo-level task grid
create_task_grid <- function(sampling_summary, sample_dir) {
  task_grid <- sampling_summary |>
    distinct(species, survey_abbrev, mpa_trend, ar1_scenario,
             time_scenario, plan, file, n_replicates) |>
    mutate(
      sample_file = file.path(sample_dir, file)
    ) |>
    select(species, survey_abbrev, mpa_trend, ar1_scenario,
           time_scenario, plan, sample_file, n_replicates)

  return(task_grid)
}

#' Execute parallel fitting with progress reporting
execute_parallel_fitting <- function(task_grid, results_dir,
                                    hist_path, n_reps_to_fit = 50,
                                    evaluation_years = EVALUATION_YEARS,
                                    .formula = catch_prop ~ 0 + fyear + restricted:year_covariate) {

  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")

  summary_stats <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(task_grid))

    furrr::future_pmap_dfr(
      task_grid,
      function(species, survey_abbrev, mpa_trend, ar1_scenario,
               time_scenario, plan, sample_file, n_replicates, ...) {

        combo <- tibble(
          species = species,
          survey_abbrev = survey_abbrev,
          mpa_trend = mpa_trend,
          ar1_scenario = ar1_scenario,
          time_scenario = time_scenario,
          plan = plan
        )

        result_summary <- fit_parameter_combo(
          combo = combo,
          sample_file = sample_file,
          results_dir = results_dir,
          hist_path = hist_path,
          n_replicates = n_reps_to_fit,
          evaluation_years = evaluation_years,
          .formula = .formula
        )

        p()

        result_summary |>
          mutate(
            species = species,
            survey_abbrev = survey_abbrev,
            mpa_trend = mpa_trend,
            ar1_scenario = ar1_scenario,
            time_scenario = time_scenario,
            plan = plan
          )
      },
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = c("fit_parameter_combo", "get_hist_data",
                   "fit_simulation", "extract_trend_estimate",
                   "combine_hist_sim_data", "create_error_row",
                   "generate_result_filename", "sp_to_hyphens",
                   "clean_family_name", "summarise_sanity",
                   "results_dir", "hist_path", "n_reps_to_fit",
                   "evaluation_years", "p", ".formula"),
        packages = c("dplyr", "sdmTMB")
      )
    )
  })

  return(summary_stats)
}

#' Create summary catalog from results
create_summary_catalog <- function(results_dir) {
  result_files <- list.files(results_dir, pattern = "_results\\.rds$",
                             recursive = TRUE, full.names = TRUE)

  if (length(result_files) == 0) {
    message("No result files found")
    return(tibble())
  }

  catalog <- purrr::map_dfr(result_files, function(fpath) {
    results <- readRDS(fpath)

    summary_row <- results |>
      group_by(eval_year) |>
      summarise(
        n_replicates = n(),
        n_converged = sum(converged, na.rm = TRUE),
        n_errors = sum(!is.na(error_msg)),
        mean_estimate = mean(estimate, na.rm = TRUE),
        sd_estimate = sd(estimate, na.rm = TRUE),
        power = mean(ci_lower > 0, na.rm = TRUE),
        .groups = "drop"
      )

    param_row <- results |>
      slice(1) |>
      select(species, survey_abbrev, sim_mpa_trend,
             sim_ar1_scenario, sim_time_scenario, sampling_plan)

    tidyr::crossing(param_row, summary_row) |>
      mutate(file = basename(fpath))
  })

  saveRDS(catalog, file.path(results_dir, "power-results-summary.rds"))
  message("\nSummary catalog saved: power-results-summary.rds")

  return(catalog)
}

#' Combine all individual result files into one dataset
combine_all_results <- function(results_dir) {
  result_files <- list.files(results_dir, pattern = "_results\\.rds$",
                             recursive = TRUE, full.names = TRUE)

  if (length(result_files) == 0) {
    warning("No result files found to combine")
    return(tibble())
  }

  message("Combining ", length(result_files), " result files...")

  all_results <- purrr::map_dfr(result_files, function(fpath) {
    tryCatch({
      readRDS(fpath)
    }, error = function(e) {
      warning("Failed to load: ", fpath, " - ", e$message)
      tibble()
    })
  })

  message("Total rows: ", nrow(all_results))
  message("Unique species: ", length(unique(all_results$species)))
  message("Unique plans: ", length(unique(all_results$sampling_plan)))
  if (length(result_files) > 0 && nrow(all_results) > 0) {
    message("Replicates per combo: ", round(nrow(all_results) / length(result_files), 1))
  }

  output_file <- file.path(results_dir, "all-fitted-results.rds")
  saveRDS(all_results, output_file)
  message("Saved combined results: ", output_file)

  return(all_results)
}

# =============================================================================
# Main Workflow
# =============================================================================

message("\n=== Power Analysis: Model Fitting ===")
tictoc::tic("Starting power analysis")
future::plan(future::sequential)
setup_parallel(USE_PARALLEL, N_WORKERS)

sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))

task_grid <- create_task_grid(sampling_summary, sample_dir) |>
  filter(species == "yelloweye rockfish",
         survey_abbrev == "HBLL OUT N")

message("Parameter combinations: ", nrow(task_grid))
message("Replicates per combination: ", N_REPLICATES)
message("Total models to fit: ", nrow(task_grid) * N_REPLICATES)
message("Parallel workers: ", if (is.null(N_WORKERS)) "auto" else N_WORKERS)

message("\n=== Executing Parallel Fitting ===")
message("Evaluation years: ", paste(EVALUATION_YEARS, collapse = ", "))
summary_stats <- execute_parallel_fitting(
  task_grid = task_grid,
  results_dir = results_dir,
  hist_path = hist_path,
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
if (nrow(summary_stats) > 0 && sum(summary_stats$n_total) > 0) {
  message("Mean convergence rate: ",
          round(100 * sum(summary_stats$n_new - summary_stats$n_errors) / sum(summary_stats$n_new), 1), "%")
}

message("\n=== Summary ===")
if (nrow(summary_stats) > 0) {
  print(summary_stats |> select(species, plan, n_new, n_total, n_errors))
}
message("\nResults saved to: ", results_dir)

future::plan(future::sequential)
tictoc::toc()

# test <- readRDS(file.path(results_dir, "all-fitted-results.rds"))

# glimpse(test)
