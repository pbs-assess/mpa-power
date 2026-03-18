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
    fit_spatiotemporal = NA_character_,
    sigma_O = NA_real_,
    sigma_E = NA_real_,
    range = NA_real_
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
combine_hist_sim_data <- function(sim_data, hist_data, eval_year) {
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

  if (any(grepl("HBLL", survey_type))) {
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

# TODO: move this into extract_trend estimate?
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

#' Process one parameter combo (all replicates)
fit_parameter_combo <- function(combo,
                               .formula = catch_prop ~ 0 + fyear + restricted + restricted:year_covariate,
                               sample_file, component_surveys, sample_dir, results_dir,
                               hist_path, n_replicates = 50,
                               evaluation_years = EVALUATION_YEARS,
                               sampling_summary = NULL,
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

  if (!is.null(component_surveys)) {
    sample_files <- sampling_summary |>
      filter(
        species == combo$species,
        survey_abbrev %in% component_surveys,
        mpa_trend == combo$mpa_trend,
        ar1_scenario == combo$ar1_scenario,
        time_scenario == combo$time_scenario,
        plan == combo$plan
      ) |>
      pull(file)

    sampled_data_all <- purrr::map_dfr(sample_files, ~readRDS(file.path(sample_dir, .x)))
    hist_data <- purrr::map_dfr(component_surveys, ~get_hist_data(combo$species, .x, hist_path, hist_cache_env))
  } else {
    sampled_data_all <- readRDS(sample_file)
    hist_data <- get_hist_data(combo$species, combo$survey_abbrev,
                               hist_path, hist_cache_env)
  }

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

        if (!is.null(component_surveys)) {
          combined_data <- combined_data |> mutate(survey_abbrev = combo$survey_abbrev)
        }

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
        re_pars <- extract_re_pars(fit)

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
          # converged = trend_results$sanity=="ok",
          sanity = trend_results$sanity,
          error_msg = trend_results$error_msg,
          fit_formula = if (!is.null(fit$error)) NA_character_ else deparse1(formula(fit)),
          fit_family = if (!is.null(fit$error)) NA_character_ else clean_family_name(fit),
          fit_spatial = if (!is.null(fit$error)) NA_character_ else fit$spatial,
          fit_spatiotemporal = if (!is.null(fit$error)) NA_character_ else fit$spatiotemporal,
          sigma_O = re_pars$sigma_O,
          sigma_E = re_pars$sigma_E,
          range = re_pars$range
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
      sample_file = file.path(sample_dir, file),
      component_surveys = list(NULL)
    ) |>
    select(species, survey_abbrev, mpa_trend, ar1_scenario,
           time_scenario, plan, sample_file, component_surveys, n_replicates)

  return(task_grid)
}

#' Add combined survey tasks to task grid
add_combined_survey_tasks <- function(task_grid, sampling_summary, sample_dir) {
  survey_combinations <- list(
    "HBLL" = c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
  )

  combined_tasks <- purrr::map_dfr(names(survey_combinations), function(combined_name) {
    component_surveys <- survey_combinations[[combined_name]]

    sampling_summary |>
      filter(survey_abbrev %in% component_surveys) |>
      distinct(species, mpa_trend, ar1_scenario, time_scenario, plan, n_replicates) |>
      mutate(
        survey_abbrev = combined_name,
        component_surveys = list(component_surveys),
        sample_file = NA_character_
      )
  })

  bind_rows(task_grid, combined_tasks)
}

#' Execute parallel fitting with progress reporting
execute_parallel_fitting <- function(task_grid, results_dir,
                                    hist_path, sample_dir, sampling_summary,
                                    n_reps_to_fit = 50,
                                    evaluation_years = EVALUATION_YEARS,
                                    .formula = catch_prop ~ 0 + fyear + restricted:year_covariate) {

  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")

  summary_stats <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(task_grid))

    furrr::future_pmap_dfr(
      task_grid,
      function(species, survey_abbrev, mpa_trend, ar1_scenario,
               time_scenario, plan, sample_file, component_surveys, n_replicates, ...) {

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
          component_surveys = component_surveys,
          sample_dir = sample_dir,
          results_dir = results_dir,
          hist_path = hist_path,
          n_replicates = n_reps_to_fit,
          evaluation_years = evaluation_years,
          sampling_summary = sampling_summary,
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
                   "fit_simulation", "extract_trend_estimate", "extract_re_pars",
                   "combine_hist_sim_data", "create_error_row",
                   "generate_result_filename", "sp_to_hyphens",
                   "clean_family_name", "summarise_sanity",
                   "results_dir", "hist_path", "sample_dir", "sampling_summary",
                   "n_reps_to_fit", "evaluation_years", "p", ".formula"),
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
        # n_converged = sum(sanity=="ok", na.rm = TRUE),
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