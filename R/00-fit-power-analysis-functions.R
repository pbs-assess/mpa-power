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

#' Apply filters to sampling summary
#'
#' Filter sampling summary tibble by various parameters. Used to create
#' subsets of the full analysis for testing or focused execution.
#'
#' @param sampling_summary Sampling summary tibble
#' @param species_filter Character vector of species names or NULL
#' @param survey_filter Character vector of survey abbreviations or NULL
#' @param mpa_trend_filter Numeric vector of MPA trend values or NULL
#' @param ar1_filter Character vector of AR1 scenario names or NULL
#' @param time_filter Character vector of time scenario names or NULL
#' @param plan_filter Character vector of plan names or NULL
#'
#' @return Filtered sampling summary tibble
apply_filters_to_sampling_summary <- function(sampling_summary,
                                             species_filter = NULL,
                                             survey_filter = NULL,
                                             mpa_trend_filter = NULL,
                                             ar1_filter = NULL,
                                             time_filter = NULL,
                                             plan_filter = NULL) {
  result <- sampling_summary

  if (!is.null(species_filter)) {
    result <- result |> filter(species %in% species_filter)
    message("Filtering to species: ", paste(species_filter, collapse = ", "))
  }

  if (!is.null(survey_filter)) {
    result <- result |> filter(survey_abbrev %in% survey_filter)
    message("Filtering to surveys: ", paste(survey_filter, collapse = ", "))
  }

  if (!is.null(mpa_trend_filter)) {
    result <- result |> filter(mpa_trend %in% mpa_trend_filter)
    message("Filtering to MPA trends: ", paste(mpa_trend_filter, collapse = ", "))
  }

  if (!is.null(ar1_filter)) {
    result <- result |> filter(ar1_scenario %in% ar1_filter)
    message("Filtering to AR1 scenarios: ", paste(ar1_filter, collapse = ", "))
  }

  if (!is.null(time_filter)) {
    result <- result |> filter(time_scenario %in% time_filter)
    message("Filtering to time scenarios: ", paste(time_filter, collapse = ", "))
  }

  if (!is.null(plan_filter)) {
    result <- result |> filter(plan %in% plan_filter)
    message("Filtering to plans: ", paste(plan_filter, collapse = ", "))
  }

  message("Filtered sampling summary: ", nrow(result), " files (from ", nrow(sampling_summary), " total)")

  return(result)
}

#' Generate result filename matching sampling data convention
generate_result_filename <- function(species, survey_abbrev, mpa_trend,
                                    ar1_scenario, time_scenario, plan, replicate) {
  plan_slug <- gsub("[^a-zA-Z0-9]", "-", plan) |>
    gsub("-+", "-", x = _) |>
    tolower()

  paste0(
    survey_abbrev, "_",
    "mpa", round(mpa_trend, digits = 3), "_",
    ar1_scenario, "_",
    time_scenario, "_",
    plan_slug, "_",
    "rep", sprintf("%03d", replicate),
    "_results.rds"
  ) |>
    gsub(" ", "-", x = _)
}

#' Generate per-fit diagnostic log filename
generate_fit_log_filename <- function(species, survey_abbrev, mpa_trend,
                                      ar1_scenario, time_scenario, plan,
                                      replicate, eval_year) {
  generate_result_filename(
    species = species,
    survey_abbrev = survey_abbrev,
    mpa_trend = mpa_trend,
    ar1_scenario = ar1_scenario,
    time_scenario = time_scenario,
    plan = plan,
    replicate = replicate
  ) |>
    sub(
      "_results\\.rds$",
      paste0("_eval", eval_year, "_fit.log"),
      x = _
    )
}

#' Run a fit while capturing warnings/messages/console output to a per-fit log
run_fit_with_logging <- function(expr, combo, rep, eval_year, log_dir = NULL) {
  if (is.null(log_dir)) {
    return(force(expr))
  }

  sp_dir <- file.path(log_dir, sp_to_hyphens(combo$species))
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

  log_path <- file.path(
    sp_dir,
    generate_fit_log_filename(
      species = combo$species,
      survey_abbrev = combo$survey_abbrev,
      mpa_trend = combo$mpa_trend,
      ar1_scenario = combo$ar1_scenario,
      time_scenario = combo$time_scenario,
      plan = combo$plan,
      replicate = rep,
      eval_year = eval_year
    )
  )

  if (file.exists(log_path)) {
    unlink(log_path)
  }

  output_path <- tempfile("fit-output-", fileext = ".log")
  message_path <- tempfile("fit-message-", fileext = ".log")
  output_con <- file(output_path, open = "wt")
  message_con <- file(message_path, open = "wt")
  warnings_seen <- character()
  messages_seen <- character()
  error_seen <- NULL
  start_time <- Sys.time()
  base_output_sinks <- sink.number(type = "output")
  base_message_sinks <- sink.number(type = "message")

  cleanup_capture <- function() {
    while (sink.number(type = "message") > base_message_sinks) {
      sink(type = "message")
    }
    while (sink.number(type = "output") > base_output_sinks) {
      sink(type = "output")
    }
    try(close(output_con), silent = TRUE)
    try(close(message_con), silent = TRUE)
  }

  on.exit(cleanup_capture(), add = TRUE)

  sink(output_con, type = "output")
  sink(message_con, type = "message")

  result <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages_seen <<- c(messages_seen, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      error_seen <<- conditionMessage(e)
      NULL
    }
  )

  cleanup_capture()

  output_lines <- if (file.exists(output_path)) readLines(output_path, warn = FALSE) else character()
  message_lines <- if (file.exists(message_path)) readLines(message_path, warn = FALSE) else character()
  unlink(c(output_path, message_path))

  metadata_lines <- c(
    paste0("timestamp: ", format(start_time, tz = Sys.timezone(), usetz = TRUE)),
    paste0("species: ", combo$species),
    paste0("survey_abbrev: ", combo$survey_abbrev),
    paste0("mpa_trend: ", combo$mpa_trend),
    paste0("ar1_scenario: ", combo$ar1_scenario),
    paste0("time_scenario: ", combo$time_scenario),
    paste0("plan: ", combo$plan),
    paste0("replicate: ", rep),
    paste0("eval_year: ", eval_year)
  )

  sections <- c(metadata_lines)

  if (length(unique(messages_seen)) > 0) {
    sections <- c(
      sections,
      "",
      "[captured_messages]",
      unique(messages_seen)
    )
  }

  if (length(unique(warnings_seen)) > 0) {
    sections <- c(
      sections,
      "",
      "[captured_warnings]",
      unique(warnings_seen)
    )
  }

  if (length(message_lines) > 0) {
    sections <- c(
      sections,
      "",
      "[message_stream]",
      message_lines
    )
  }

  if (length(output_lines) > 0) {
    sections <- c(
      sections,
      "",
      "[output_stream]",
      output_lines
    )
  }

  if (!is.null(error_seen)) {
    sections <- c(
      sections,
      "",
      "[error]",
      error_seen
    )
  }

  if (length(sections) > length(metadata_lines)) {
    writeLines(sections, log_path)
  }

  if (!is.null(error_seen)) {
    stop(error_seen)
  }

  result
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
                     paste0(sp_hyp, "-", survey_hyp, ".rds"))
  if (!file.exists(fname)) {
    stop("Historical data not found: ", fname)
  }

  hist_data <- readRDS(fname)
  assign(key, hist_data, envir = cache_env)
  hist_data
}

#' Prepare combined historical/simulated data for power fitting
prepare_power_fit_data <- function(dat) {
  dat |>
    mutate(
      fyear_value = if_else(historical, as.character(year), "future"),
      fyear = factor(fyear_value),
      future_step = as.integer(historical == FALSE),
      future_year_covariate = if_else(historical, 0, year_covariate)
    )
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

  # Get unique surveys present in simulated data
  sim_surveys <- unique(sim_data_prep$survey_abbrev)

  # Filter historical data to only include those surveys
  hist_data_filtered <- hist_data |>
    filter(survey_abbrev %in% sim_surveys)

  combined_data <- bind_rows(hist_data_filtered, sim_data_prep) |>
    select(replicate, survey_abbrev, X, Y, block_id, restricted, historical,
           year, year_covariate, last_sampled_year,
           catch_count, hook_count, offset,
           plan, sim_mpa_trend, sim_ar1_scenario, sim_time_scenario) |>
    mutate(
      replicate = ifelse(historical, 0, replicate),
      catch_prop = catch_count / hook_count
    ) |>
    prepare_power_fit_data() |>
    filter(year <= eval_year)

  if (nrow(combined_data) == 0) {
    stop("No data remaining after filtering to eval_year = ", eval_year)
  }

  return(combined_data)
}

#' Fit sdmTMB model to sampled data
fit_simulation <- function(dat,
                           formula = catch_prop ~ 0 + fyear + restricted + year_covariate +
                             restricted:future_step + restricted:future_year_covariate,
                           spatial = "on",
                           spatiotemporal = "iid",
                           family = betabinomial(link = "cloglog"),
                           cutoff = 10,
                           control = sdmTMBcontrol(collapse_spatial_variance = TRUE,
                                                             multiphase = TRUE, profile = FALSE, newton_loops = 1L),
                           silent = TRUE, ...) {

  survey_type <- unique(dat$survey_abbrev)

  if (any(grepl("HBLL", survey_type))) {
    weights <- dat$hook_count
    offset <- NULL
  } else {
    weights <- NULL
    offset <- dat$offset
  }

  # f <- here::here("data-generated", "mesh-cache", "coastwide-mesh.rds")
  # if (!file.exists(f)) {
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  #   saveRDS(mesh, file = f)
  # } else {
  #   old_mesh <- readRDS(f)
  #   mesh <- make_mesh(dat, xy_cols = c("X", "Y"), mesh = old_mesh$mesh)
  # }

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
    control = control,
    ...
  )

  fit <-
    tryCatch({
      do.call(sdmTMB, params)
    }, error = function(e) {
      list(error = TRUE, message = e$message)
    })

  # d <- fit$data
  # ggplot(d, aes(year, catch_count / hook_count, colour = factor(restricted))) +
  #   geom_point()
  #
  # x <- d |> group_by(year, restricted) |>
  #   summarise(m = mean(catch_count/hook_count))
  # ggplot(x, aes(year, m, colour = factor(restricted))) +
  #   geom_point() +
  #   geom_vline(xintercept = 2024)
  #
  # ggplot(d, aes(X, Y, colour = factor(restricted), size = catch_count/hook_count)) +
  #   geom_point(pch = 21) + facet_wrap(~year)
  #
  # browser()
  fit
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
extract_trend_estimate <- function(fit, trend_param = "restricted:future_year_covariate") {
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

#' Process one parameter combo (specified replicates)
#'
#' Saves results as individual replicate files for atomic writes and parallel safety.
fit_single_replicate <- function(combo,
                                 rep,
                                 component_surveys,
                                 evaluation_years_filter = NULL,
                                 save_fits = FALSE,
                                 fits_dir = NULL,
                                 log_dir = NULL,
                                 .formula = catch_prop ~ 0 + fyear + restricted + year_covariate +
                                   restricted:future_step + restricted:future_year_covariate,
                                 .trend_param = "restricted:future_year_covariate",
                                 sample_dir, results_dir,
                                 hist_path,
                                 evaluation_years = EVALUATION_YEARS,
                                 sampling_summary = NULL,
                                 hist_cache_env = new.env(parent = emptyenv())) {

  evaluation_years_to_use <- if (!is.null(evaluation_years_filter)) {
    intersect(evaluation_years, evaluation_years_filter)
  } else {
    evaluation_years
  }

  if (length(evaluation_years_to_use) == 0) {
    message("  No evaluation years to fit after filtering")
    return(tibble(replicate = rep, n_new = 0, n_total = 0, n_errors = 0))
  }

  if (!is.null(component_surveys)) {
    hist_data <- purrr::map_dfr(component_surveys, ~get_hist_data(combo$species, .x, hist_path, hist_cache_env))
  } else {
    hist_data <- get_hist_data(combo$species, combo$survey_abbrev,
                               hist_path, hist_cache_env)
  }

  result_file <- generate_result_filename(
    combo$species, combo$survey_abbrev, combo$mpa_trend,
    combo$ar1_scenario, combo$time_scenario, combo$plan,
    replicate = rep
  )
  result_path <- file.path(results_dir, sp_to_hyphens(combo$species), result_file)

  if (file.exists(result_path)) {
    existing_results <- readRDS(result_path)
    completed_eval_years <- existing_results |>
      filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
      pull(eval_year)

    eval_years_to_fit <- setdiff(evaluation_years_to_use, completed_eval_years)

    if (length(eval_years_to_fit) == 0) {
      message("  Rep ", rep, ": complete, skipping")
      return(tibble(replicate = rep, n_new = 0, n_total = nrow(existing_results), n_errors = 0))
    }

    message("  Rep ", rep, ": ", length(completed_eval_years), " done, ",
            length(eval_years_to_fit), " remaining")
  } else {
    existing_results <- NULL
    eval_years_to_fit <- evaluation_years_to_use
    message("  Rep ", rep, ": fitting ", length(eval_years_to_fit), " eval years")
  }

  if (!is.null(component_surveys)) {
    sampled_data_rep <- load_sampled_data(
      species = combo$species,
      survey_abbrev = component_surveys,
      mpa_trend = combo$mpa_trend,
      ar1_scenario = combo$ar1_scenario,
      time_scenario = combo$time_scenario,
      plan = combo$plan,
      replicates = rep,
      sampling_summary = sampling_summary,
      sample_dir = sample_dir
    )
  } else {
    sampled_data_rep <- load_sampled_data(
      species = combo$species,
      survey_abbrev = combo$survey_abbrev,
      mpa_trend = combo$mpa_trend,
      ar1_scenario = combo$ar1_scenario,
      time_scenario = combo$time_scenario,
      plan = combo$plan,
      replicates = rep,
      sampling_summary = sampling_summary,
      sample_dir = sample_dir
    )
  }

  if (nrow(sampled_data_rep) == 0) {
    error_results <- purrr::map_dfr(eval_years_to_fit, function(ey) {
      create_error_row(combo, rep, ey, "Missing replicate in sampled data")
    })

    combined_results <- if (!is.null(existing_results)) {
      bind_rows(error_results, existing_results) |>
        distinct(eval_year, .keep_all = TRUE) |>
        arrange(eval_year)
    } else {
      error_results
    }

    dir.create(dirname(result_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(combined_results, result_path)

    return(tibble(replicate = rep, n_new = nrow(error_results),
                 n_total = nrow(combined_results), n_errors = nrow(error_results)))
  }

  new_results <- purrr::map_dfr(eval_years_to_fit, function(eval_year) {
    tryCatch({
      combined_data <- combine_hist_sim_data(sampled_data_rep, hist_data, eval_year)

      if (!is.null(component_surveys)) {
        combined_data <- combined_data |> mutate(survey_abbrev = combo$survey_abbrev)
      }

      fit <- run_fit_with_logging(
        expr = fit_simulation(
          dat = combined_data,
          formula = .formula,
          spatial = "on",
          spatiotemporal = "iid",
          family = betabinomial(link = "cloglog"),
          cutoff = 20,
          control = sdmTMBcontrol(collapse_spatial_variance = TRUE, multiphase = FALSE,
                                  profile = FALSE, newton_loops = 1L),
          silent = TRUE
        ),
        combo = combo,
        rep = rep,
        eval_year = eval_year,
        log_dir = log_dir
      )

      if (save_fits && !is.null(fits_dir)) {
        dir.create(file.path(fits_dir, sp_to_hyphens(combo$species)),
                  showWarnings = FALSE, recursive = TRUE)
        fit_file <- file.path(fits_dir, sp_to_hyphens(combo$species), paste0(
          sp_to_hyphens(combo$species), "_",
          combo$survey_abbrev, "_",
          "rep", sprintf("%03d", rep), "_",
          "eval", eval_year, "_fit.rds"
        ))
        saveRDS(list(
          fit = fit,
          data = combined_data,
          combo = combo,
          replicate = rep,
          eval_year = eval_year
        ), fit_file)
      }

      trend_results <- extract_trend_estimate(fit, .trend_param)
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

  combined_results <- if (!is.null(existing_results)) {
    bind_rows(new_results, existing_results) |>
      distinct(eval_year, .keep_all = TRUE) |>
      arrange(eval_year)
  } else {
    new_results |> arrange(eval_year)
  }

  dir.create(dirname(result_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(combined_results, result_path)

  tibble(
    replicate = rep,
    n_new = nrow(new_results),
    n_total = nrow(combined_results),
    n_errors = sum(!is.na(new_results$error_msg))
  )
}

fit_parameter_combo <- function(combo,
                               replicates_to_fit,
                               component_surveys,
                               evaluation_years_filter = NULL,
                               save_fits = FALSE,
                               fits_dir = NULL,
                               log_dir = NULL,
                               .formula = catch_prop ~ 0 + fyear + restricted + year_covariate +
                                 restricted:future_step + restricted:future_year_covariate,
                               .trend_param = "restricted:future_year_covariate",
                               sample_dir, results_dir,
                               hist_path,
                               evaluation_years = EVALUATION_YEARS,
                               sampling_summary = NULL,
                               hist_cache_env = new.env(parent = emptyenv())) {
  results_summary <- purrr::map_dfr(replicates_to_fit, function(rep) {
    fit_single_replicate(
      combo = combo,
      rep = rep,
      component_surveys = component_surveys,
      evaluation_years_filter = evaluation_years_filter,
      save_fits = save_fits,
      fits_dir = fits_dir,
      log_dir = log_dir,
      .formula = .formula,
      .trend_param = .trend_param,
      sample_dir = sample_dir,
      results_dir = results_dir,
      hist_path = hist_path,
      evaluation_years = evaluation_years,
      sampling_summary = sampling_summary,
      hist_cache_env = hist_cache_env
    )
  })

  # Return summary across all replicates
  tibble(
    n_new = sum(results_summary$n_new),
    n_total = sum(results_summary$n_total),
    n_errors = sum(results_summary$n_errors)
  )
}

#' Create combo-level task grid from sampling summary
#'
#' Groups the flattened sampling summary (with individual replicate files)
#' to create a task grid where each row represents a unique parameter
#' combination with a list of available replicates.
#'
#' @param sampling_summary Sampling summary tibble with replicate column
#' @param sample_dir Directory containing sampled data files
#'
#' @return Task grid tibble with replicates_available list column
create_task_grid <- function(sampling_summary, sample_dir) {
  task_grid <- sampling_summary |>
    # Group by combo (without replicate) to get list of available replicates
    group_by(species, survey_abbrev, mpa_trend, ar1_scenario,
             time_scenario, plan) |>
    summarise(
      replicates_available = list(sort(unique(replicate))),
      n_replicates = n(),
      .groups = "drop"
    ) |>
    mutate(
      sample_file = NA_character_,  # Not used with individual files
      component_surveys = list(NULL)
    ) |>
    select(species, survey_abbrev, mpa_trend, ar1_scenario,
           time_scenario, plan, replicates_available, n_replicates,
           sample_file, component_surveys)

  return(task_grid)
}

#' Expand combo-level task grid to one row per replicate
expand_task_grid_by_replicate <- function(task_grid, replicate_filter = NULL) {
  purrr::pmap_dfr(task_grid, function(species, survey_abbrev, mpa_trend,
                                      ar1_scenario, time_scenario, plan,
                                      replicates_available, n_replicates,
                                      sample_file, component_surveys, ...) {
    reps_to_fit <- if (!is.null(replicate_filter)) {
      intersect(replicates_available, replicate_filter)
    } else {
      replicates_available
    }

    if (length(reps_to_fit) == 0) {
      return(tibble())
    }

    tibble(
      species = species,
      survey_abbrev = survey_abbrev,
      mpa_trend = mpa_trend,
      ar1_scenario = ar1_scenario,
      time_scenario = time_scenario,
      plan = plan,
      replicate = reps_to_fit,
      component_surveys = list(component_surveys)
    )
  })
}

#' Add combined survey tasks to task grid
#'
#' Creates tasks for combined survey analyses (e.g., "HBLL" combining
#' OUT N, OUT S, and INS N). Groups by parameter combo to get list of
#' available replicates across all component surveys.
#'
#' @param task_grid Existing task grid from create_task_grid()
#' @param sampling_summary Sampling summary tibble
#' @param sample_dir Directory containing sampled data files
#'
#' @return Task grid with combined survey rows added
add_combined_survey_tasks <- function(task_grid, sampling_summary, sample_dir) {
  survey_combinations <- list(
    "HBLL" = c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
  )

  combined_tasks <- purrr::map_dfr(names(survey_combinations), function(combined_name) {
    component_surveys <- survey_combinations[[combined_name]]

    sampling_summary |>
      filter(survey_abbrev %in% component_surveys) |>
      # Group to get available replicates across component surveys
      group_by(species, mpa_trend, ar1_scenario, time_scenario, plan) |>
      summarise(
        # Only include replicates available in ALL component surveys
        replicates_available = list({
          reps_by_survey <- split(replicate, survey_abbrev)
          reps_common <- Reduce(intersect, reps_by_survey)
          sort(unique(reps_common))
        }),
        n_replicates = length(unlist(replicates_available)),
        .groups = "drop"
      ) |>
      mutate(
        survey_abbrev = combined_name,
        component_surveys = list(component_surveys),
        sample_file = NA_character_
      ) |>
      select(species, survey_abbrev, mpa_trend, ar1_scenario,
             time_scenario, plan, replicates_available, n_replicates,
             sample_file, component_surveys)
  })

  bind_rows(task_grid, combined_tasks)
}

#' Execute parallel fitting with progress reporting
execute_parallel_fitting <- function(task_grid, results_dir,
                                    hist_path, sample_dir, sampling_summary,
                                    replicate_filter = NULL,
                                    evaluation_years_filter = NULL,
                                    save_fits = FALSE,
                                    fits_dir = NULL,
                                    log_dir = NULL,
                                    evaluation_years = EVALUATION_YEARS,
                                    .formula = catch_prop ~ 0 + fyear + restricted + year_covariate +
                                      restricted:future_step + restricted:future_year_covariate,
                                    .trend_param = "restricted:future_year_covariate") {

  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")

  # Capture filters for parallel workers
  replicate_filter_captured <- replicate_filter
  evaluation_years_filter_captured <- evaluation_years_filter
  save_fits_captured <- save_fits
  fits_dir_captured <- fits_dir
  log_dir_captured <- log_dir
  replicate_task_grid <- expand_task_grid_by_replicate(task_grid, replicate_filter_captured)

  if (nrow(replicate_task_grid) == 0) {
    message("No replicate tasks to fit after filtering")
    return(tibble())
  }

  replicate_summary <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(replicate_task_grid))

    furrr::future_pmap_dfr(
      replicate_task_grid,
      function(species, survey_abbrev, mpa_trend, ar1_scenario,
               time_scenario, plan, replicate, component_surveys, ...) {

        combo <- tibble(
          species = species,
          survey_abbrev = survey_abbrev,
          mpa_trend = mpa_trend,
          ar1_scenario = ar1_scenario,
          time_scenario = time_scenario,
          plan = plan
        )

        result_summary <- fit_single_replicate(
          combo = combo,
          rep = replicate,
          component_surveys = component_surveys,
          evaluation_years_filter = evaluation_years_filter_captured,
          save_fits = save_fits_captured,
          fits_dir = fits_dir_captured,
          log_dir = log_dir_captured,
          sample_dir = sample_dir,
          results_dir = results_dir,
          hist_path = hist_path,
          evaluation_years = evaluation_years,
          sampling_summary = sampling_summary,
          .formula = .formula,
          .trend_param = .trend_param
        )

        p()

        result_summary |>
          mutate(
            species = species,
            survey_abbrev = survey_abbrev,
            mpa_trend = mpa_trend,
            ar1_scenario = ar1_scenario,
            time_scenario = time_scenario,
            plan = plan,
            replicate = replicate
          )
      },
      .options = furrr::furrr_options(
        seed = TRUE,
        scheduling = Inf,
        globals = c("fit_single_replicate", "get_hist_data",
                   "fit_simulation", "extract_trend_estimate", "extract_re_pars",
                   "combine_hist_sim_data", "create_error_row",
                   "generate_result_filename", "sp_to_hyphens",
                   "generate_fit_log_filename", "run_fit_with_logging",
                   "clean_family_name", "summarise_sanity", "load_sampled_data",
                   "results_dir", "hist_path", "sample_dir", "sampling_summary",
                   "evaluation_years_filter_captured", "save_fits_captured",
                   "fits_dir_captured", "log_dir_captured", "evaluation_years", "p", ".formula", ".trend_param",
                   "prepare_power_fit_data"),
        packages = c("dplyr", "sdmTMB")
      )
    )
  })

  summary_stats <- replicate_summary |>
    group_by(species, survey_abbrev, mpa_trend, ar1_scenario, time_scenario, plan) |>
    summarise(
      n_new = sum(n_new),
      n_total = sum(n_total),
      n_errors = sum(n_errors),
      .groups = "drop"
    )

  return(summary_stats)
}

#' Create summary catalog from results
create_summary_catalog <- function(results_dir) {
  result_files <- list.files(results_dir, pattern = "_rep[0-9]{3}_results\\.rds$",
                             recursive = TRUE, full.names = TRUE)

  if (length(result_files) == 0) {
    message("No result files found")
    return(tibble())
  }

  message("Reading ", length(result_files), " replicate result files...")

  # Load all results and combine
  all_results <- purrr::map_dfr(result_files, function(fpath) {
    tryCatch({
      readRDS(fpath)
    }, error = function(e) {
      warning("Failed to load: ", fpath, " - ", e$message)
      tibble()
    })
  })

  if (nrow(all_results) == 0) {
    message("No results to summarize")
    return(tibble())
  }

  # Group by parameter combination and eval_year to create summary
  catalog <- all_results |>
    group_by(species, survey_abbrev, sim_mpa_trend,
             sim_ar1_scenario, sim_time_scenario, sampling_plan, eval_year) |>
    summarise(
      n_replicates = n(),
      n_errors = sum(!is.na(error_msg)),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      power = mean(ci_lower > 0, na.rm = TRUE),
      .groups = "drop"
    )

  saveRDS(catalog, file.path(results_dir, "power-results-summary.rds"))
  message("\nSummary catalog saved: power-results-summary.rds")

  return(catalog)
}

#' Combine all individual result files into one dataset
combine_all_results <- function(results_dir) {
  result_files <- list.files(results_dir, pattern = "_rep[0-9]{3}_results\\.rds$",
                             recursive = TRUE, full.names = TRUE)

  if (length(result_files) == 0) {
    message("No result files found to combine")
    return(tibble())
  }

  message("Combining ", length(result_files), " replicate result files...")

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

  output_file <- file.path(results_dir, "all-fitted-results.rds")
  saveRDS(all_results, output_file)
  message("Saved combined results: ", output_file)

  return(all_results)
}
