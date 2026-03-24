source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

library(dplyr)
library(purrr)

hist_path <- here::here("data-generated", "cleaned-species-data")
sample_dir <- here::here("data-generated", "sampled-data")
sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))

USE_PARALLEL <- TRUE
N_WORKERS <- 5L

map_fn <- setup_parallel(USE_PARALLEL, N_WORKERS)

# Match one scenario from 04-fit-simulation.R as closely as possible.
species <- "yelloweye rockfish"
survey_abbrev <- "HBLL" # use "HBLL" for the combined-survey fit
exp(log(c(1.05, 1.1, 1.25, 1.5)) / 25)
log(c(1.05, 1.1, 1.25, 1.5)) / 25
mpa_trend <- 1.016
ar1_scenario <- "fitted_AR1"
time_scenario <- "twenty-five_years"
plan <- "status quo"
replicate <- 1
eval_years <- c(2030, 2034, 2038, 2042, 2046)
# SIMULATED_LOCATION_MODE <- "all_sampled"
# SIMULATED_LOCATION_MODE <- "historical_locations_only"
# MODEL_MODE <- "baseline"
MODEL_MODE <- "future_step"

if (MODEL_MODE == "baseline") {
  FORMULA <- catch_prop ~ 0 + fyear + restricted*year_covariate
  TARGET_TERMS <- "restricted:year_covariate"
} else if (MODEL_MODE == "future_step") {
  FORMULA <- catch_prop ~ 0 + fyear + restricted + year_covariate +
    restricted:future_step +
    restricted:future_year_covariate
  TARGET_TERMS <- c(
    "restricted:future_step",
    "restricted:future_year_covariate"
  )
} else {
  stop("Unknown MODEL_MODE: ", MODEL_MODE)
}

component_surveys <- if (survey_abbrev == "HBLL") {
  c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
} else {
  NULL
}

hist_cache <- new.env(parent = emptyenv())

filter_simulated_locations <- function(dat, mode = c("all_sampled", "historical_locations_only")) {
  mode <- match.arg(mode)

  if (mode == "all_sampled") {
    return(dat)
  }

  if (!"historical_location" %in% names(dat)) {
    stop("`historical_location` column is missing from sampled data.")
  }

  dat |>
    filter(historical_location == 1)
}

prepare_fit_data <- function(dat, model_mode = c("baseline", "future_step")) {
  model_mode <- match.arg(model_mode)

  dat <- dat |>
    mutate(
      fyear_value = if_else(historical, as.character(year), "future"),
      fyear = factor(fyear_value)
    )

  if (model_mode == "future_step") {
    dat |>
      mutate(
        future_step = as.integer(historical == FALSE),
        future_year_covariate = if_else(historical, 0, year_covariate)
      )
  } else {
    dat
  }
}

sampled_data <- load_sampled_data(
  species = species,
  survey_abbrev = if (is.null(component_surveys)) survey_abbrev else component_surveys,
  mpa_trend = mpa_trend,
  ar1_scenario = ar1_scenario,
  time_scenario = time_scenario,
  plan = plan,
  replicates = replicate,
  sampling_summary = sampling_summary,
  sample_dir = sample_dir
)

sampled_data <- filter_simulated_locations(sampled_data, SIMULATED_LOCATION_MODE)

hist_data <- if (!is.null(component_surveys)) {
  purrr::map_dfr(component_surveys, ~get_hist_data(species, .x, hist_path, hist_cache))
} else {
  get_hist_data(species, survey_abbrev, hist_path, hist_cache)
}

fit_one_year <- function(eval_year) {
  combined_data <- combine_hist_sim_data(sampled_data, hist_data, eval_year)
  combined_data <- prepare_fit_data(combined_data, MODEL_MODE)

  if (!is.null(component_surveys)) {
    combined_data <- combined_data |>
      mutate(survey_abbrev = survey_abbrev)
  }

  fit <- fit_simulation(
    dat = combined_data,
    formula = FORMULA,
    spatial = "on",
    spatiotemporal = "iid",
    family = betabinomial(link = "cloglog"),
    cutoff = 20,
    control = sdmTMBcontrol(
      collapse_spatial_variance = TRUE,
      multiphase = FALSE,
      profile = FALSE,
      newton_loops = 0L
    ),
    silent = TRUE
    # time_varying = ~1,
    # extra_time = seq(min(combined_data$year), max(combined_data$year)),
    # time_varying_type = "ar1"
  )

  list(
    eval_year = eval_year,
    data = combined_data,
    fit = fit
  )
}

fits_by_year <- if (USE_PARALLEL) {
  map_fn(eval_years, fit_one_year, .options = furrr::furrr_options(seed = TRUE))
} else {
  map_fn(eval_years, fit_one_year)
}

names(fits_by_year) <- as.character(eval_years)

down_curve <- purrr::map_dfr(fits_by_year, function(x) {
  fit <- x$fit

  if (!is.null(fit$error) && isTRUE(fit$error)) {
    return(tibble(
      eval_year = x$eval_year,
      term = TARGET_TERMS,
      estimate = NA_real_,
      std.error = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      significant = NA,
      error_msg = fit$message
    ))
  }

  tidy(fit, conf.int = TRUE) |>
    filter(term %in% TARGET_TERMS) |>
    transmute(
      eval_year = x$eval_year,
      term,
      estimate,
      std.error,
      conf.low,
      conf.high,
      significant = conf.low > 0 | conf.high < 0,
      error_msg = NA_character_
    )
})

down_curve |> filter(term == "restricted:future_year_covariate")

combined_plot_data <- combine_hist_sim_data(sampled_data, hist_data, max(eval_years))
combined_plot_data <- prepare_fit_data(combined_plot_data, MODEL_MODE)

combined_plot_data |>
  group_by(year, restricted) |>
  summarise(
    mean_cp = mean(catch_prop, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  group_by(restricted) |>
  # mutate(mean_cp = mean_cp - mean(mean_cp)) |>
  # filter(year > 2023) |>
  ggplot(aes(year, mean_cp, colour = factor(restricted))) + geom_path() +
  geom_point() +
  geom_vline(xintercept = hist_data$year |> max()) +
  geom_vline(xintercept = eval_years, lty = 2) +
  ggtitle(paste(SIMULATED_LOCATION_MODE, MODEL_MODE, sep = " | ")) +
  scale_y_log10()

# combined_data |>
#   group_by(year, restricted) |>
#   summarise(mean_cp = mean(catch_prop, na.rm = TRUE), .groups = "drop") |>
#   tidyr::pivot_wider(names_from = restricted, values_from = mean_cp) |>
#   mutate(diff = `1` - `0`) |>
#   ggplot(aes(year, diff)) + geom_path()
