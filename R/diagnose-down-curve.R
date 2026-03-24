source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

library(dplyr)
library(purrr)

hist_path <- here::here("data-generated", "cleaned-species-data")
sample_dir <- here::here("data-generated", "sampled-data")
sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))

# Match one scenario from 04-fit-simulation.R as closely as possible.
species <- "yelloweye rockfish"
survey_abbrev <- "HBLL INS N" # use "HBLL" for the combined-survey fit
mpa_trend <- 1.004
ar1_scenario <- "fitted_AR1"
time_scenario <- "twenty-five_years"
plan <- "status quo"
replicate <- 4
eval_years <- c(2030, 2034, 2038, 2042, 2046)

FORMULA <- catch_prop ~ 0 + fyear + restricted*year_covariate

component_surveys <- if (survey_abbrev == "HBLL") {
  c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
} else {
  NULL
}

hist_cache <- new.env(parent = emptyenv())

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

hist_data <- if (!is.null(component_surveys)) {
  purrr::map_dfr(component_surveys, ~get_hist_data(species, .x, hist_path, hist_cache))
} else {
  get_hist_data(species, survey_abbrev, hist_path, hist_cache)
}

fits_by_year <- purrr::map(eval_years, function(eval_year) {
  combined_data <- combine_hist_sim_data(sampled_data, hist_data, eval_year)

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
      profile = TRUE,
      newton_loops = 0L
    ),
    silent = T
  )

  list(
    eval_year = eval_year,
    data = combined_data,
    fit = fit
  )
})

names(fits_by_year) <- as.character(eval_years)

down_curve <- purrr::map_dfr(fits_by_year, function(x) {
  fit <- x$fit

  if (!is.null(fit$error) && isTRUE(fit$error)) {
    return(tibble(
      eval_year = x$eval_year,
      term = "restricted:year_covariate",
      estimate = NA_real_,
      std.error = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      significant = NA,
      error_msg = fit$message
    ))
  }

  tidy(fit, conf.int = TRUE) |>
    filter(term == "restricted:year_covariate") |>
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

print(down_curve)

down_curve
