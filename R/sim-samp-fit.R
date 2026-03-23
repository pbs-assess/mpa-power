source('R/00-utils.R')

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs())
source("R/00-run_ssf.R")

# sp_list <- c(
#   "yelloweye rockfish", "lingcod", "quillback rockfish",
#   "canary rockfish", "silvergray rockfish",
#   "north pacific spiny dogfish", "pacific halibut"
# )
# mpa_rate_list <- c("5%", "10%", "25%", "50%")

sp_list <- c("yelloweye rockfish", "lingcod", "quillback rockfish")
mpa_rate_list <- c("5%", "10%", "25%", "50%")
# mpa_rate_list <- c("25%")

replicates <- 1:15
fit_eval_years <- c(2030, 2034, 2038, 2042, 2046)

# tag_list <- c("svc-rates-sd-0.01", "no-svc-rates")
tag_list <- c("no-svc-rates")
for (sp in sp_list) {
  for (mpa_rate in mpa_rate_list) {
    for (tag in tag_list) {
      run_ssf(sp, mpa_rate, tag = tag, reps = replicates, eval_years = fit_eval_years, parallel = T)
    }
  }
}

# check_data_prep <- prep_full_timeseries("yelloweye rockfish",
#   sampling_plan = "status_quo", rep_num = 1,
#   sample_dir = here::here("data-generated", "2-sampled-data-svc-rates-sd-0.01"),
#   hist_dir = here::here("data-generated", "cleaned-species-data"))

if (FALSE) {
  # Check that years are correct
  p1 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = year, y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "Year", y = "Catch proportion", colour = "Historical")

  # Check that fyear is correct
  p2 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = as.numeric(fyear), y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "fyear", y = "Catch proportion", colour = "Historical")
  #
  # # Check that year_post_imp is correct
  p3 <- ggplot(data = check_data_prep) +
    geom_point(aes(x = year_post_imp, y = catch_prop, colour = factor(historical)), shape = 21) +
    scale_colour_manual(values = c("orange", "dodgerblue")) +
    facet_wrap(~ survey_abbrev, scales = "free_y") +
    labs(x = "Post-implementation year", y = "Catch proportion", colour = "Historical")
  (p1 / p2 / p3) + plot_annotation(title = "Data alignment / simulation + historical comparison") +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  # # ---
  #
  # meep()
}
analysis_species <- sp_list
analysis_mpa_rates <- mpa_rate_list
# analysis_tag <- "svc-rates-sd-0.01"
analysis_tag <- "no-svc-rates"
sampling_plan <- "status_quo"
eval_years <- c(2030, 2034, 2038, 2042, 2046)

rate_lookup <- tibble(
  mpa_rate = c("5%", "10%", "25%", "50%"),
  lambda_val = c(1.05, 1.10, 1.25, 1.50)^(1 / 25)
)

rate_colours <- c(
  "5%" = "#4C78A8",
  "10%" = "#59A14F",
  "25%" = "#F28E2B",
  "50%" = "#E15759"
)

summarise_power <- function(power_df,
  by = c("species", "mpa_rate", "sampling_plan", "eval_year")) {
  power_df |>
    group_by(!!!syms(by)) |>
    summarise(
      n_reps = n(),
      n_converged = sum(converged),
      convergence_rate = n_converged / n_reps,
      n_significant = sum(significant & converged),
      power = dplyr::if_else(n_converged > 0, n_significant / n_converged, NA_real_),
      power_allreps = n_significant / n_reps,
      type_s_error = dplyr::if_else(n_significant > 0,
        sum(!sign_correct & significant & converged) / n_significant,
        NA_real_),
      type_m_error = mean(ratio_to_true[significant & converged], na.rm = TRUE),
      mean_estimate = mean(estimate[converged], na.rm = TRUE),
      true_effect = first(true_effect),
      mean_bias = mean(estimate[converged] - true_effect),
      .groups = "drop"
    )
}

load_power_results <- function(species_vec, mpa_rate_vec, sampling_plan, results_dir, analysis_tag) {
  combo_grid <- tidyr::crossing(
    species = species_vec,
    mpa_rate = mpa_rate_vec
  ) |>
    left_join(rate_lookup, by = "mpa_rate")

  purrr::pmap_dfr(combo_grid, function(species, mpa_rate, lambda_val) {
    result_files <- list.files(
      results_dir,
      pattern = paste0(
        "^", sp_to_hyphens(species),
        "-", round(lambda_val, 4),
        "-", sampling_plan,
        "-rep[0-9]{3}_eval[0-9]{4}\\.rds$"
      ),
      full.names = TRUE
    )

    if (length(result_files) == 0) {
      warning("No results found for species=", species,
        ", mpa_rate=", mpa_rate,
        ", sampling_plan=", sampling_plan,
        ", tag=", analysis_tag)
      return(tibble())
    }

    purrr::map_dfr(result_files, readRDS) |>
      mutate(
        species = species,
        mpa_rate = mpa_rate,
        lambda_val = lambda_val,
        mpa_trend = log(lambda_val),
        sampling_plan = sampling_plan,
        tag = analysis_tag
      )
  })
}

results_dir <- here::here("data-generated", paste0("4-results-", analysis_tag))

message("Reading results from ", basename(results_dir),
  "\n - species: ", paste(analysis_species, collapse = ", "),
  "\n - rates: ", paste(analysis_mpa_rates, collapse = ", "),
  "\n - sampling plan: ", sampling_plan,
  "\n - replicates: ", min(replicates), " to ", max(replicates),
  "\n - eval years: ", paste(eval_years, collapse = ", "))

all_fitted_results <- load_power_results(
  species_vec = analysis_species,
  mpa_rate_vec = analysis_mpa_rates,
  sampling_plan = sampling_plan,
  results_dir = results_dir,
  analysis_tag = analysis_tag
)

if (nrow(all_fitted_results) == 0) {
  stop("No fitted result files matched the requested species/rate/tag/plan settings.")
}

result_check <- all_fitted_results |>
  arrange(species, mpa_rate, replicate, eval_year) |>
  slice(1)
glimpse(result_check)

power_df0 <- all_fitted_results |>
  filter(eval_year %in% eval_years, replicate %in% replicates) |>
  mutate(
    true_effect = mpa_trend,
    significant = !(ci_lower < 0 & ci_upper > 0),
    sign_correct = estimate * true_effect > 0,
    ratio_to_true = abs(estimate) / abs(true_effect)
  )

power_df <- power_df0 |>
  summarise_power(by = c("species", "mpa_rate", "sampling_plan", "eval_year")) |>
  mutate(
    species = factor(species, levels = analysis_species),
    mpa_rate = factor(mpa_rate, levels = analysis_mpa_rates)
  )

p1 <- power_df |>
  ggplot(aes(x = eval_year, y = power, colour = mpa_rate, group = mpa_rate)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 2.2) +
  facet_wrap(~ species, ncol = 2) +
  scale_colour_manual(values = rate_colours[analysis_mpa_rates], drop = FALSE) +
  scale_x_continuous(breaks = eval_years) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = paste("Power by species and recovery rate:", analysis_tag),
    subtitle = paste("Sampling plan:", sampling_plan),
    x = "Evaluation year",
    y = "Power",
    colour = "Recovery over 25 years"
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p1
