library(dplyr)
library(ggplot2)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

summarise_power <- function(power_df,
  by = c("species", "survey_abbrev", "mpa_effect_label", "eval_year")) {
  power_df |>
    group_by(!!!syms(by)) |>
    summarise(
      mpa_effect_label = first(mpa_effect_label),
      n_reps = n(),
      n_converged = sum(converged),
      convergence_rate = n_converged / n_reps,
      n_significant = sum(significant & converged),
      n_significant_signed = sum(significant & converged & sign_correct),
      power = n_significant / n_converged,
      power_signed = n_significant_signed / n_converged,
      power_allreps = n_significant / n_reps,
      type_s_error = sum(!sign_correct & significant & converged) / n_significant,
      type_m_error = mean(ratio_to_true[significant & converged & sign_correct], na.rm = TRUE),
      mean_estimate = mean(estimate[significant & converged & sign_correct], na.rm = TRUE),
      true_effect = first(true_effect),
      mean_bias = mean(estimate[converged] - true_effect),
      .groups = "drop"
    ) |>
    mutate(
      mpa_effect_pct = round(100 * true_effect, 2)
    )
}


# Load power results
results_dir <- here::here("data-generated", "power-results")
combined_results <- combine_all_results(results_dir)

all_fitted_results0 <- readRDS(file.path(results_dir, "all-fitted-results.rds"))

rates <- unique(all_fitted_results0$sim_mpa_trend)
rate_percents <- c("5%", "10%", "25%", "50%")
rates_lu <- data.frame(
  mpa_effect_label = factor(rate_percents, levels = rate_percents),
  sim_mpa_trend = round(c(exp(log(c(1.05, 1.10, 1.25, 1.5)) / 25)), 3),
  true_effect = log(c(1.05, 1.10, 1.25, 1.5)) / 25
)

all_fitted_results <- all_fitted_results0 |>
  filter(sampling_plan %in% c("historical survey-year bootstrap", "historical survey-year bootstrap - no MPA every 2nd survey")) |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  mutate(converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  left_join(rates_lu)# |>
  # filter(sampling_plan == "historical survey-year bootstrap")
  # filter(sampling_plan == "status quo")

all_fitted_results |> glimpse()
filter(all_fitted_results, !converged) |>
  distinct(sanity, error_msg, converged) |>
  glimpse()
table(all_fitted_results$sanity)

# Calculate replicate-level metrics
# ----------------------------------
power_df0 <- all_fitted_results |>
  mutate(
    significant = !(ci_lower < 0 & ci_upper > 0), # Significance: CI doesn't include 0
    sign_correct = estimate * true_effect > 0, # more robust than assuming positive effect (e.g., lingcod)
    # ratio_to_true = exp(estimate * 25) - exp(true_effect * 25)
    ratio_to_true = (exp(estimate * 25) - 1) / (exp(true_effect * 25) - 1)
  )

# power_df0 |> glimpse()

# Calculate scenario-level metrics
# ------------------------------------------------------------------------------
combo <- c("species", "mpa_effect_label", "sampling_plan", "eval_year")

power_df <- summarise_power(power_df0, by = combo)
spp_levels <- power_df |>
  filter(mpa_effect_label == "25%") |>
  group_by(species) |>
  slice(which.max(power_signed)) |>
  arrange(power_signed) |>
  pull(species)
# Order species by increasing max power at 25% MPA effect size
power_df <- power_df |> mutate(species = factor(species, levels = spp_levels))

saveRDS(power_df, here::here('data-generated', 'power-df-historical-sampling.rds'))
