library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(ggrepel)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-fit-power-analysis-functions.R"))

# FIXME: missing!
# hbll_ecp_encounter_cpue <- readRDS(here::here("data-generated", "overlays", "hbll-ecp-encounter-cpue.rds"))
# glimpse(hbll_ecp_encounter_cpue)

presentation <- FALSE
if (presentation) {
  fig_dir <- here::here("figures", "presentations", "2026-03-10-survey-design-workshop")
} else {
  fig_dir <- here::here("figures")
}

supp_dir <- here::here("figures", "supplementary")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# Gelman and Carlin 2014
# - d_rep = estimate from replicate
# - d_true = true effect
# - SE_rep = standard error of estimate from replicate
# - alpha = 0.05 (significance level)
# Power = P(d_rep)
# Type S(ign) error rate: probability that d_rep has the incorrect sign given that it is significantly different from 0
# Type M(agnitude) error rate (exaggeration ratio): probability that the ratio of d_rep to the true effect size is greater than 1

# From code example:
# exaggeration <- mean(abs(estimate)[significant])/ d_true

# Okabe-Ito palette
trend_colours <- c(
  "5%"    = "#56B4E9",  # sky blue
  "10%" = "#E69F00",  # golden orange
  "25%" = "#D55E00",  # vermillion
  "50%" = "#2D7A2D"  # dark green
)

summarise_power <- function(power_df,
  by = c("species", "survey_abbrev", "sim_mpa_trend", "sim_ar1_scenario",
    "sampling_plan", "eval_year")) {
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
      type_m_error = mean(ratio_to_true[significant & converged], na.rm = TRUE),
      mean_estimate = mean(estimate[converged], na.rm = TRUE),
      true_effect = first(true_effect),
      mean_bias = mean(estimate[converged] - true_effect),
      .groups = "drop"
    ) |>
    mutate(
      mpa_effect_pct = round(100 * true_effect, 2)
    )
}

# Load results
# ------------------------------------------------------------------------------
# results_dir <- here::here("data-generated", "no-hc-power-results")

results_dir <- here::here("data-generated", "power-results")
combined_results <- combine_all_results(results_dir)

all_fitted_results0 <- readRDS(file.path(results_dir, "all-fitted-results.rds")) #|>
# select(-converged)

# rr <- readRDS(here::here('data-generated', 'recovery-rates-lambda.rds'))
# out_ye_rates <- filter(rr, species == "yelloweye rockfish") |> pull(lambda) |> round(3)

# all_fitted_results0 |>
#   filter(survey_abbrev == "HBLL") |>
#   pull(replicate) |>
#   max()

rates <- unique(all_fitted_results0$sim_mpa_trend)
rate_percents <- c("5%", "10%", "25%", "50%")
rates_lu <- data.frame(
  mpa_effect_label = factor(rate_percents, levels = rate_percents),
  sim_mpa_trend = round(c(exp(log(c(1.05, 1.10, 1.25, 1.5)) / 25)), 3)
)



all_fitted_results <- all_fitted_results0 |>
  filter(sampling_plan == "historical survey-year bootstrap") |>
  # filter(species == "yelloweye rockfish") |>
  # filter(species != "yelloweye rockfish" |
  #  (species == "yelloweye rockfish" & sim_mpa_trend %in% out_ye_rates)) |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
  # filter(species == "lingcod") |>
  mutate(mpa_trend = round(sim_mpa_trend, 3),
    converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  # mutate(sampling_plan = factor(sampling_plan, levels = c("status quo", "MPAs every 4 years")),
    # sim_ar1_scenario = factor(sim_ar1_scenario, levels = c("no_AR1", "moderate_AR1", "high_AR1"))) |>
  # group_by(species) |>
  left_join(rates_lu) |>
  # mutate(
  #   mpa_effect_label = case_when(
  #     n_distinct(sim_mpa_trend) == 1 ~ "single rate",
  #     n_distinct(sim_mpa_trend) == 2 ~ c("low", "high")[dense_rank(sim_mpa_trend)],
  #     TRUE ~ c("low", "moderate", "high")[dense_rank(sim_mpa_trend)]
  #   )
  # ) |>
  ungroup()

# Species ordered by increasing lowest mpa_trend used
spp_levels <- all_fitted_results |>
  distinct(species, mpa_trend) |>
  arrange(mpa_trend) |>
  pull(species) |>
  unique()

all_fitted_results |> glimpse()
filter(all_fitted_results, !converged) |>
  distinct(sanity, error_msg, converged) |>
  glimpse()
table(all_fitted_results$sanity)

# Calculate replicate-level metrics
# ----------------------------------
power_df0 <- all_fitted_results |>
  mutate(
    true_effect = log(sim_mpa_trend),
    significant = !(ci_lower < 0 & ci_upper > 0), # Significance: CI doesn't include 0
    sign_correct = estimate * true_effect > 0, # more robust than assuming positive effect (e.g., lingcod)
    ratio_to_true = abs(estimate) / abs(true_effect)
  )

# power_df0 |> glimpse()

# ignore these for now - sanity was set too strict for no Newton loops:
# power_df0 <- mutate(power_df0, converged = if_else(sanity == "gradients; all", TRUE, converged))

# Calculate scenario-level metrics
# ------------------------------------------------------------------------------
combo <- c("species", "survey_abbrev",
  "sim_mpa_trend", "sim_ar1_scenario",
  "sampling_plan",
  "eval_year")

power_df <- power_df0 |>
  summarise_power(
    by = c("species", "survey_abbrev", "sim_mpa_trend",
      "sim_ar1_scenario", "sampling_plan", "eval_year")
  )

# power_df |>
#   ggplot() +
#   aes(x = interaction(eval_year, species), y = convergence_rate,
#     colour = mpa_effect_pct, shape = sim_ar1_scenario) +
#   geom_point() +
#   facet_wrap(survey_abbrev ~ sampling_plan, ncol = 5) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   # ggtitle("Convergence rate by plan and MPA trend (all species pooled)")

# ------------------------------------------------------------------------------
# Supplementary figure to support choice of replicate number (stability of power analysis results)
# ------------------------------------------------------------------------------
# filter_species <- "lingcod"
filter_species <- c("yelloweye rockfish", "lingcod", "quillback rockfish")
filter_species <- c("yelloweye rockfish")
filter_species <- c("quillback rockfish")
# filter_species <- c("pacific halibut")
filter_survey <- "HBLL"
filter_ar1 <- "fitted_AR1"

# Caterpillar plot
combo <- c("species", "survey_abbrev",
  "sim_mpa_trend", "sim_ar1_scenario",
  "sampling_plan", "eval_year"
)

# Cumulative power plot - to check stability of power analysis results ---------
# Add replicate count per combo so we only sample up to each combo's n_reps
power_df0_n <- power_df0 |>
  add_count(!!!syms(combo), name = "combo_n_reps")

samples <- map_dfr(1:max(power_df$n_reps), \(x) {
  power_df0_n |>
    filter(combo_n_reps >= x) |>
    group_by(!!!syms(combo)) |>
    slice_sample(n = x, replace = FALSE) |>
    summarise_power(by = combo) |>
    mutate(n_samps = x)
})

# ------------------------------------------------------------------------------
# Main power plot
filter_species <- c("lingcod", "pacific spiny dogfish", "quillback rockfish", "yelloweye rockfish", "pacific halibut")
dat <- power_df |>
  filter(
    # species %in% filter_species,
    survey_abbrev == filter_survey
  ) |>
  mutate(sim_ar1_scenario = "fitted_AR1") |> # FIXME - should not need after next run
  mutate(sampling_plan = factor(sampling_plan, levels = c("status quo", "MPAs every 4 years")),
    sim_ar1_scenario = factor(sim_ar1_scenario, levels = c("no_AR1", "fitted_AR1")),
    species = factor(species, levels = spp_levels),
    # mpa_effect_label = factor(mpa_effect_label, levels = c("low", "moderate", "high"))
  )

# dat |>
#   filter(power > 0.8) |>
#   arrange(species, eval_year, sampling_plan) |>
#   distinct(eval_year, species, sim_ar1_scenario, sampling_plan,
#     mpa_effect_label, mpa_effect_pct)

year_threshold <- dat |>
  filter(sampling_plan == "status quo", sim_ar1_scenario == "fitted_AR1", power >= 0.8) |>
  group_by(species, mpa_effect_label) |>
  slice_min(eval_year, n = 1) |>
  select(species, mpa_effect_label, year_80pct_power = eval_year) |>
  ungroup()

power_summary <- dat |>
  filter(sampling_plan == "status quo") |> #, sim_ar1_scenario == "fitted_AR1") |>
  select(species, mpa_effect_label, eval_year, power_signed, n_converged, n_reps) |>
  pivot_wider(
    names_from = eval_year,
    values_from = power_signed,
    names_prefix = "power_"
  ) |>
  left_join(year_threshold, by = c("species", "mpa_effect_label")) |>
  mutate(year_80pct_power = replace_na(as.character(year_80pct_power), ">2046")) |>
  arrange(species, mpa_effect_label)

saveRDS(dat, here::here("data-generated", "power-results-df.rds"))

# Power plot ------
# g <- dat |>
#   ggplot() +
#   aes(x = eval_year, y = power_signed, colour = factor(mpa_effect_label),
#     group = interaction(mpa_effect_label, sim_ar1_scenario, sampling_plan)) +
#   # geom_line(aes(linetype = sampling_plan)) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
#   # geom_label(data = dat |> filter(eval_year == 2048), aes(label = (mpa_effect_label)), x = 2048) +
#   facet_grid(cols = vars(stringr::str_to_title(species)), rows = vars(mpa_effect_label), scales = "free_y") +
#   scale_colour_manual(values = trend_colours) +
#   labs(colour = "MPA trend", linetype = "Sampling plan") +
#   theme(legend.position = "bottom") + xlab("Evaluation year") + ylab("Correctly signed power") +
#   scale_y_continuous(limits = c(-0.005, 1.005), expand = expansion(mult = c(0, 00))) +
#   theme(
#     legend.position = "bottom",
#     panel.spacing.y = unit(0.5, "lines")
#   )
# print(g)
