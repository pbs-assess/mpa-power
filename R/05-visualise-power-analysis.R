library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-fit-power-analysis-functions.R"))

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

trend_colours <-
  c("low" = "#46327E",
  "moderate" = "#1F9E89",
  "high" = "#2D7A2D")

# Okabe-Ito palette
trend_colours <- c(
  "low"    = "#56B4E9",  # sky blue
  "moderate" = "#E69F00",  # golden orange
  "high"   = "#D55E00"   # vermillion
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
      power = n_significant / n_converged,
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

rr <- readRDS(here::here('data-generated', 'recovery-rates-lambda.rds'))
out_ye_rates <- filter(rr, species == "yelloweye rockfish") |> pull(lambda) |> round(3)

all_fitted_results0 |>
  filter(survey_abbrev == "HBLL") |>
  pull(replicate) |>
  max()

all_fitted_results <- all_fitted_results0 |>
  filter(species != "yelloweye rockfish" |
   (species == "yelloweye rockfish" & sim_mpa_trend %in% out_ye_rates)) |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
  # filter(species == "lingcod") |>
  mutate(mpa_trend = round(sim_mpa_trend, 3),
         converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  mutate(sampling_plan = factor(sampling_plan, levels = c("status quo", "MPAs every 4 years")),
         sim_ar1_scenario = factor(sim_ar1_scenario, levels = c("no_AR1", "moderate_AR1", "high_AR1"))) |>
  group_by(species) |>
  mutate(
    mpa_effect_label = case_when(
      n_distinct(sim_mpa_trend) == 1 ~ "single rate",
      n_distinct(sim_mpa_trend) == 2 ~ c("low", "high")[dense_rank(sim_mpa_trend)],
      TRUE ~ c("low", "moderate", "high")[dense_rank(sim_mpa_trend)]
    )
  ) |>
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
unique(all_fitted_results$sanity)

# Calculate replicate-level metrics
# ----------------------------------
power_df0 <- all_fitted_results |>
  mutate(
    true_effect = log(sim_mpa_trend),
    significant = !(ci_lower < 0 & ci_upper > 0), # Significance: CI doesn't include 0
    sign_correct = estimate * true_effect > 0, # more robust than assuming positive effect (e.g., lingcod)
    ratio_to_true = abs(estimate) / abs(true_effect)
  )

power_df0 |> glimpse()

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
filter_species <- "yelloweye rockfish"
filter_survey <- "HBLL"
filter_ar1 <- "fitted_AR1"

# Caterpillar plot
combo <- c("species", "survey_abbrev",
           "sim_mpa_trend", "sim_ar1_scenario",
           "sampling_plan", "eval_year"
)

# Bias check on estimate -------------------------------------------------------
power_df0 |>
  filter(
    species == filter_species &
    survey_abbrev == filter_survey #&
    # sim_ar1_scenario == filter_ar1
  ) |>
  group_by(!!!syms(combo)) |>
  mutate(id = row_number()) |>
  mutate(combo = paste(!!!syms(combo), collapse = " - ")) |>
  filter(replicate <= max(replicate)) |>
ggplot(data = ) +
  geom_point(aes(x = id, y = estimate)) +
  geom_hline(aes(yintercept = true_effect), colour = "red") +
  # geom_point(aes(x = id, y = ratio_to_true)) +
  # geom_hline(aes(yintercept = 1), colour = "red") +
  facet_grid(rows = vars(sampling_plan, sim_mpa_trend), cols = vars(eval_year)) +
  ggtitle(paste(filter_species, filter_survey, filter_ar1, sep = " - "))
ggsave(file.path(supp_dir, "bias-check-on-estimate-lingcod.pdf"), width = 9, height = 6.5)
ggsave(file.path(supp_dir, "bias-check-on-estimate-lingcod.png"), width = 9, height = 6.5)

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
samples |>
  filter(
    survey_abbrev == filter_survey &
    sampling_plan == "status quo" #&
    # species == filter_species &
    # sim_ar1_scenario == filter_ar1
  ) |>
  mutate(species = factor(species, levels = spp_levels)) |>
ggplot(data = _) +
  aes(x = n_samps, y = power, colour = mpa_effect_pct) +
  geom_point() +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(sampling_plan ~ interaction(eval_year, species), ncol = 5)
ggsave(file.path(supp_dir, "cumulative-power-plot-all-species.pdf"), width = 15, height = 13)
ggsave(file.path(supp_dir, "cumulative-power-plot-all-species.png"), width = 15, height = 13)

# Plot of raw estimates and confidence intervals --------------------------------
# power_df0 |>
#   filter(converged) |>
#   filter(sim_ar1_scenario == "no_AR1") |>
#   mutate(
#       # Create combined status variable
#       status = case_when(
#         !converged ~ "Not converged",
#         converged & significant ~ "Significant",
#         converged & !significant ~ "Not significant"
#       ),
#       status = factor(status, levels = c("Not significant", "Significant", "Not converged"))
#     ) |>
#   ggplot() +
#   geom_hline(aes(yintercept = log(sim_mpa_trend)), colour = "grey50", alpha = 0.5) +
#   aes(x = replicate, y = estimate, colour = status) +
#   geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width = 0.1), shape = 21) +
#   scale_colour_manual(
#   values = c(
#     "Not significant" = "grey30",      # Dark grey (not black for better contrast)
#     "Significant" = "#D55E00",         # Okabe-Ito orange (colorblind safe)
#     "Not converged" = "grey80"         # Light grey
#   )) +
#   facet_grid(rows = vars(eval_year), cols = vars(sim_mpa_trend, sampling_plan))



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
        mpa_effect_label = factor(mpa_effect_label, levels = c("low", "moderate", "high")))

dat |>
  filter(power > 0.8) |>
  arrange(species, eval_year, sampling_plan) |>
  distinct(eval_year, species, sim_ar1_scenario, sampling_plan,
    mpa_effect_label, mpa_effect_pct)


year_threshold <- dat |>
  filter(sampling_plan == "status quo", sim_ar1_scenario == "fitted_AR1", power >= 0.8) |>
  group_by(species, mpa_effect_label) |>
  slice_min(eval_year, n = 1) |>
  select(species, mpa_effect_label, year_80pct_power = eval_year) |>
  ungroup()

power_summary <- dat |>
  filter(sampling_plan == "status quo") |> #, sim_ar1_scenario == "fitted_AR1") |>
  select(species, mpa_effect_label, eval_year, power, n_converged, n_reps) |>
  pivot_wider(
    names_from = eval_year,
    values_from = power,
    names_prefix = "power_"
  ) |>
  left_join(year_threshold, by = c("species", "mpa_effect_label")) |>
  mutate(year_80pct_power = replace_na(as.character(year_80pct_power), ">2046")) |>
  arrange(species, mpa_effect_label)

sampling_comparison <- dat |>
  filter(
    sim_ar1_scenario == "fitted_AR1",
    eval_year %in% c(2038, 2046)  # Key evaluation years
  ) |>
  select(species, mpa_effect_label, eval_year, sampling_plan, power) |>
  pivot_wider(
    names_from = sampling_plan,
    values_from = power,
    names_prefix = "power_"
  ) |>
  mutate(
    power_difference = `power_status quo` - `power_MPAs every 4 years`
  ) |>
  arrange(species, mpa_effect_label, eval_year)

# View(power_summary)
# View(sampling_comparison)

dat |>
  ggplot() +
  aes(x = eval_year, y = power, colour = mpa_effect_label,
    group = interaction(mpa_effect_label, sim_ar1_scenario, sampling_plan)) +
  geom_line(aes(linetype = sampling_plan)) +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  # geom_label(data = dat |> filter(eval_year == 2048), aes(label = (mpa_effect_label)), x = 2048) +
  facet_grid(cols = vars(species), rows = vars(sim_ar1_scenario), scales = "free_y") +
  ggtitle("Power") +
  scale_colour_manual(values = trend_colours) +
  labs(colour = "MPA trend", linetype = "Sampling plan") +
  theme(legend.position = "bottom")
if (presentation) {
  ggsave(file.path(fig_dir, "power.png"), width = 12, height = 5)
} else {
  ggsave(file.path(fig_dir, "power.png"), width = 14, height = 5)
}

# Type M error rate
dat |>
  ggplot() +
  aes(x = eval_year, y = type_m_error, colour = mpa_effect_label,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_line(aes(linetype = sampling_plan)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = trend_colours) +
  facet_grid(cols = vars(species), rows = vars(sim_ar1_scenario)) +
  ggtitle("Type M error rate (exaggeration ratio)") +
  labs(colour = "MPA trend", linetype = "Sampling plan") +
  theme(legend.position = "bottom")
if (presentation) {
  ggsave(file.path(fig_dir, "type-m-error-rate.png"), width = 12, height = 5)
} else {
  ggsave(file.path(fig_dir, "type-m-error-rate.png"), width = 14, height = 5)
}

# Type S error rate
dat |>
  ggplot() +
  aes(x = eval_year, y = type_s_error, colour = mpa_effect_label,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_point() +
  facet_grid(cols = vars(species), rows = vars(sim_ar1_scenario))

# Convergence rate plot ---------------------------------------------------------
# FIXME these plots are broken/not sure if we even need
# @supp
dat |>
  ggplot() +
  aes(x = interaction(eval_year, species, sampling_plan), y = convergence_rate, colour = mpa_effect_pct,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_point() +
  facet_grid(rows = vars(species), cols = vars(sim_ar1_scenario)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

filter(power_df0, species == "yelloweye rockfish") |>
  tidyr::unite(combo, !!!syms(combo), sep = " - ", remove = FALSE) |>
  group_by(!!!syms(combo), combo) |>
  summarise(spatial_collapsed = sum(fit_spatial == "off") / max(replicate),
            spatiotemporal_collapsed = sum(fit_spatiotemporal == "off") / max(replicate)) |>
ggplot() +
  aes(x = eval_year, y = spatiotemporal_collapsed, colour = sim_mpa_trend,
    group = interaction(sim_mpa_trend, sim_ar1_scenario, sampling_plan)) +
  geom_line(aes(linetype = sampling_plan)) +
  geom_point() +
  facet_wrap(sim_ar1_scenario ~ sampling_plan)

power_df |>
  filter(species == filter_species, survey_abbrev == filter_survey) |>
  mutate(sampling_plan = factor(sampling_plan, levels = c("status quo", "MPAs every 4 years")),
        sim_ar1_scenario = factor(sim_ar1_scenario, levels = c("no_AR1", "moderate_AR1", "high_AR1"))) |>
  ggplot() +
  aes(x = eval_year, y = power, colour = mpa_effect_pct,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_line(aes(linetype = sampling_plan)) +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  facet_grid(rows = vars(species), cols = vars(sim_ar1_scenario))


power_df0 |>
  filter(converged) |>
  ggplot(aes(x = ratio_to_true, fill = factor(mpa_trend))) +
  geom_histogram() +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Blues")[c(5, 6, 7, 8)]) +
  facet_grid(survey_abbrev ~ plan + mpa_trend, scales = "free")

power_df |>
  # mutate(plan = stringr::str_replace(plan, "status quo ", ""),
  #         plan = as.numeric(plan)) |>
  ggplot(aes(x = mpa_trend, y = power, colour = factor(plan))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey30") +
  scale_y_continuous(labels = scales::percent_format()) +
  # scale_colour_manual(values = RColorBrewer::brewer.pal(9, "Blues")[c(3, 5, 7, 9)]) +
  labs(x = "MPA effect (% per year)", y = "Power",
        colour = "Sampling intensity"#,
      #  caption = "Dashed line = 80% power threshold"
        ) +
  theme(legend.position = "top") +
  facet_grid(plan ~ survey_abbrev)
ggsave(file.path(pres_dir, "power-by-mpa-trend-sampling-intensity.png"), width = 9.5, height = 3.8)


power_df |>
    mutate(
      # Convert to percentage increase per year
      mpa_effect_pct = 100 * (mpa_trend - 1),
      # Clean up plan names for facet labels
      plan_clean = stringr::str_replace_all(plan, c(
        "status-quo-" = "",
        "status-quo" = "Status Quo",
        "-" = " ",
        "effort" = "Effort",
        "every" = "Every",
        "years" = "Years",
        "no sampling in mpas" = "No Sampling in MPAs"
      )),
      plan_clean = stringr::str_to_title(plan_clean)
    ) |>
    ggplot(aes(x = mpa_effect_pct, y = power)) +
    geom_line(linewidth = 1.5, colour = "#4575b4") +
    geom_point(size = 4, colour = "#4575b4") +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50", linewidth = 0.8) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_x_continuous(breaks = c(1.5, 3, 5)) +
    labs(
      title = "Yelloweye rockfish",
      # subtitle = "Yelloweye Rockfish across different survey and sampling designs",
      x = "MPA Recovery Rate (% increase per year)",
      y = "Power"#,
      # caption = "Dashed line = 80% power threshold"
    ) +
    facet_grid(plan_clean ~ survey_abbrev) +
    gfplot::theme_pbs(base_size = 16) #+
    # theme(
    #   strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
    #   strip.text.x = element_text(face = "bold", size = 10),
    #   panel.grid.minor = element_blank(),
    #   plot.title = element_text(face = "bold", size = 13),
    #   plot.subtitle = element_text(size = 10, colour = "grey30"),
    #   plot.caption = element_text(hjust = 0, colour = "grey50")
    # )
ggsave(file.path(pres_dir, "power-by-mpa-trend-sampling-intensity.png"),
width = 20.42553, height = 12.01053)



power_df |>
  mutate(plan = stringr::str_replace(plan, "status quo ", ""),
          plan = as.numeric(plan)) |>
  ggplot(aes(x = plan, y = type_S, colour = factor(mpa_trend))) +
  facet_wrap(~ survey_abbrev) +
  geom_line() +
  geom_point(size = 2, position = position_dodge(width = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(9, "Blues")[c(5, 7, 9)])

power_df |>
  mutate(plan = stringr::str_replace(plan, "status quo ", ""),
          plan = as.numeric(plan)) |>
  ggplot(aes(x = plan, y = type_m, colour = factor(mpa_trend))) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(9, "Blues")[c(5, 7, 9)]) +
  facet_wrap(~ survey_abbrev)

