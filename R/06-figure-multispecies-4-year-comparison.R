library(dplyr)
library(ggplot2)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_SPECIES <- c("yelloweye rockfish", "lingcod", "quillback rockfish")
TARGET_SURVEY <- "HBLL"
TARGET_AR1 <- "fitted_AR1"
TARGET_MPA_EFFECT_LABEL <- "25%"
TARGET_EVAL_YEARS <- c(2038, 2042, 2046)
TARGET_EVAL_YEARS <- c(2042, 2046)
TARGET_PLANS <- c(
  "historical survey-year bootstrap",
  "historical survey-year bootstrap - no MPA every 2nd survey"
)

summarise_power <- function(power_df,
  by = c("species", "survey_abbrev", "sampling_plan", "mpa_effect_label", "eval_year")) {
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
      true_effect = first(true_effect),
      .groups = "drop"
    ) |>
    mutate(
      mpa_effect_pct = round(100 * true_effect, 2)
    )
}

results_dir <- here::here("data-generated", "power-results")
all_fitted_results0 <- readRDS(file.path(results_dir, "all-fitted-results.rds"))

rates_lu <- tibble::tibble(
  mpa_effect_label = factor(c("5%", "10%", "25%", "50%"), levels = c("5%", "10%", "25%", "50%")),
  sim_mpa_trend = round(exp(log(c(1.05, 1.10, 1.25, 1.5)) / 25), 3),
  true_effect = log(c(1.05, 1.10, 1.25, 1.5)) / 25
)

all_fitted_results <- all_fitted_results0 |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
  mutate(sim_mpa_trend = round(sim_mpa_trend, 3)) |>
  mutate(converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  left_join(rates_lu, by = "sim_mpa_trend")

plot_results <- all_fitted_results |>
  filter(
    # replicate %in% 1:100,
    species %in% TARGET_SPECIES,
    survey_abbrev == TARGET_SURVEY,
    sim_ar1_scenario == TARGET_AR1,
    mpa_effect_label == TARGET_MPA_EFFECT_LABEL,
    eval_year %in% TARGET_EVAL_YEARS,
    sampling_plan %in% TARGET_PLANS
  )

missing_plans <- setdiff(TARGET_PLANS, unique(plot_results$sampling_plan))
if (length(missing_plans) > 0) {
  stop(
    "Missing required sampling plans in all-fitted-results.rds: ",
    paste(missing_plans, collapse = ", "),
    call. = FALSE
  )
}

missing_species <- setdiff(TARGET_SPECIES, unique(plot_results$species))
if (length(missing_species) > 0) {
  stop(
    "Missing required species in all-fitted-results.rds: ",
    paste(missing_species, collapse = ", "),
    call. = FALSE
  )
}

power_df0 <- plot_results |>
  mutate(
    significant = !(ci_lower < 0 & ci_upper > 0),
    sign_correct = estimate * true_effect > 0,
    ratio_to_true = (exp(estimate * 25) - 1) / (exp(true_effect * 25) - 1)
  )

power_df <- summarise_power(power_df0)

power_plot <- power_df |>
  mutate(
    sampling_plan_label = recode(
      sampling_plan,
      "historical survey-year bootstrap" = "Status quo",
      "historical survey-year bootstrap - no MPA every 2nd survey" = "MPAs every 4 years"
    ),
    sampling_plan_label = factor(
      sampling_plan_label,
      levels = c("Status quo", "MPAs every 4 years")
    ),
    species_label = tools::toTitleCase(species),
    species_label = factor(
      species_label,
      levels = tools::toTitleCase(c("yelloweye rockfish", "lingcod", "quillback rockfish"))
    ),
    mpa_effect_label = factor(mpa_effect_label, levels = levels(rates_lu$mpa_effect_label))
  ) |>
  ggplot(aes(
    x = eval_year,
    y = power_signed,
    linetype = sampling_plan_label,
    shape = sampling_plan_label,
    group = sampling_plan_label
  )) +
  geom_hline(yintercept = 0.8, linetype = "solid", colour = "grey85") +
  geom_hline(yintercept = 0.6, linetype = "solid", colour = "grey85") +
  geom_hline(yintercept = 0.4, linetype = "solid", colour = "grey85") +
  geom_hline(yintercept = 0.2, linetype = "solid", colour = "grey85") +
  geom_line() +
  geom_point() +
  facet_wrap(~species_label) +
  scale_x_continuous(breaks = TARGET_EVAL_YEARS) +
  scale_linetype_manual(values = c("Status quo" = "solid", "MPAs every 4 years" = "dashed")) +
  scale_shape_manual(values = c("Status quo" = 19, "MPAs every 4 years" = 21)) +
  scale_y_continuous(limits = c(0, 1.), expand = expansion(mult = c(0, 0)), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Evaluation year",
    y = "Correctly signed power",
    linetype = "Sampling plan",
    shape = "Sampling plan",
    subtitle = "25% MPA increase over 25 years"
  ) +
  theme(legend.position = "top")

error_df <- power_df0 |>
  group_by(species, survey_abbrev, sampling_plan, eval_year) |>
  summarise(
    type_m = mean(ratio_to_true[significant & converged & sign_correct], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    sampling_plan_label = recode(
      sampling_plan,
      "historical survey-year bootstrap" = "Status quo",
      "historical survey-year bootstrap - no MPA every 2nd survey" = "MPAs every 4 years"
    ),
    sampling_plan_label = factor(
      sampling_plan_label,
      levels = c("Status quo", "MPAs every 4 years")
    ),
    species_label = tools::toTitleCase(species),
    species_label = factor(
      species_label,
      levels = tools::toTitleCase(c("yelloweye rockfish", "lingcod", "quillback rockfish"))
    )
  )

error_plot <- error_df |>
  tidyr::pivot_longer(c(type_m), names_to = "error_type", values_to = "value") |>
  mutate(error_type = recode(error_type, "type_m" = "Type M (exaggeration ratio)")) |>
  ggplot(aes(
    x = eval_year,
    y = value,
    linetype = sampling_plan_label,
    shape = sampling_plan_label,
    group = sampling_plan_label
  )) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ species_label) +
  scale_x_continuous(breaks = TARGET_EVAL_YEARS) +
  scale_y_continuous(limits = c(1, NA), expand = expansion(mult = c(0, .1))) +
  scale_linetype_manual(values = c("Status quo" = "solid", "MPAs every 4 years" = "dashed")) +
  scale_shape_manual(values = c("Status quo" = 19, "MPAs every 4 years" = 21)) +
  labs(
    x = "Evaluation year",
    y = "Type M (exageration ratio)",
    linetype = "Sampling plan",
    shape = "Sampling plan"
  ) +
  theme(legend.position = "top")

combined_plot <- patchwork::wrap_plots(power_plot, error_plot, ncol = 1, axes = "collect", guides = "collect", axis_titles = "collect") &
  theme(legend.position = "top") &
  theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(0.2, "lines"))
print(combined_plot)

ggsave(
  file.path(fig_dir, "multispecies-power-bootstrap-vs-no-mpa-every-2nd-survey-25pct.png"),
  plot = combined_plot,
  width = 6,
  height = 5
)

# Convergence diagnostics: check if lower convergence at 2038 inflates power_signed
power_df |>
  filter(eval_year == 2038) |>
  mutate(sampling_plan_label = recode(
    sampling_plan,
    "historical survey-year bootstrap" = "Status quo",
    "historical survey-year bootstrap - no MPA every 2nd survey" = "MPAs every 4 years"
  )) |>
  select(species, sampling_plan_label, n_reps, n_converged, convergence_rate, power_signed, power_allreps) |>
  arrange(species, sampling_plan_label) |>
  print(n = Inf)


patchwork::wrap_plots(
  power_plot + theme(legend.position = "right"),
  error_plot + theme(legend.position = "right"), ncol = 1, axes = "collect", guides = "collect", axis_titles = "collect") &
  # theme(legend.position = "right") &
  theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(0.2, "lines"))
print(combined_plot)

power_plot + theme(legend.position = "right")
ggsave(file.path(fig_dir, "presentations", "2026-05-05-CSAS-meeting",
"multispecies-power-bootstrap-vs-no-mpa-every-2nd-survey-25pct-legend-right.png"),
  width = 9.7, height = 3.7
)
