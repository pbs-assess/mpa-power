library(dplyr)
library(ggplot2)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_SPECIES <- "quillback rockfish"
TARGET_SURVEY <- "HBLL"
TARGET_AR1 <- "fitted_AR1"
TARGET_MPA_EFFECT_LABEL <- "50%"
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
    species == TARGET_SPECIES,
    survey_abbrev == TARGET_SURVEY,
    sim_ar1_scenario == TARGET_AR1,
    mpa_effect_label == TARGET_MPA_EFFECT_LABEL,
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

power_df0 <- plot_results |>
  mutate(
    significant = !(ci_lower < 0 & ci_upper > 0),
    sign_correct = estimate * true_effect > 0
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
    mpa_effect_label = factor(mpa_effect_label, levels = levels(rates_lu$mpa_effect_label))
  ) |>
  ggplot(aes(
    x = eval_year,
    y = power_signed,
    linetype = sampling_plan_label,
    group = sampling_plan_label
  )) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  scale_linetype_manual(values = c(
    "Status quo" = "solid",
    "MPAs every 4 years" = "dashed"
  )) +
  scale_x_continuous(breaks = sort(unique(power_df$eval_year))) +
  scale_y_continuous(limits = c(-0.005, 1.005), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Evaluation year",
    y = "Correctly signed power",
    linetype = "Sampling plan",
    title = "Quillback rockfish",
    subtitle = "50% MPA increase over 25 years"
  ) +
  theme(legend.position = "top")

print(power_plot)

ggsave(
  file.path(fig_dir, "quillback-power-bootstrap-vs-no-mpa-every-2nd-survey.png"),
  plot = power_plot,
  width = 4,
  height = 4
)
