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
TARGET_EVAL_YEAR <- 2042
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
    species %in% TARGET_SPECIES,
    survey_abbrev == TARGET_SURVEY,
    sim_ar1_scenario == TARGET_AR1,
    mpa_effect_label == TARGET_MPA_EFFECT_LABEL,
    eval_year == TARGET_EVAL_YEAR,
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
    species_label = tools::toTitleCase(species),
    species_label = factor(
      species_label,
      levels = tools::toTitleCase(c("yelloweye rockfish", "lingcod", "quillback rockfish"))
    ),
    mpa_effect_label = factor(mpa_effect_label, levels = levels(rates_lu$mpa_effect_label))
  ) |>
  ggplot(aes(
    x = species_label,
    y = power_signed,
    shape = sampling_plan_label,
    group = sampling_plan_label
  )) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  scale_shape_manual(values = c(
    "Status quo" = 19,
    "MPAs every 4 years" = 21
  )) +
  scale_y_continuous(limits = c(0, 1.), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Species",
    y = "Correctly signed power",
    shape = "Sampling plan",
    title = "Power at evaluation year 2042",
    subtitle = "25% MPA increase over 25 years"
  ) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.1), axis.title.x.bottom = element_blank())

print(power_plot)

ggsave(
  file.path(fig_dir, "multispecies-power-bootstrap-vs-no-mpa-every-2nd-survey-25pct-2042.png"),
  plot = power_plot,
  width = 5,
  height = 4
)
