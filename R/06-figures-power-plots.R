library(dplyr)
library(ggplot2)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-utils.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

fig_dir <- here::here("figures")

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

all_fitted_results0 <- readRDS(file.path(results_dir, "all-fitted-results.rds")) #|>

rates <- unique(all_fitted_results0$sim_mpa_trend)
rate_percents <- c("5%", "10%", "25%", "50%")
rates_lu <- data.frame(
  mpa_effect_label = factor(rate_percents, levels = rate_percents),
  sim_mpa_trend = round(c(exp(log(c(1.05, 1.10, 1.25, 1.5)) / 25)), 3),
  true_effect = log(c(1.05, 1.10, 1.25, 1.5)) / 25
)

all_fitted_results <- all_fitted_results0 |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  mutate(converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  left_join(rates_lu)

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
combo <- c("species", "mpa_effect_label", "eval_year")

power_df <- summarise_power(power_df0, by = combo)
spp_levels <- power_df |>
  filter(mpa_effect_label == "25%") |>
  group_by(species) |>
  slice(which.max(power_signed)) |>
  arrange(power_signed) |>
  pull(species)
# Order species by increasing max power at 25% MPA effect size
power_df <- power_df |> mutate(species = factor(species, levels = spp_levels))

# d <- filter(power_df, species == "yelloweye rockfish", mpa_effect_label == "10%", eval_year == "2030")
# # names(d)

# # exp(d$mean_estimate * 25) * 100 - exp(d$true_effect * 25) * 100

# exp(d$mean_estimate * 25)
# exp(d$true_effect * 25)
# exp(d$mean_estimate * 25) / exp(d$true_effect * 25)
# (exp(d$mean_estimate * 25) *100) / (exp(d$true_effect * 25) * 100)

# (exp(d$mean_estimate * 25) *100) / (exp(d$true_effect * 25) * 100)

# ((exp(d$mean_estimate * 25) * 100) - 100) /
# ((exp(d$true_effect * 25) * 100) - 100)

# ((exp(d$mean_estimate * 25)) - 1) / ((exp(d$true_effect * 25)) - 1)

# 10 * 24


# exp(d$mean_estimate * 25) - exp(d$true_effect * 25)
#


# ------------------------------------------------------------------------------
# Main power plot
year_threshold <- power_df |>
  filter(power >= 0.8) |>
  group_by(species, eval_year, mpa_effect_label) |>
  slice_min(eval_year, n = 1) |>
  select(species, mpa_effect_label, year_80pct_power = eval_year) |>
  ungroup()

# Power plot ------
power_df |>
  ggplot() +
  aes(x = eval_year, y = power_signed, colour = mpa_effect_label) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ species, labeller = as_labeller(stringr::str_to_title),
    nrow = 2) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  scale_x_continuous(breaks = unique(power_df$eval_year)) +
  labs(colour = "Recovery over 25 years") +
  theme(legend.position = "top",
        panel.spacing = unit(1, "lines")) +
  labs(x = "Evaluation year", y = "Correctly signed power") +
  scale_y_continuous(limits = c(-0.005, 1.005), expand = expansion(mult = c(0, 00)))
ggsave(file.path(fig_dir, "main-power-plot.png"), width = 6.2, height = 5)

# Type M error plot ------------------------------------------------------------
power_df |>
  ggplot() +
  aes(x = eval_year, y = type_m_error, colour = mpa_effect_label) +
  geom_point() +
  geom_line() +
  facet_wrap(~ species, labeller = as_labeller(stringr::str_to_title),
    nrow = 2) +
  scale_colour_viridis_d(option = "plasma", end = 0.85) +
  scale_x_continuous(breaks = unique(power_df$eval_year)) +
  geom_hline(yintercept = 1) +
  labs(colour = "Recovery over 25 years") +
  theme(legend.position = "top",
        panel.spacing = unit(1, "lines")) +
  scale_y_log10(limits = c(1, NA), expand = expansion(mult = c(0, 0.05)), breaks = c(1, 2, 5, 10, 30)) +
  labs(x = "Evaluation year", y = "Multiplicative magnitude error\non the 25-year percent increase")
ggsave(file.path(fig_dir, "type-m-error-plot.png"), width = 6.2, height = 5)

# Type S error plot ------------------------------------------------------------
# Current option
power_df |>
  mutate(
    species = factor(stringr::str_to_title(species), levels = stringr::str_to_title(spp_levels)),
    has_error = type_s_error > 0
  ) |>
  ggplot(aes(x = factor(eval_year), y = type_s_error, colour = mpa_effect_label)) +
    geom_segment(aes(xend = factor(eval_year), yend = 0, alpha = has_error)) +
    geom_point(aes(size = has_error, alpha = has_error)) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 2), guide = "none") +
    facet_grid(rows = vars(mpa_effect_label), cols = vars(species)) +
    scale_colour_viridis_d(option = "plasma", end = 0.85, guide = "none") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.05),
                        expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Evaluation year", y = "Type S error rate") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.spacing = unit(0.5, "lines")
    )
ggsave(file.path(fig_dir, "type-s-error-plot.png"), width = 6.2, height = 5)


# ------------------------------------------------------------------------------
# Supplementary figure to support choice of replicate number (stability of power analysis results)
# ------------------------------------------------------------------------------
# filter_species <- c("yelloweye rockfish", "lingcod", "quillback rockfish")
filter_species <- c("yelloweye rockfish")
filter_species <- c("lingcod")

# Bias check on estimate -------------------------------------------------------
power_df0 |>
  # filter(species %in% filter_species) |>
  filter(mpa_effect_label == "25%") |>
  group_by(!!!syms(combo)) |>
  mutate(id = row_number()) |>
  mutate(combo = paste(!!!syms(combo), collapse = " - ")) |>
  filter(replicate <= max(replicate)) |>
  ggplot() +
  geom_point(aes(x = id, y = estimate)) +
  geom_hline(aes(yintercept = true_effect), colour = "red") +
  facet_grid(rows = vars(species), cols = vars(eval_year)) +
  labs(x = "Replicate", y = "Estimate") +
  ggtitle(paste0("25% recovery over 25 years"))
ggsave(file.path(fig_dir, paste0("bias-check-on-estimate-", sp_to_hyphens(filter_species), ".png")),
  width = 6.2, height = 4.6)
# ggsave(file.path(supp_dir, "bias-check-on-estimate-lingcod.png"), width = 9, height = 6.5)

# Cumulative power plot - to check stability of power analysis results ---------
# Add replicate count per combo so we only sample up to each combo's n_reps
power_df0_n <- power_df0 |>
  add_count(!!!syms(combo), name = "combo_n_reps")
glimpse(power_df0_n)

samples <- purrr::map_dfr(1:max(power_df0_n$combo_n_reps), \(x) {
  power_df0_n |>
    filter(combo_n_reps >= x) |>
    group_by(!!!syms(combo)) |>
    slice_sample(n = x, replace = FALSE) |>
    summarise_power(by = combo) |>
    mutate(n_samps = x)
})
samples |>
  # filter(species %in% filter_species) |>
  # mutate(species = factor(species, levels = spp_levels)) |>
  ggplot(data = _) +
  aes(x = n_samps, y = power_signed, colour = factor(mpa_effect_label)) +
  geom_point(size = 0.8) +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
   scale_colour_viridis_d(option = "plasma", end = 0.85) +
  facet_grid(cols = vars(eval_year), rows = vars(species), labeller = as_labeller(stringr::str_to_title)) +
  labs(x = "Number of replicates", y = "Power", colour = "Recovery over 25 years") +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "cumulative-power-plot-all-species.png"), width = 8, height = 9.5)
