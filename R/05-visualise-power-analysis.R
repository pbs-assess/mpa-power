library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(ggrepel)

theme_set(gfplot::theme_pbs())

source(here::here("R", "00-fit-power-analysis-functions.R"))

hbll_ecp_encounter_cpue <- readRDS(here::here("data-generated", "overlays", "hbll-ecp-encounter-cpue.rds"))
glimpse(hbll_ecp_encounter_cpue)

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
  # filter(species == "yelloweye rockfish") |>
  # filter(species != "yelloweye rockfish" |
  #  (species == "yelloweye rockfish" & sim_mpa_trend %in% out_ye_rates)) |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish")) |>
  filter(is.na(error_msg) | error_msg != "Missing replicate in sampled data") |>
  # filter(species == "lingcod") |>
  mutate(mpa_trend = round(sim_mpa_trend, 3),
         converged = ifelse(sanity == "ok", TRUE, FALSE)) |>
  mutate(sampling_plan = factor(sampling_plan, levels = c("status quo", "MPAs every 4 years")),
         sim_ar1_scenario = factor(sim_ar1_scenario, levels = c("no_AR1", "moderate_AR1", "high_AR1"))) |>
  group_by(species) |>
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
filter_species <- c("yelloweye rockfish", "lingcod", "quillback rockfish")
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
    species %in% filter_species &
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
# ggsave(file.path(supp_dir, "bias-check-on-estimate-lingcod.pdf"), width = 9, height = 6.5)
# ggsave(file.path(supp_dir, "bias-check-on-estimate-lingcod.png"), width = 9, height = 6.5)

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
    # species %in% filter_species &
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
# ggsave(file.path(supp_dir, "cumulative-power-plot-all-species.pdf"), width = 15, height = 13)
# ggsave(file.path(supp_dir, "cumulative-power-plot-all-species.png"), width = 15, height = 13)

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
        # mpa_effect_label = factor(mpa_effect_label, levels = c("low", "moderate", "high"))
        )

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

# Power plot ------
dat |>
  ggplot() +
  aes(x = eval_year, y = power, colour = mpa_effect_label,
    group = interaction(mpa_effect_label, sim_ar1_scenario, sampling_plan)) +
  geom_line(aes(linetype = sampling_plan)) +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  # geom_label(data = dat |> filter(eval_year == 2048), aes(label = (mpa_effect_label)), x = 2048) +
  facet_grid(cols = vars(species), rows = vars(mpa_effect_label), scales = "free_y") +
  ggtitle("Power") +
  scale_colour_manual(values = trend_colours) +
  labs(colour = "MPA trend", linetype = "Sampling plan") +
  theme(legend.position = "bottom")
if (presentation) {
  ggsave(file.path(fig_dir, "power.png"), width = 12, height = 5)
} else {
  ggsave(file.path(fig_dir, "power.png"), width = 7.6, height = 7.6)
}

# Type M error rate ------
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
# if (presentation) {
#   ggsave(file.path(fig_dir, "type-m-error-rate.png"), width = 12, height = 5)
# } else {
  ggsave(file.path(fig_dir, "type-m-error-rate.png"), width = 7.6, height = 4.8)
# }

# Type S error rate
dat |>
  ggplot() +
  aes(x = interaction(eval_year, mpa_effect_label), y = type_s_error, colour = mpa_effect_label,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_point() +
  geom_segment(aes(xend = interaction(eval_year, mpa_effect_label), yend = 0)) +
  scale_colour_manual(values = trend_colours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels = \(x) sub("\\..*", "", x)) +  # keep only year
  facet_grid(cols = vars(species), rows = vars(sampling_plan))


# Comparing power with encounter rate and CPUE ---------------------------------
# encounter_cpue_df <- hbll_ecp_encounter_cpue |>
#   mutate(total_encounter_rate = (encounters_inside + encounters_outside) / (n_sets_inside + n_sets_outside),
#          total_cpue_rate = (total_catch_inside + total_catch_outside) / (n_sets_inside + n_sets_outside)) |>
#   filter(species_common_name %in% unique(power_df$species)) |>
#   select(species_common_name, encounter_rate_inside, encounter_rate_outside, cpue_inside, cpue_outside, total_encounter_rate, total_cpue_rate) |>
#   mutate(across(encounter_rate_inside:total_cpue_rate, \(x) round(x, 2)))

# test <- left_join(power_df, encounter_cpue_df, by = c("species" = "species_common_name"))
# glimpse(test)

# test_p_fn <- function(x_val, title) {
#   filter(test) |>
#     ggplot() +
#     aes(x = !!sym(x_val), y = power, colour = mpa_effect_label) +
#     geom_point() +
#     geom_line(aes(linetype = sampling_plan)) +
#     scale_colour_manual(values = trend_colours) +
#     facet_grid(rows = vars(sampling_plan), cols = vars(eval_year)) +
#     ggtitle(title)
# }

# p1 <- test_p_fn("encounter_rate_inside", "Encounter rate inside MPAs")
# p2 <- test_p_fn("encounter_rate_outside", "Encounter rate outside MPAs")
# p3 <- test_p_fn("cpue_inside", "CPUE inside MPAs")
# p4 <- test_p_fn("cpue_outside", "CPUE outside MPAs")
# p5 <- test_p_fn("total_encounter_rate", "Total encounter rate")
# p6 <- test_p_fn("total_cpue_rate", "Total CPUE rate")

# p1 / p3

# p1 / p3 / p5 / p6

# p2 / p4 / p5 / p6

test_p_fn_xaxis_all <- function(dat, x_val, xtitle) {
  # Create axis data for all year-species combinations
  all_axis_data <- dat |>
    distinct(species, eval_year, !!sym(x_val)) |>
    arrange(eval_year, !!sym(x_val))

  breaks_vals <- all_axis_data |> pull(!!sym(x_val))
  # Create labels with species name, year, and numeric value
  labels_vals <- sprintf("%s  %.2f",
                        all_axis_data |> pull(species) |> gsub("rockfish", "", x = _),
                        breaks_vals)

  ggplot(dat) +
    aes(x = !!sym(x_val), y = power, colour = mpa_effect_label) +
    geom_point() +
    geom_line(aes(linetype = sampling_plan,
                  group = interaction(mpa_effect_label, sampling_plan))) +
    scale_colour_manual(values = trend_colours) +
    scale_x_continuous(
      breaks = breaks_vals,
      labels = labels_vals,
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    facet_grid(rows = vars(sampling_plan), cols = vars(eval_year)) +
    labs(colour = "MPA trend", linetype = "Sampling plan",
         x = gsub("_", " ", x_val)) +
    xlab(xtitle) +
    ylab("Power") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      legend.position = "bottom"
    )
}

p1a <- test |> filter(eval_year %in% c(2030, 2038, 2046)) |>
  test_p_fn_xaxis_all("encounter_rate_inside", "Encounter Rate Inside MPAs")
p3a <- test |> filter(eval_year %in% c(2030, 2038, 2046)) |>
  test_p_fn_xaxis_all("cpue_inside", "CPUE Inside MPAs")

# Display variants
(p1a / p3a) + plot_layout(guides = "collect") & theme(legend.position = "top") &
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(file.path(fig_dir, "power-by-mpa-trend-sampling-intensity.png"),
  width = 8.7, height = 10.3)

# Convergence rate plot ---------------------------------------------------------
# @supp

dat |>
  ggplot() +
  aes(x = interaction(eval_year, mpa_effect_label), y = convergence_rate, colour = mpa_effect_label,
    group = interaction(mpa_effect_pct, sim_ar1_scenario, sampling_plan)) +
  geom_point() +
  geom_segment(aes(xend = interaction(eval_year, mpa_effect_label), yend = 0.9)) +
  scale_colour_manual(values = trend_colours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels = \(x) sub("\\..*", "", x)) +  # keep only year
  facet_grid(cols = vars(species), rows = vars(sampling_plan))


# ENCOUNTER RATE AND CPUE PLOTS ------------------------------------------------
make_hbll_ecp_tigure <- function(df, fill_limits = c(0, NA), padding = 0.5, digits = 2L,
                                show_species_labels = TRUE) {
  df$species_common_name <- stringr::str_to_title(df$species_common_name)
  df$species_common_name <- gsub(" North Pacific", "", df$species_common_name)
  sp <- df |>
    filter(zone == "Inside MPAs") |>
    arrange(val) |>
    pull(species_common_name) |>
    unique()
  if (length(sp) == 0) sp <- unique(df$species_common_name)
  df$species_common_name <- factor(df$species_common_name, levels = sp)
  df$zone <- factor(df$zone, levels = c("Inside MPAs", "Outside MPAs"))
  df$txt <- round(df$val, digits)
  g <- df |>
    ggplot(aes(x = zone, y = species_common_name)) +
    geom_tile(aes(fill = val), colour = "white") +
    geom_text(aes(label = txt), size = ggplot2::rel(3), hjust = 0.5, vjust = 0.5) +
    scale_fill_viridis_c(
      limits = fill_limits, begin = 0.15, end = 1, alpha = 0.9, option = "D", direction = 1
    ) +
    xlab("") +
    ylab("") +
    coord_cartesian(
      expand = FALSE,
      xlim = range(as.numeric(df$zone)) + c(-padding, padding),
      ylim = range(as.numeric(df$species_common_name)) + c(-padding - 0.5, padding + 0.5),
      clip = "off"
    ) +
    gfplot::theme_pbs() +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"))) +
    theme(
      plot.background = element_rect(fill = NA),
      plot.margin = margin(t = 0, r = -2, b = 0, l = -2),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(colour = "grey10"),
      axis.text.y = if (show_species_labels) element_text(hjust = 1.03) else element_blank(),
      legend.margin = margin(t = 0, r = 0.1, b = -15, l = 0),
      legend.position = "top",
      legend.box = "horizontal"
    ) +
    scale_x_discrete(position = "top")
  g
}

glimpse(hbll_ecp_encounter_cpue)

# Long data for tigures: encounter rate and CPUE by zone
hbll_ecp_encounter_long <-
  hbll_ecp_encounter_cpue |>
  select(species_common_name, encounter_rate_inside, encounter_rate_outside) |>
  tidyr::pivot_longer(
    cols = c(encounter_rate_inside, encounter_rate_outside),
    names_to = "zone", values_to = "val"
  ) |>
  mutate(zone = if_else(zone == "encounter_rate_inside", "Inside MPAs", "Outside MPAs"))

hbll_ecp_cpue_long <-
  hbll_ecp_encounter_cpue |>
  select(species_common_name, cpue_inside, cpue_outside) |>
  tidyr::pivot_longer(
    cols = c(cpue_inside, cpue_outside),
    names_to = "zone", values_to = "val"
  ) |>
  mutate(zone = if_else(zone == "cpue_inside", "Inside MPAs", "Outside MPAs"))

tig_spp <- hbll_ecp_encounter_long |>
  group_by(species_common_name) |>
  mutate(max_val = max(val)) |>
  filter(max_val > 0.4) |>
  pull(species_common_name) |>
  unique()

hbll_ecp_encounter_tigure <- hbll_ecp_encounter_long |>
  group_by(species_common_name) |>
  mutate(max_val = max(val)) |>
  filter(max_val > 0.4) |>
  make_hbll_ecp_tigure(
  fill_limits = c(0, 1), digits = 2L, show_species_labels = TRUE
  ) + labs(fill = "Encounter rate")

hbll_ecp_cpue_tigure <- hbll_ecp_cpue_long |> filter(species_common_name %in% tig_spp) |>
  make_hbll_ecp_tigure(
  fill_limits = c(0, NA), digits = 2L, show_species_labels = FALSE
) + labs(fill = "CPUE (catch per set)")

hbll_ecp_tigure_combined <-
  hbll_ecp_encounter_tigure + hbll_ecp_cpue_tigure +
  patchwork::plot_annotation(
    title = "HBLL E-CP species: inside vs outside MPAs (per unit effort)",
    theme = theme(plot.title = element_text(size = 11))
  )
hbll_ecp_tigure_combined

