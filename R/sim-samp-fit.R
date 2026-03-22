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

sp_list <- c("yelloweye rockfish")
mpa_rate_list <- c("10%")
# mpa_rate_list <- c("25%")

replicates <- 1:2

tag_list <- c("svc-rates-sd-0.01", "no-svc-rates")
for (sp in sp_list) {
  for (mpa_rate in mpa_rate_list) {
    for (tag in tag_list) {
      run_ssf(sp, mpa_rate, tag = tag, reps = replicates)
    }
  }
}

check_data_prep <- prep_full_timeseries("yelloweye rockfish",
  sampling_plan = "status_quo", rep_num = 1,
  sample_dir = here::here("data-generated", "2-sampled-data-svc-rates-sd-0.01"),
  hist_dir = here::here("data-generated", "cleaned-species-data"))

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

result_check <- readRDS(file.path(results_dir,
  paste0(sp_to_hyphens(sp),
    "-", round(lambda_val, 4),
    "-", sampling_plan,
    "-rep", sprintf("%03d", 2),
    "_", "eval", 2046, ".rds"
  )
))
glimpse(result_check)

# 4. POWER VISUALISATION -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs())

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
      # mpa_effect_label = first(mpa_effect_label),
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

get_cumul_power <- function(power_df, combo) {
  power_df0_n <- power_df0 |>
    add_count(!!!syms("eval_year"), name = "combo_n_reps")
  samples <- purrr::map_dfr(1:max(power_df$n_reps), \(x) {
    power_df0_n |>
      filter(combo_n_reps >= x) |>
      group_by(!!!syms(combo)) |>
      slice_sample(n = x, replace = FALSE) |>
      summarise_power(by = combo) |>
      mutate(n_samps = x)
  })
  # samples |> mutate(species = factor(species, levels = spp_levels))
}

# fit_dir <- here::here("data-generated", "test-fits-svc-rates")
# results_dir <- here::here("data-generated", "test-results-svc-rates")
# results_dir <- here::here("data-generated", "test-results-no-svc-rates")
# eval_years <- c(2030, 2034, 2038, 2042, 2046)

# fit_files <- list.files(fit_dir)
# test_fit <- readRDS(file.path(fit_dir, fit_files[2]))
# # sdmTMB::sanity(test_fit)
# results_files <- list.files(results_dir)
# test_res <- readRDS(file.path(results_dir, results_files[2]))

# lambda_val <- readRDS(file.path("data-generated", "recovery-rates-lambda.rds")) |>
#   filter(species == .env$species)
#   filter(case == "50% rate") |>
#   pull(lambda)

message("Reading results from ", basename(results_dir),
  "\n - lambda: ", round(lambda_val, 4),
  "\n - sampling plan: ", sampling_plan,
  "\n - replicate: ", min(replicates), " to ", max(replicates),
  "\n - eval year: ", paste(eval_years, collapse = ", "))

# res_tag <- "no-svc-rates"
res_tag <- "svc-rates-sd-0.01"
res_dir <- here::here("data-generated", paste0("4-fitted-results-", res_tag))
results_files <- list.files(results_dir,
  pattern = paste0(sp_to_hyphens(sp),
  "-", round(lambda_val, 4),
  "-", sampling_plan,
   "-rep[0-9]{3}_eval[0-9]{4}\\.rds$"),
  full.names = TRUE)

all_fitted_results <- purrr::map_dfr(results_files, function(f) {
  readRDS(f)
}) |>
  mutate(mpa_trend = log(lambda_val),
         plan = sampling_plan,
         tag = res_tag,
         species = sp
         )

# combo <- c("species", "mpa_trend", "tag", "plan", "eval_year")

combo <- c("eval_year")

# Replicate-level metrics
power_df0 <- all_fitted_results |>
  mutate(
    true_effect = mpa_trend,
    significant = !(ci_lower < 0 & ci_upper > 0), # Significance: CI doesn't include 0
    sign_correct = estimate * true_effect > 0, # more robust than assuming positive effect (e.g., lingcod)
    ratio_to_true = abs(estimate) / abs(true_effect)
  )

power_df <- power_df0 |>
  summarise_power(by = combo)


# Bias check on estimate -------------------------------------------------------
b1 <- power_df0 |>
  filter(replicate <= max(replicate)) |>
ggplot(data = _) +
  geom_point(aes(x = replicate, y = estimate)) +
  geom_hline(aes(yintercept = log(lambda_val)), colour = "red") +
  facet_grid(. ~ eval_year) +
  ggtitle("Bias check")
# b1

# Cumulative power plot --------------------------------------------------------
cpower <- get_cumul_power(power_df, combo)
c1 <- cpower |>
ggplot(data = ) +
  aes(x = n_samps, y = power, colour = mpa_effect_pct) +
  geom_point() +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  facet_grid(. ~ eval_year) +
  ggtitle("Cumulative power")
# c1
# Main power plot --------------------------------------------------------------
p1 <- power_df |>
ggplot(data = ) +
  aes(x = eval_year, y = power, colour = mpa_effect_pct) +
  geom_point() +
  geom_line(aes(group = mpa_effect_pct)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  ggtitle("Power")
# p1

(p1 / b1 / c1) + plot_annotation(title = tag)
