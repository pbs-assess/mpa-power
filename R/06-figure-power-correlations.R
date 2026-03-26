library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(mgcv)
library(rpart)
library(ranger)

power_df <- readRDS("data-generated/power-results-df.rds")
theta <- readRDS("data-generated/fit_characteristics.rds")
ar1_theta <- readRDS("data-generated/ar1-parameters.rds") |>
  rename(survey = survey_abbrev)

theta <- left_join(theta, ar1_theta, by = c("species", "survey"))

# `power_df` is at the combined HBLL level, while `theta` is survey-specific.
# Aggregate conditioning-model characteristics to species-level summaries before joining.
theta_species <- theta |>
  group_by(species) |>
  summarise(
    n_surveys = n(),
    across(
      c(
        rho_V,
        sigma_V,
        encountered_count_per_year,
        encountered_rate_avg,
        encountered_rate_restricted,
        encountered_count_per_year_restricted,
        range,
        phi,
        sigma_O,
        sigma_E
      ),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(
    restricted_share_mean = encountered_count_per_year_restricted_mean /
      encountered_count_per_year_mean
  )

theta_species$species[theta_species$species == "north pacific spiny dogfish"] <- "pacific spiny dogfish"

diag_df <- power_df |>
  left_join(theta_species, by = "species") |>
  mutate(
    eval_year = factor(eval_year),
    mpa_effect_label = factor(mpa_effect_label),
    mpa_effect_pct = factor(mpa_effect_pct)
  )
  # filter(!mpa_effect_label %in% c("5%"))
  # filter(!eval_year %in% c(2030))
  # filter(!species %in% c("quillback rockfish"))

temp <- diag_df[,c("encountered_rate_restricted_mean", "encountered_rate_avg_mean", "encountered_count_per_year_restricted_mean", "sigma_O_mean", "sigma_E_mean", "phi_mean", "range_mean", "sigma_V_mean")] |>
  distinct()
GGally::ggpairs(temp)
ggsave("figures/pairs-covariates-power.png", width = 10, height = 10)

make_corr_plot <- function(column, .xlab = NULL) {
  corr_df <- diag_df |>
    group_by(mpa_effect_label, eval_year) |>
    summarise(
      corr = if (sum(!is.na({{column}}) & !is.na(power_signed)) >= 2) {
        cor({{column}}, power_signed, use = "complete.obs")
      } else {
        NA_real_
      },
      label = sprintf("r = %.2f", corr),
      .groups = "drop"
    )

  g <- diag_df |>
    mutate(species = stringr::str_to_title(species)) |>
    ggplot(aes({{column}}, power_signed, colour = species)) + geom_point() +
    facet_grid(mpa_effect_label~eval_year, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    geom_smooth(
      method = "lm", se = FALSE, formula = y ~ x,
      mapping = aes(group = NULL, colour = NULL), colour = "grey50", lwd = 0.7
    ) +
    geom_text(
      data = corr_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.1,
      size = 3,
      vjust = 1.3,
      colour = "grey50"
    ) +
    ylab("Power (correctly signed)") +
    labs(colour = "Species") +
    theme(legend.position = "top")
  if (!is.null(.xlab)) g <- g + xlab(.xlab)
  g
}

# make_corr_plot(encountered_count_per_year_mean)
make_corr_plot(encountered_rate_restricted_mean, "Mean countered rate in MPA network")
# make_corr_plot(encountered_rate_avg_mean)
ggsave("figures/power-encounter-rate-correlation.png", width = 7, height = 6)

make_corr_plot(sigma_V_mean, "Mean interannual variability")
ggsave("figures/power-sigmaV-correlation.png", width = 7, height = 6)

make_corr_plot(log(phi_mean), "Log observation error precision")
ggsave("figures/power-phi-correlation.png", width = 7, height = 6)

make_corr_plot(sigma_O_mean, "Spatial random field SD")
ggsave("figures/power-sigmaO-correlation.png", width = 7, height = 6)

make_corr_plot(sigma_E_mean, "Spatiotemporal random field SD")
ggsave("figures/power-sigmaE-correlation.png", width = 7, height = 6)



# diag_df |>
#   ggplot(aes(rho_V_mean, power_signed, colour = species)) + geom_point() +
#   facet_grid(mpa_effect_label~eval_year, scales = "free_y") +
#   scale_color_brewer(palette = "Dark2") +
#   geom_smooth(method = "lm", se = FALSE, mapping = aes(group = NULL, colour = NULL), colour = "grey50", lwd = 0.7)





#
# numeric_predictors <- c(
#   "encountered_count_per_year_mean",
#   "encountered_rate_avg_mean",
#   "encountered_rate_restricted_mean",
#   "encountered_count_per_year_restricted_mean",
#   "restricted_share_mean",
#   "range_mean",
#   "phi_mean",
#   "sigma_O_mean",
#   "sigma_E_mean"
# )
#
# cor_summary <- purrr::map_dfr(numeric_predictors, function(v) {
#   tibble(
#     variable = v,
#     spearman = cor(diag_df[[v]], diag_df$power_signed, method = "spearman", use = "complete.obs"),
#     pearson = cor(diag_df[[v]], diag_df$power_signed, method = "pearson", use = "complete.obs")
#   )
# }) |>
#   arrange(desc(abs(spearman)))
#
# print(cor_summary, n = Inf)
#
# # A small GAM is useful here because the signal is likely driven by combinations,
# # not by single predictors on their own.
# gam_df <- diag_df |>
#   select(
#     power_signed,
#     eval_year,
#     mpa_effect_label,
#     encountered_rate_restricted_mean,
#     restricted_share_mean,
#     phi_mean,
#     sigma_O_mean
#   ) |>
#   drop_na()
#
# power_gam <- mgcv::gam(
#   power_signed ~
#     eval_year +
#     mpa_effect_label +
#     s(encountered_rate_restricted_mean, k = 4) +
#     s(phi_mean, k = 4) +
#     ti(encountered_rate_restricted_mean, sigma_O_mean, k = c(4, 4)),
#   data = gam_df,
#   method = "REML"
# )
# par(mfrow = c(2, 2))
# plot(power_gam)
#
# print(summary(power_gam))
#
# tree_df <- diag_df |>
#   select(
#     power_signed,
#     eval_year,
#     mpa_effect_label,
#     encountered_count_per_year_mean,
#     encountered_rate_restricted_mean,
#     restricted_share_mean,
#     range_mean,
#     phi_mean,
#     sigma_O_mean,
#     sigma_E_mean
#   ) |>
#   drop_na()
#
# power_tree <- rpart::rpart(
#   power_signed ~ .,
#   data = tree_df,
#   method = "anova",
#   control = rpart.control(cp = 0.01, minsplit = 8, maxdepth = 4)
# )
#
# printcp(power_tree)
# print(power_tree$variable.importance)
#
# rf_df <- diag_df |>
#   select(
#     power_signed,
#     eval_year,
#     mpa_effect_label,
#     encountered_rate_restricted_mean,
#     range_mean,
#     phi_mean,
#     sigma_O_mean,
#     sigma_E_mean,
#     restricted_share_mean
#   ) |>
#   drop_na()
#
# power_rf <- ranger(
#   power_signed ~ eval_year + mpa_effect_label + encountered_rate_restricted_mean +
#     range_mean + phi_mean + sigma_O_mean + sigma_E_mean + restricted_share_mean,
#   data = rf_df,
#   importance = "permutation",
#   num.trees = 1000,
#   seed = 42
# )
#
# print(power_rf)
# print(sort(power_rf$variable.importance, decreasing = TRUE))
#
# rf_importance <- tibble(
#   variable = names(power_rf$variable.importance),
#   importance = unname(power_rf$variable.importance)
# ) |>
#   arrange(desc(importance))
#
# rf_predict_mean <- function(model, data) {
#   mean(predict(model, data = data)$predictions, na.rm = TRUE)
# }
#
# make_pdp_1d <- function(model, data, var_name, n = 40) {
#   grid_vals <- seq(min(data[[var_name]], na.rm = TRUE), max(data[[var_name]], na.rm = TRUE), length.out = n)
#
#   tibble(x = grid_vals) |>
#     mutate(
#       variable = var_name,
#       y = purrr::map_dbl(x, function(v) {
#         new_data <- data
#         new_data[[var_name]] <- v
#         rf_predict_mean(model, new_data)
#       })
#     )
# }
#
# make_pdp_2d <- function(model, data, var1_name, var2_name, n = 30) {
#   grid <- tidyr::expand_grid(
#     x = seq(min(data[[var1_name]], na.rm = TRUE), max(data[[var1_name]], na.rm = TRUE), length.out = n),
#     y = seq(min(data[[var2_name]], na.rm = TRUE), max(data[[var2_name]], na.rm = TRUE), length.out = n)
#   )
#
#   grid |>
#     mutate(
#       var1 = var1_name,
#       var2 = var2_name,
#       z = purrr::pmap_dbl(list(x, y), function(x, y) {
#         new_data <- data
#         new_data[[var1_name]] <- x
#         new_data[[var2_name]] <- y
#         rf_predict_mean(model, new_data)
#       })
#     )
# }
#
# make_binned_heatmap_data <- function(data, x_var, y_var, n_bins = 4) {
#   x_bin_name <- paste0(x_var, "_bin")
#   y_bin_name <- paste0(y_var, "_bin")
#
#   data |>
#     filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]), !is.na(power_signed)) |>
#     mutate(
#       !!x_bin_name := ntile(.data[[x_var]], n_bins),
#       !!y_bin_name := ntile(.data[[y_var]], n_bins)
#     ) |>
#     group_by(eval_year, mpa_effect_label, .data[[x_bin_name]], .data[[y_bin_name]]) |>
#     summarise(
#       mean_power_signed = mean(power_signed, na.rm = TRUE),
#       n = n(),
#       .groups = "drop"
#     ) |>
#     rename(x_bin = all_of(x_bin_name), y_bin = all_of(y_bin_name)) |>
#     mutate(x_var = x_var, y_var = y_var)
# }
#
# plot_binned_heatmap <- function(heatmap_df, x_label, y_label, title = NULL) {
#   ggplot(heatmap_df, aes(x_bin, y_bin, fill = mean_power_signed)) +
#     geom_tile() +
#     facet_grid(eval_year ~ mpa_effect_label) +
#     scale_fill_viridis_c() +
#     labs(
#       x = x_label,
#       y = y_label,
#       fill = "Mean signed power",
#       title = title
#     )
# }
#
# pdp_range <- make_pdp_1d(power_rf, rf_df, "range_mean")
# pdp_restricted_share <- make_pdp_1d(power_rf, rf_df, "restricted_share_mean")
# pdp_sigma_O <- make_pdp_1d(power_rf, rf_df, "sigma_O_mean")
# pdp_encountered_restricted <- make_pdp_1d(power_rf, rf_df, "encountered_rate_restricted_mean")
# pdp_surface <- make_pdp_2d(power_rf, rf_df, "sigma_O_mean", "encountered_rate_restricted_mean")
#
# importance_plot <- ggplot(rf_importance, aes(reorder(variable, importance), importance)) +
#   geom_col(fill = "grey35") +
#   coord_flip() +
#   labs(
#     x = NULL,
#     y = "Permutation importance",
#     title = "Random forest variable importance"
#   )
#
# pdp_range_plot <- ggplot(pdp_range, aes(x, y)) +
#   geom_line(linewidth = 1) +
#   labs(
#     x = "range_mean",
#     y = "Mean RF predicted power",
#     title = "Partial dependence: range_mean"
#   )
#
# pdp_restricted_share_plot <- ggplot(pdp_restricted_share, aes(x, y)) +
#   geom_line(linewidth = 1) +
#   labs(
#     x = "restricted_share_mean",
#     y = "Mean RF predicted power",
#     title = "Partial dependence: restricted_share_mean"
#   )
#
# pdp_sigma_O_plot <- ggplot(pdp_sigma_O, aes(x, y)) +
#   geom_line(linewidth = 1) +
#   labs(
#     x = "sigma_O_mean",
#     y = "Mean RF predicted power",
#     title = "Partial dependence: sigma_O_mean"
#   )
#
# pdp_encountered_restricted_plot <- ggplot(pdp_encountered_restricted, aes(x, y)) +
#   geom_line(linewidth = 1) +
#   labs(
#     x = "encountered_rate_restricted_mean",
#     y = "Mean RF predicted power",
#     title = "Partial dependence: encountered_rate_restricted_mean"
#   )
#
# pdp_surface_plot <- ggplot(pdp_surface, aes(x, y, fill = z)) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   labs(
#     x = "sigma_O_mean",
#     y = "encountered_rate_restricted_mean",
#     fill = "Mean RF\npredicted power",
#     title = "2D partial dependence surface"
#   )
#
# print(importance_plot)
# print(pdp_range_plot)
# print(pdp_restricted_share_plot)
# print(pdp_sigma_O_plot)
# print(pdp_encountered_restricted_plot)
# print(pdp_surface_plot)
#
# ggplot(diag_df, aes(encountered_rate_restricted_mean, power_signed, colour = mpa_effect_label)) +
#   geom_point(size = 2, alpha = 0.8) +
#   facet_wrap(~ eval_year) +
#   labs(
#     x = "Restricted encounter rate (species mean across surveys)",
#     y = "Signed power",
#     colour = "MPA effect"
#   )
#
# ggplot(diag_df, aes(phi_mean, power_signed, colour = mpa_effect_label)) +
#   geom_point(size = 2, alpha = 0.8) +
#   facet_grid(mpa_effect_label~ eval_year) +
#   labs(
#     x = "Phi",
#     y = "Signed power",
#     colour = "MPA effect"
#   )
#
#
# ggplot(diag_df, aes(sigma_O_mean, power_signed, colour = mpa_effect_label)) +
#   geom_point(size = 2, alpha = 0.8) +
#   facet_wrap(~ eval_year) +
#   labs(
#     x = "Spatial SD (species mean across surveys)",
#     y = "Signed power",
#     colour = "MPA effect"
#   )
#
# # This is the combination plot to focus on.
# ggplot(
#   diag_df,
#   aes(encountered_rate_restricted_mean, sigma_O_mean, colour = power_signed)
# ) +
#   geom_point(alpha = 0.9, size = 2.5) +
#   facet_grid(eval_year ~ mpa_effect_label) +
#   scale_colour_viridis_c() +
#   labs(
#     x = "Restricted encounter rate (species mean across surveys)",
#     y = "Spatial SD (species mean across surveys)",
#     colour = "Signed power"
#   )
#
# heatmap_df <- diag_df |>
#   mutate(
#     encounter_bin = ntile(encountered_rate_restricted_mean, 4),
#     sigma_bin = ntile(sigma_O_mean, 4)
#   ) |>
#   group_by(eval_year, mpa_effect_label, encounter_bin, sigma_bin) |>
#   summarise(
#     mean_power_signed = mean(power_signed, na.rm = TRUE),
#     n = n(),
#     .groups = "drop"
#   )
#
# ggplot(heatmap_df, aes(encounter_bin, sigma_bin, fill = mean_power_signed)) +
#   geom_tile() +
#   facet_grid(eval_year ~ mpa_effect_label) +
#   scale_fill_viridis_c() +
#   labs(
#     x = "Restricted encounter-rate quartile",
#     y = "Spatial SD quartile",
#     fill = "Mean signed power"
#   )
#
# heatmap_range_share <- make_binned_heatmap_data(
#   diag_df,
#   x_var = "range_mean",
#   y_var = "restricted_share_mean"
# )
#
# heatmap_range_encounter <- make_binned_heatmap_data(
#   diag_df,
#   x_var = "range_mean",
#   y_var = "encountered_rate_restricted_mean"
# )
#
# heatmap_sigma_share <- make_binned_heatmap_data(
#   diag_df,
#   x_var = "sigma_O_mean",
#   y_var = "restricted_share_mean"
# )
#
# print(
#   plot_binned_heatmap(
#     heatmap_range_share,
#     x_label = "Range quartile",
#     y_label = "Restricted-share quartile",
#     title = "Range × restricted share"
#   )
# )
#
# print(
#   plot_binned_heatmap(
#     heatmap_range_encounter,
#     x_label = "Range quartile",
#     y_label = "Restricted encounter-rate quartile",
#     title = "Range × restricted encounter rate"
#   )
# )
#
# print(
#   plot_binned_heatmap(
#     heatmap_sigma_share,
#     x_label = "Spatial SD quartile",
#     y_label = "Restricted-share quartile",
#     title = "Spatial SD × restricted share"
#   )
# )
#
