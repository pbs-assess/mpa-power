library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(mgcv)
library(rpart)
library(ranger)
theme_set(gfplot::theme_pbs())

source("R/05-make-power-df.R")

power_df <- readRDS("data-generated/power-results-df.rds")
theta <- readRDS(fit_characteristics_file)
ar1_theta <- readRDS(ar1_parameters_file) |>
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

# temp <- diag_df[,c("encountered_rate_restricted_mean", "encountered_rate_avg_mean", "encountered_count_per_year_restricted_mean", "sigma_O_mean", "sigma_E_mean", "phi_mean", "range_mean", "sigma_V_mean")] |>
#   distinct()
# GGally::ggpairs(temp)
# ggsave("figures/pairs-covariates-power.png", width = 10, height = 10)

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
make_corr_plot(encountered_rate_restricted_mean, "Mean encountered rate in MPA network")
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

# modelling:

d <- diag_df
# glimpse(d)
#
# range(d$power_signed)
# library(glmmTMB)
#
#
# fit <- glmmTMB(power_signed ~ scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) + scale(sigma_O_mean) + (1|mpa_effect_label) + (1|eval_year), family = ordbeta(), data = d)
# broom.mixed::tidy(fit, conf.int = TRUE) |> as.data.frame()
#
# fit_rs <- glmmTMB(power_signed ~ scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) + (scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) | mpa_effect_label) + (scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) | eval_year), family = gaussian(), data = d)
#
# summary(fit_rs)
# broom.mixed::tidy(fit_rs, conf.int = TRUE) |> as.data.frame()
#
# dd <- filter(d, mpa_effect_label)
# fit_rs <- glmmTMB(power_signed ~ scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) + (scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) | mpa_effect_label) + (scale(encountered_count_per_year_restricted_mean) + scale(sigma_V_mean) + scale(phi_mean) + scale(sigma_E_mean) | eval_year), family = gaussian(), data = d)
#
# summary(fit_rs)
# broom.mixed::tidy(fit_rs, conf.int = TRUE) |> as.data.frame()


