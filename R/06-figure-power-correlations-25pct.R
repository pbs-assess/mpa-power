library(dplyr)
library(ggplot2)
library(tidyr)

source("R/05-make-power-df.R")

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
  filter(
    mpa_effect_label == "25%",
    !eval_year %in% c(2030, 2034)
  ) |>
  mutate(
    eval_year = factor(eval_year, levels = c(2038, 2042, 2046)),
    species = stringr::str_to_title(species)
  )

plot_df <- diag_df |>
  transmute(
    species,
    eval_year,
    power_signed,
    `Mean encountered rate\nin MPA network` = encountered_rate_restricted_mean,
    `Mean interannual\nvariability SD` = sigma_V_mean,
    `Spatiotemporal\nrandom field SD` = sigma_E_mean
  ) |>
  pivot_longer(
    cols = c(
      `Mean encountered rate\nin MPA network`,
      `Mean interannual\nvariability SD`,
      `Spatiotemporal\nrandom field SD`
    ),
    names_to = "metric_label",
    values_to = "metric_value"
  ) |>
  mutate(
    metric_label = factor(
      metric_label,
      levels = c(
        "Mean encountered rate\nin MPA network",
        "Mean interannual\nvariability SD",
        "Spatiotemporal\nrandom field SD"
      )
    )
  )

corr_df <- plot_df |>
  group_by(metric_label, eval_year) |>
  summarise(
    corr = if (sum(!is.na(metric_value) & !is.na(power_signed)) >= 2) {
      cor(metric_value, power_signed, use = "complete.obs")
    } else {
      NA_real_
    },
    label = ifelse(is.na(corr), "r = NA", sprintf("r = %.2f", corr)),
    .groups = "drop"
  )

g <- ggplot(plot_df, aes(metric_value, power_signed, colour = species)) +
  geom_point() +
  facet_grid(eval_year ~ metric_label, scales = "free") +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(n.breaks = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    formula = y ~ x,
    mapping = aes(group = NULL, colour = NULL),
    colour = "grey50",
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_text(
    data = corr_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.3,
    size = 3,
    colour = "grey50"
  ) +
  labs(
    x = "Conditioning model property value",
    y = "Power (correctly signed)",
    colour = "Species"
  ) +
  gfplot::theme_pbs() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

print(g)

ggsave(
  "figures/power-correlations-25pct.png",
  g,
  width = 5,
  height = 5
)
