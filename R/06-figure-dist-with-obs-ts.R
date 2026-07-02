source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-power-analysis-functions.R"))

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(patchwork)

# Made in R/06-figure-spp-distributions.R
dist_plot <- readRDS(file.path("data-generated", "fig-ye-dist.rds"))

hist_path <- cleaned_data_dir
sampling_summary <- readRDS(file.path(sample_dir, "sampling-summary.rds"))

survey_abbrev <- "HBLL OUT N"
plan <- "historical survey-year bootstrap"
ar1_scenario <- "fitted_AR1"
time_scenario <- "twenty-five_years"
eval_years <- c(2030, 2034, 2038, 2042, 2046)

species_to_plot <- sampling_summary |>
  filter(
    species == "yelloweye rockfish",
    survey_abbrev == .env$survey_abbrev,
    plan == .env$plan,
    ar1_scenario == .env$ar1_scenario,
    time_scenario == .env$time_scenario
    # !grepl("halibut", species, ignore.case = TRUE)
  ) |>
  distinct(species) |>
  arrange(species) |>
  pull(species)

# Start with a basic multi-species example grid. Expand trend/replicate next.
example_grid <- tidyr::expand_grid(
  species = species_to_plot,
  mpa_trend = 1.009,
  replicate = 1:3
)

hist_cache <- new.env(parent = emptyenv())

split_species_label <- function(x) {
  dplyr::if_else(
    x == "Pacific Spiny Dogfish",
    "Pacific Spiny\nDogfish",
    stringr::str_replace(x, "^(\\S+)\\s+", "\\1\n")
  )
}

format_species_label <- function(x) {
  x |>
    stringr::str_to_title() |>
    dplyr::recode("North Pacific Spiny Dogfish" = "Pacific Spiny Dogfish") |>
    split_species_label()
}

species_labels <- format_species_label(species_to_plot)

trend_to_total_increase <- function(mpa_trend, n_years = 25) {
  mpa_trend^n_years - 1
}

format_panel_label <- function(mpa_trend, replicate, n_years = 25) {
  total_increase <- trend_to_total_increase(mpa_trend, n_years = n_years) * 100
  paste0(
    "Rep ", replicate,
    " | ", round(total_increase), "% increase"
  )
}

panel_levels <- example_grid |>
  distinct(mpa_trend, replicate) |>
  mutate(panel = format_panel_label(mpa_trend, replicate)) |>
  pull(panel)

build_obs_ts_plot_data <- function(species,
                                   survey_abbrev,
                                   mpa_trend,
                                   replicate,
                                   ar1_scenario = "fitted_AR1",
                                   time_scenario = "twenty-five_years",
                                   plan = "historical survey-year bootstrap",
                                   eval_years = c(2030, 2034, 2038, 2042, 2046),
                                   hist_path = here::here("data-generated", "cleaned-species-data"),
                                   sample_dir = here::here("data-generated", "sampled-data"),
                                   sampling_summary = readRDS(file.path(sample_dir, "sampling-summary.rds")),
                                   hist_cache = new.env(parent = emptyenv())) {

  sampled_data <- load_sampled_data(
    species = species,
    survey_abbrev = survey_abbrev,
    mpa_trend = mpa_trend,
    ar1_scenario = ar1_scenario,
    time_scenario = time_scenario,
    plan = plan,
    replicates = replicate,
    sampling_summary = sampling_summary,
    sample_dir = sample_dir
  )

  hist_data <- get_hist_data(
    species = species,
    survey_abbrev = survey_abbrev,
    hist_path = hist_path,
    cache_env = hist_cache
  )

  combine_hist_sim_data(
    sim_data = sampled_data,
    hist_data = hist_data,
    eval_year = max(eval_years)
  ) |>
    group_by(year, restricted, historical) |>
    summarise(
      mean_cp = mean(catch_prop, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(
      restricted_f = factor(restricted, levels = c(0, 1), labels = c("Outside MPA", "Inside MPA")),
      species_f = factor(format_species_label(species), levels = species_labels),
      replicate = replicate,
      mpa_trend = mpa_trend,
      panel = factor(format_panel_label(mpa_trend, replicate), levels = panel_levels)
    )
}

obs_ts_plot_data <- purrr::pmap_dfr(
  example_grid,
  \(species, mpa_trend, replicate) {
    build_obs_ts_plot_data(
      species = species,
      survey_abbrev = survey_abbrev,
      mpa_trend = mpa_trend,
      replicate = replicate,
      ar1_scenario = ar1_scenario,
      time_scenario = time_scenario,
      plan = plan,
      eval_years = eval_years,
      hist_path = hist_path,
      sample_dir = sample_dir,
      sampling_summary = sampling_summary,
      hist_cache = hist_cache
    )
  }
) |>
  mutate(rep_panel = paste("Rep", replicate),
         strip_text_top = "25% increase")

hist_breaks <- obs_ts_plot_data |>
  filter(historical) |>
  group_by(species_f, panel) |>
  summarise(last_hist_year = max(year), .groups = "drop")

ts_plot <- ggplot(
  obs_ts_plot_data,
  aes(year, mean_cp * 100, colour = restricted_f, linetype = historical)
) +
  geom_path() +
  geom_point(size = 1.2) +
  geom_vline(
    data = hist_breaks,
    aes(xintercept = last_hist_year),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = eval_years,
    linetype = 2,
    alpha = 0.2
  ) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "solid"), guide = "none") +
  scale_y_log10() +
  facet_grid(rep_panel ~ strip_text_top, scales = "free_y") +
  labs(
    # title = paste("Observed TS comparisons |", survey_abbrev),
    x = NULL,
    y = "Mean observed count per 100 hooks",
    colour = NULL
  ) +
  gfplot::theme_pbs() +
  scale_colour_manual(
    values = c(
      "Outside MPA" = "grey60",
      "Inside MPA" = RColorBrewer::brewer.pal(5, "Set2")[2]
    )
  ) +
  theme(legend.position = "top")

(dist_plot +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position = "top",
        legend.margin = margin(t = 0, r = 0, b = -3, l = 0)
        ) +
  guides(fill = "none", colour = "none") +
  guides(fill = guide_colorbar(title.position = "top",
    barheight = unit(0.02, "npc"),
    barwidth = unit(0.28, "npc"),
    title = "Expected catch per 100 hooks"))) +
ts_plot +
  theme(strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(r = 8))) +
  guides(colour = guide_legend(
      override.aes = list(
        size = 1.3,        # bigger points
        linewidth = 0.6    # thicker lines
      )
    )) +
    theme(
      legend.margin = margin(t = -3, r = 0, b = 0, l = 0),
      legend.spacing.x = unit(1, "cm"),  # more space between legend items
      legend.text.position = "top",
      legend.text = element_text(margin = margin(l = 5, r = 5, b = 8), size = 9)  # margin around text
    ) +
plot_annotation(tag_levels = "a", tag_prefix = "", tag_suffix = ")") &
theme(plot.tag.position = c(0, 0.99), plot.tag = element_text(size = 10))

ggsave("figures/ye-dist-with-obs-ts.png", width = 5.9, height = 6)


ts_plot +
  (obs_ts_plot_data |>
  filter(replicate %in% c(1, 3)) |>
  mutate(rep_panel = gsub("3", "2", rep_panel))) +
  theme(axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        axis.text = element_text(size = 14),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        legend.margin = margin(t = -3, r = 0, b = 20, l = 0),
        legend.text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list( size = 2, linewidth = 0.9)))
ggsave("figures/presentations/2026-05-05-CSAS-meeting/ts-plot.png", width = 9, height = 5.3)
