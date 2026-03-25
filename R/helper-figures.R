# Figures to help describe the dimensions considered in the power analysis

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

library(ggsidekick)
library(tidyr)
library(patchwork)
library(dplyr)

sample_dir <- here::here("data-generated", "sampled-data")

presentation <- FALSE

if (presentation) {
  angle <- 0
  fig_dir <- here::here("figures", "presentations", "2026-03-10-survey-design-workshop")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  auto_coord <- function() gfplot::coord_sf_auto(hbll_grid_sf |> rotate_sf(angle = angle), buffer = 0)
} else {
  angle <- -40
  fig_dir <- here::here("figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  auto_coord <- function() gfplot::coord_sf_auto(display_mpa |> rotate_sf(angle = angle), buffer = 0)
}

# Allocation values for HBLL survey regions
hbll_region_levels <- c("HBLL OUT N", "HBLL OUT S", "HBLL INS N")
hbll_allocations <- readRDS(here::here("data-generated", "hbll-allocations.rds"))
hbll_allocations |>
  group_by(survey_abbrev) |>
  summarise(allocation = sum(allocation)) |>
  ungroup() |>
  mutate(survey_abbrev = factor(survey_abbrev, levels = hbll_region_levels)) |>
  filter(survey_abbrev != "HBLL INS S") |>
  mutate(text = paste0(allocation, " (", survey_abbrev, ")")) |>
  arrange(survey_abbrev) |>
  pull(text) |>
  paste(collapse = ", ")

hbll_grid <- gfdata::load_survey_blocks(type = "XY") |>
  filter(survey_abbrev %in% hbll_region_levels)
hbll_grid_sf <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(survey_abbrev %in% hbll_region_levels)
hbll_grid_crs <- st_crs(hbll_grid_sf)
hbll_outline <- hbll_grid_sf |> st_union()
hbll_region_outline <- hbll_grid_sf |> group_by(survey_abbrev) |> st_union()

fit_dir <- here::here("data-generated", "fits")

spp <- c("lingcod", "north pacific spiny dogfish", "yelloweye rockfish", "quillback rockfish", "pacific halibut")
survey_cols <- c("HBLL OUT N" = "#0072B2", "HBLL OUT S" = "#D55E00", "HBLL INS N" = "#009E73")
display_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa-500m.rds"))

historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

historical_n_sets <- historical_locations |>
  group_by(restricted) |>
  st_drop_geometry() |>
  summarise(n = n())

# Label positions: centroid of each survey region (so labels stay inside plot)
hbll_labels <- hbll_grid_sf |>
  group_by(survey_abbrev) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid() |>
  mutate(label = as.character(survey_abbrev))

hbll_labels |> st_transform(crs = 4326)

if (presentation) {
  hbll_labels <- tribble(
    ~X, ~Y, ~survey_abbrev,
    -130, 53.9, "HBLL OUT N",
    -130, 50, "HBLL OUT S",
    -126.7, 51.2, "HBLL INS N"
  ) |>
    mutate(survey_abbrev = factor(survey_abbrev, levels = hbll_region_levels)) |>
    st_as_sf(coords = c("X", "Y"), crs = 4326)
} else {
  hbll_labels <- tribble(
    ~X, ~Y, ~survey_abbrev,
    -129.9, 53.6, "HBLL OUT N",
    -127.3, 51.2, "HBLL OUT S",
    -126, 51, "HBLL INS N"
  ) |>
    mutate(survey_abbrev = factor(survey_abbrev, levels = hbll_region_levels)) |>
    st_as_sf(coords = c("X", "Y"), crs = 4326)
}

# st_bbox(display_mpa |> st_transform(crs = 4326))
#       xmin       ymin       xmax       ymax
# -134.08115   50.09211 -124.78667   55.41001

# Overall map of survey regions
ne_coast <- rnaturalearth::ne_states(country = c("canada", "united states of america"), returnclass = "sf") |>
  st_transform(crs = 3005) |>
  filter(name %in% c("British Columbia", "Washington", "Alaska"))

p1 <- ggplot() +
  geom_sf(data = ne_coast |> rotate_sf(angle = angle), fill = NA, linewidth = 0.1) +
  geom_sf(data = display_mpa |> rotate_sf(angle = angle), fill = "grey85", linewidth = 0.15) +
  geom_sf(data = hbll_grid_sf |> rotate_sf(angle = angle),
    aes(fill = survey_abbrev, colour = survey_abbrev),
    alpha = 0.5, linewidth = 0) +
  scale_fill_manual(values = survey_cols) +
  scale_colour_manual(values = survey_cols) +
  guides(fill = "none", colour = "none") +
  auto_coord() +
  theme(axis.title = element_blank())

p_out <- p1 +
  geom_sf(data = historical_locations |> filter(restricted == 0) |> rotate_sf(angle = angle),
    shape = 21, size = 0.4) +
  auto_coord() +
  theme(axis.title = element_blank()) +
  ggtitle(paste0("a) Survey sets outside MPA (n = ", historical_n_sets$n[1], ")"))
p_in <- p1 +
  geom_sf(data = historical_locations |> filter(restricted == 1) |> rotate_sf(angle = angle),
    shape = 21, size = 0.4) +
  geom_sf_label(data = hbll_labels |> rotate_sf(angle = angle),
    aes(label = survey_abbrev, fill = survey_abbrev),
    size = ifelse(presentation, 5, 3),
    hjust = 0, alpha = 0.7, colour = "black", linewidth = 0) +
  auto_coord() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle(paste0("b) Survey sets inside MPA (n = ", historical_n_sets$n[2], ")"))

p_out + p_in
if (presentation) {
  ggsave(file.path(fig_dir, "hbll-mpa-historical-overlay.png"), width = 13.4, height = 7.8)
} else {
  ggsave(file.path(fig_dir, "hbll-mpa-historical-overlay.png"), width = 7.5, height = 7.8)
}

# (1) Plot predicted distributions each species (last year sampled)
sp_fits <- purrr::map(spp, load_cached_species, fit_dir = fit_dir)

hbll_restricted <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  st_drop_geometry() |>
  select(survey_abbrev, block_id, X, Y, restricted)

spp_preds <- sp_fits |>
  purrr::flatten() |>
  purrr::imap(\(fit, name) {
    species <- fit$data$species_common_name |> unique()
    pred <- predict_hbll(fit, hbll_restricted)
    mutate(pred, species = species)
  })

preds_sf <- spp_preds |>
  bind_rows() |>
  left_join(hbll_grid_sf, y = _)

preds_sf |> rotate_sf(angle = angle) |>
  group_by(survey_abbrev)  |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(species = factor(species, levels = spp)) |>
  mutate(species = forcats::fct_recode(species, "pacific spiny dogfish" = "north pacific spiny dogfish")) |>
  ggplot() +
  geom_sf(data = pacea::bc_coast |> rotate_sf(angle = angle), fill = "grey90", linewidth = 0.1) +
  geom_sf(aes(fill = exp(est) * 100, colour = exp(est) * 100)) +
  scale_fill_viridis_c(trans = "fourth_root_power",
    breaks = c(1, 5, 20, 80)) +
  scale_colour_viridis_c(trans = "fourth_root_power",
    breaks = c(1, 5, 20, 80)) +
  gfplot::coord_sf_auto(hbll_grid_sf |> rotate_sf(angle = angle), buffer = 0) +
  facet_wrap(~ species, ncol = ifelse(presentation, 5, 5),
    labeller = labeller(species = function(x) tools::toTitleCase(tolower(x)))) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = ifelse(presentation, 20, 8))) +
  labs(
    fill = "Mean expected count\nper 100 hooks",
    colour = "Mean expected count\nper 100 hooks"
  )
if (presentation) {
ggsave(file.path(fig_dir, "predicted-distributions.png"), width = 19, height = 6)
} else {
ggsave(file.path(fig_dir, "predicted-distributions.png"), width = 7, height = 5.7)
}

#Predicted mean catch per 100 hooks for five groundfish species from spatial GLMMs
# fitted to Hard Bottom Longline (HBLL) survey data
# (Outside North, Outside South, Inside North).
# Predictions shown for the most recent year sampled in each survey region;
# British Columbia coast shown for reference."

# (2) Plot the underlying static spatial field for each species
# This is one component of the simulation
preds_sf |>
group_by(survey_abbrev)  |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(species = factor(species, levels = spp)) |>
  mutate(species = forcats::fct_recode(species, "pacific spiny dogfish" = "north pacific spiny dogfish")) |>
  ggplot() +
  geom_sf(data = pacea::bc_coast |> rotate_sf(angle = angle), fill = "grey90", linewidth = 0.1) +
  geom_sf(aes(fill = omega_s, colour = omega_s)) +
  geom_sf(data = hbll_outline, fill = "transparent", linewidth = 0.15) +
  scale_fill_gradient2() +
  scale_colour_gradient2() +
  facet_wrap(~ species, ncol = 5) +
  theme(legend.position = "bottom") +
  auto_coord()

if (presentation) {
  preds_sf |>
  filter(species == "yelloweye rockfish") |>
group_by(survey_abbrev)  |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(species = factor(species, levels = spp)) |>
  mutate(species = forcats::fct_recode(species, "pacific spiny dogfish" = "north pacific spiny dogfish")) |>
  ggplot() +
  geom_sf(data = pacea::bc_coast |> rotate_sf(angle = angle), fill = "grey90", linewidth = 0.1) +
  geom_sf(aes(fill = omega_s, colour = omega_s)) +
  geom_sf(data = hbll_outline, fill = "transparent", linewidth = 0.15) +
  scale_fill_gradient2() +
  scale_colour_gradient2() +
  theme(legend.position = "bottom") +
  auto_coord() +
  ggtitle("Yelloweye Rockfish")
  ggsave(file.path(fig_dir, "spatial-random-field-yelloweye.png"), width = 3, height = 4)
} else {
  ggsave(file.path(fig_dir, "spatial-random-field.png"), width = 7, height = 5.7)
}

# For each simulation replicate, spatial random effects were drawn from a single
# realisation from their conditional posterior given the conditioning model fixed effects.
# This spatial field is constant over time within each replicate, but varies
# between replicates to reflect some of the uncertainty associated with this spatial random field.


# (3) Sampling scenario comparison
sample_f <- list.files(file.path(sample_dir, "yelloweye-rockfish"),
  pattern = ".*mpa1\\.018.*", full.names = TRUE)

sampled_data <- purrr::map_dfr(sample_f, readRDS)
plan_levels <- c("status quo", "MPAs at 5 year intervals")

# sampled_data |>
#   mutate(plan = factor(plan, levels = plan_levels)) |>
#   mutate(plan = forcats::fct_recode(plan, "Status quo" = "status quo", "MPAs at 4 year intervals" = "MPAs at 5 year intervals")) |>
#   filter(replicate == 1) |>
#   filter(year <= 2028) |>
#   XY_to_sf(crs_to = st_crs(hbll_grid_crs)) |>
#   rotate_sf(angle = angle) |>
# ggplot() +
#   geom_sf(data = ne_coast |> rotate_sf(angle = angle), fill = "white", linewidth = 0.05) +
#   geom_sf(data = display_mpa |> rotate_sf(angle = angle), fill = "grey85", linewidth = 0.05) +
#   geom_sf(aes(colour = observed, shape = factor(restricted)), size = 1.2, alpha = 0.5) +
#   scale_shape_manual(name = "Restricted", values = c(`0` = 21, `1` = 19)) +
#   scale_colour_viridis_c(name = "Observed", breaks = scales::breaks_pretty(n = 3)) +
#   theme(legend.position = "bottom") +
#   theme(strip.text = element_text(size = ifelse(presentation, 18, 8))) +
#   facet_grid(plan ~ year) +
#   auto_coord()
# if (presentation) {
#   ggsave(file.path(fig_dir, "sampling-scenario-comparison.png"), width = 11.7, height = 7.1)
# } else {
#   ggsave(file.path(fig_dir, "sampling-scenario-comparison.png"), width = 9, height = 10.9)
# }

 # (2) From a replicate of the simulated data, plot the relative index of abundance
# for each species coastwide

# ------------------------------------------------------------------------------
# Compute relative index for a single survey region + replicate.
#
# Extracted from the one-off exploratory code below so you can adjust
# `survey_abbrev` (region) and `replicate_id` without re-editing the pipeline.
relative_index_single_survey <- function(
  species,
  survey_abbrev,
  replicate_id,
  mpa_trend = 1.016,
  ar1_scenario = "fitted_AR1",
  time_scenario = "twenty-five_years",
  plan = "historical survey-year bootstrap",
  nsim_predict = 100,
  predict_seed = 123
) {
  fit_survey_lu <- tribble(
    ~survey_abbrev, ~fit_name,
    "HBLL OUT N", "fit_ON",
    "HBLL OUT S", "fit_OS",
    "HBLL INS N", "fit_IN"
  )

  if (!survey_abbrev %in% fit_survey_lu$survey_abbrev) {
    stop(
      "Unknown survey_abbrev: ", survey_abbrev,
      "\nExpected one of: ", paste(fit_survey_lu$survey_abbrev, collapse = ", ")
    )
  }

  # Load conditioning fit for the requested survey region.
  fits_single <- load_cached_species(species, fit_dir = fit_dir)
  fit_single <- fits_single[[fit_survey_lu$fit_name[fit_survey_lu$survey_abbrev == survey_abbrev]]]

  canonical_plan <- dplyr::case_match(
    plan,
    "MPAs at 5 year intervals" ~ "MPAs every 4 years",
    "MPAs at 4 year intervals" ~ "MPAs every 4 years",
    .default = plan
  )

  sampling_summary_path <- file.path(sample_dir, "sampling-summary.rds")
  if (!file.exists(sampling_summary_path)) {
    stop("Sampling summary not found: ", sampling_summary_path)
  }

  sampling_summary <- readRDS(sampling_summary_path)

  selected_file <- sampling_summary |>
    filter(
      species == .env$species,
      survey_abbrev == .env$survey_abbrev,
      mpa_trend == round(.env$mpa_trend, digits = 3),
      ar1_scenario == .env$ar1_scenario,
      time_scenario == .env$time_scenario,
      plan == .env$canonical_plan,
      replicate == .env$replicate_id
    )

  if (nrow(selected_file) == 0) {
    available_scenarios <- sampling_summary |>
      filter(
        species == .env$species,
        survey_abbrev == .env$survey_abbrev,
        ar1_scenario == .env$ar1_scenario,
        time_scenario == .env$time_scenario
      ) |>
      distinct(plan, mpa_trend) |>
      arrange(plan, mpa_trend)

    stop(
      "No sampled file found for species=", species,
      ", survey_abbrev=", survey_abbrev,
      ", mpa_trend=", round(mpa_trend, digits = 3),
      ", ar1_scenario=", ar1_scenario,
      ", time_scenario=", time_scenario,
      ", plan=", canonical_plan,
      ", replicate=", replicate_id,
      "\nAvailable plan/trend combinations:\n",
      paste0("  - ", available_scenarios$plan, " (", available_scenarios$mpa_trend, ")", collapse = "\n")
    )
  }

  sampled_path <- file.path(sample_dir, selected_file$file[[1]])
  if (!file.exists(sampled_path)) {
    stop("Sampled data not found: ", sampled_path)
  }

  sampled_dat_single <- readRDS(sampled_path)

  first_sim_year <- min(sampled_dat_single$year)

  # Prediction grid from the requested replicate.
  # Note: file contains only this replicate, so just filter by year
  pred_grid_single <- sampled_dat_single |>
    filter(year == first_sim_year) |>
    distinct(survey_abbrev, X, Y, block_id, restricted)

  hist_years_single <- sort(unique(fit_single$data$year))

  nd_single <- pred_grid_single |>
    expand_grid(year = hist_years_single) |>
    mutate(fyear = as.factor(year))

  # Historical predictions (from the conditioning model, using the grid defined above).
  set.seed(predict_seed)
  pred_matrix_single <- predict(
    fit_single,
    newdata = nd_single,
    nsim = nsim_predict,
    type = "response"
  )
  nd_single$mu <- apply(pred_matrix_single, 1, median)

  hist_index_single <- nd_single |>
    select(year, mu) |>
    group_by(year) |>
    summarise(index = sum(mu), .groups = "drop") |>
    mutate(period = "historical")

  # Simulated index for the requested replicate.
  # Note: file contains only this replicate, no need to filter
  sim_index_single <- sampled_dat_single |>
    group_by(year, restricted) |>
    summarise(index = sum(mu), .groups = "drop") |>
    mutate(period = "simulated")

  # === Scale by last historical year ===
  last_hist_year_single <- max(hist_index_single$year)

  # Baseline for each restricted group from last historical year
  # (from predictions, not from historical data, since historical wasn't split).
  baseline_single_by_restricted <- nd_single |>
    filter(year == last_hist_year_single) |>
    group_by(restricted) |>
    summarise(baseline = sum(mu), .groups = "drop")

  total_baseline_single <- sum(baseline_single_by_restricted$baseline)

  # Scale historical (single combined line)
  hist_scaled_single <- hist_index_single |>
    mutate(relative_index = index / total_baseline_single) |>
    select(year, relative_index, period) |>
    mutate(area = "Historical (combined)")

  # Scale simulated (split by restricted, each relative to its own baseline)
  sim_scaled_single <- sim_index_single |>
    left_join(baseline_single_by_restricted, by = "restricted") |>
    mutate(relative_index = index / baseline) |>
    select(year, restricted, relative_index, period) |>
    mutate(area = ifelse(restricted == 1, "Inside MPA", "Outside MPA"))

  # Add connection point: last historical year split by restricted
  connection_point <- baseline_single_by_restricted |>
    mutate(
      year = last_hist_year_single,
      relative_index = baseline / baseline,  # = 1.0
      period = "historical",
      area = ifelse(restricted == 1, "Inside MPA", "Outside MPA")
    ) |>
    select(year, restricted, relative_index, period, area)

  bind_rows(
    hist_scaled_single,
    connection_point,
    sim_scaled_single
  )
}

example_species <- "yelloweye rockfish"
example_survey <- "HBLL OUT N"

test <- relative_index_single_survey(
  species = example_species,
  survey_abbrev = example_survey,
  replicate_id = 4,
  mpa_trend = 1.009,
  ar1_scenario = "fitted_AR1",
  time_scenario = "twenty-five_years",
  plan = "historical survey-year bootstrap"
)

ggplot(data = test) +
  aes(x = year, y = relative_index, colour = area, group = area) +
  geom_line(linewidth = 1) +
  # geom_vline(xintercept = last_hist_year_single + 0.5, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
  geom_vline(xintercept = c(2030, 2034, 2038, 2042, 2046), linetype = "dashed", alpha = 0.2) +
  scale_colour_manual(
    values = c("Historical (combined)" = "grey30",
               "Inside MPA" = "#D55E00",
               "Outside MPA" = "#0072B2")
  ) +
  labs(title = paste(stringr::str_to_title(example_species), " - ", example_survey),
        y = "Relative index",
        x = "Year",
        colour = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    plot.margin = margin(t = 3, r = 10, b = 3, l = 3)
  )

ggsave(file.path(fig_dir, "relative-index-single-survey.png"), width = 7.2, height = 4.3)

# replicate_ids <- 1
# all_replicates_indices_outn <- purrr::map_dfr(replicate_ids, \(rep) {
#   res <- relative_index_single_survey(
#     species = "yelloweye rockfish",
#     survey_abbrev = "HBLL OUT N",
#     replicate_id = rep,
#     mpa_trend = 1.021,
#     ar1_scenario = "fitted_AR1",
#     time_scenario = "twenty-five_years",
#     plan = "status quo"
#   )
#   res |>
#     dplyr::mutate(rep = rep)
# })

# all_replicates_indices_outs <- purrr::map_dfr(replicate_ids, \(rep) {
#   res <- relative_index_single_survey(
#     species = "yelloweye rockfish",
#     survey_abbrev = "HBLL OUT S",
#     replicate_id = rep,
#     mpa_trend = 1.021,
#     ar1_scenario = "fitted_AR1",
#     time_scenario = "twenty-five_years",
#     plan = "status quo"
#   )
#   res |>
#     dplyr::mutate(rep = rep)
# })

# ggplot(data = all_replicates_indices |> filter(rep == 1)) +
#   aes(x = year, y = relative_index, colour = area, group = interaction(rep, area)) +
#   # geom_line(data = all_replicates_indices |> filter(period == "historical"), linewidth = 1) +
#   # geom_line(data = all_replicates_indices |> filter(period == "simulated"), linewidth = 1, alpha = 0.5) +
#   geom_line(linewidth = 1) +
#   geom_vline(xintercept = last_hist_year_single + 0.5, linetype = "dashed", alpha = 0.5) +
#   geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
#   geom_vline(xintercept = c(2030, 2034, 2038, 2042, 2046), linetype = "dashed", alpha = 0.2) +
#   scale_colour_manual(
#     values = c("Historical (combined)" = "grey30",
#                "Inside MPA" = "#D55E00",
#                "Outside MPA" = "#0072B2")
#   ) +
#   labs(title = paste(stringr::str_to_title(sp), " - ", survey_single),
#         y = "Relative index",
#         x = "Year",
#         colour = NULL) +
#   theme(
#     legend.position = "inside",
#     legend.position.inside = c(0.15, 0.85),
#     plot.margin = margin(t = 3, r = 10, b = 3, l = 3)
#     )


example_survey <- "HBLL OUT S"
all_indices_single <- relative_index_single_survey(
  species = example_species,
  survey_abbrev = example_survey,
  replicate_id = 6,
  mpa_trend = 1.009,
  ar1_scenario = "fitted_AR1",
  time_scenario = "twenty-five_years",
  plan = "historical survey-year bootstrap"
)

ggplot(all_indices_single, aes(x = year, y = relative_index, colour = area, group = area)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
  geom_vline(xintercept = c(2030, 2034, 2038, 2042, 2046), linetype = "dashed", alpha = 0.2) +
  scale_colour_manual(
    values = c("Historical (combined)" = "grey30",
               "Inside MPA" = "#D55E00",
               "Outside MPA" = "#0072B2")
  ) +
  labs(title = paste(stringr::str_to_title(example_species), " - ", example_survey),
        y = "Relative index",
        x = "Year",
        colour = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    plot.margin = margin(t = 3, r = 10, b = 3, l = 3)
    )

if (presentation) {
  ggsave(file.path(fig_dir, "relative-index-single-survey.png"), width = 7, height = 4)
} else {
  ggsave(file.path(fig_dir, "relative-index-single-survey.png"), width = 7, height = 5.7)
}
