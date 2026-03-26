# Predicted distributions each species (last year sampled)
source(here::here("R", "00-fit-sim-functions.R"))
source(here::here("R", "00-utils.R"))

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(sf)

theme_set(gfplot::theme_pbs())


fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of species ultimately used in power analysis
spp_list <- list.files(here::here("data-generated", "power-results"))
spp_levels <- spp_list |> gsub("-", " ", x = _)

# Conditioning models for each species
fit_dir <- here::here("data-generated", "fits")
fit_files <- list.files(fit_dir, pattern = paste0(paste0(spp_list, collapse = "|"), ".*rds$"), full.names = TRUE)
sp_fits <- purrr::map(fit_files, readRDS)

# Spatial data for mapping and predictions
hbll_restricted <- readRDS(file.path("data-generated", "hbll-restricted-sf.rds")) |>
  sf::st_drop_geometry() |>
  dplyr::select(survey_abbrev, block_id, X, Y, restricted)

hbll_sf <- gfdata::load_survey_blocks(type = "polygon") |>
  filter(survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S", "HBLL INS N"))

# Predictions for each species
spp_preds <- sp_fits |>
  purrr::map_dfr(\(fit) {
    species <- fit$data$species_common_name |> unique()
    survey_abbrev <- fit$data$survey_abbrev |> unique()
    nd <- sdmTMB::replicate_df(
      dat = hbll_restricted,
      time_name = "year",
      time_values = max(fit$data$year)
    ) |>
    dplyr::mutate(fyear = as.factor(year))
    pred <- predict(fit, newdata = nd)
    dplyr::mutate(pred, species = species)
  })

# Get plot limits
mpa_500 <- readRDS(file.path("data-generated", "spatial", "simple-mpa-500m.rds")) |>
  st_transform(crs = st_crs(hbll_sf))
ylim <- c(st_bbox(mpa_500 |> rotate_sf())['ymin'], st_bbox(hbll_sf |> rotate_sf())['ymax'])


# spp_levels <- c(
#   "canary rockfish", "silvergray rockfish", "north pacific spiny dogfish",
#   "quillback rockfish", "lingcod", "yelloweye rockfish")
# Easier to use sf object for mapping
preds_sf <- spp_preds |>
  left_join(hbll_sf, by = c("survey_abbrev", "block_id")) |>
  st_as_sf() |>
  mutate(species = gsub("north ", "", species))
  # mutate(species = factor(species, levels = spp_levels)) |> # change if we order by something like max power
  # mutate(species = forcats::fct_recode(species, "pacific spiny dogfish" = "north pacific spiny dogfish"))

# Species-specific colour scale breaks
species_breaks <- list(
  "canary rockfish" = c(0.1, 1, 3),
  "silvergray rockfish" = c(0.1, 1, 5, 10),
  "pacific spiny dogfish" = c(1, 5, 20, 70),
  "quillback rockfish" = c(0.1, 2, 5, 15),
  "lingcod" = c(0.1, 2, 5, 10),
  "yelloweye rockfish" = c(1, 5, 20, 50)
)

plot_dists <- function(p_spp_sf) {
  # Get species name
  species_name <- unique(p_spp_sf$species)

  # Get species-specific breaks
  breaks_custom <- species_breaks[[species_name]]

  # Create labels that match the breaks exactly (no unnecessary decimals)
  labels_custom <- as.character(breaks_custom)

  p_spp_sf |> rotate_sf() |>
    ggplot() +
    geom_sf(data = pacea::bc_coast |> rotate_sf(), fill = "grey90", linewidth = 0.1) +
    geom_sf(aes(fill = exp(est) * 100, colour = exp(est) * 100)) +
    scale_fill_viridis_c(trans = "fourth_root_power",
      breaks = breaks_custom,
      labels = labels_custom) +
    scale_colour_viridis_c(trans = "fourth_root_power",
      breaks = breaks_custom,
      labels = labels_custom) +
    gfplot::coord_sf_auto(sf_obj = hbll_sf |> rotate_sf(), ylim = ylim) +
    facet_wrap(~ species, ncol = 6,
      labeller = labeller(species = function(x) tools::toTitleCase(tolower(x)))) +
    theme(legend.position = "bottom",
      strip.text = element_text(size = 10),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.text = element_text(size = 9)) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = 5, title = NULL)) +
    labs(fill = NULL, colour = NULL)
}

# preds_sf |> filter(species == "lingcod") |>
#   group_by(species) |>
#   group_split()  |>
#   purrr::map(plot_dists)

# dev.new(width = 9, height = 5.5)
dist_plots <- preds_sf |>
  group_by(species) |>
  group_split()  |>
  purrr::map(plot_dists)

p_grid <- patchwork::wrap_plots(dist_plots, ncol = 6) &
  theme(plot.margin = margin(2, 2, 2, 2))

p_grid <- p_grid +
  patchwork::plot_annotation(
    caption = "Expected mean catch per 100 hooks",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 10, margin = margin(t = 5))
    )
  )

p_grid

ggsave(file.path(fig_dir, "predicted-distributions.png"), width = 9, height = 5.5)