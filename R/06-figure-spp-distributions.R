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

# Easier to use sf object for mapping
preds_sf <- spp_preds |>
  left_join(hbll_sf, by = c("survey_abbrev", "block_id")) |>
  st_as_sf()

dev.new(width = 9, height = 5.)
preds_sf |> rotate_sf() |>
  mutate(species = factor(species, levels = spp_levels)) |> # change if we order by something like max power
  mutate(species = forcats::fct_recode(species, "pacific spiny dogfish" = "north pacific spiny dogfish")) |>
  ggplot() +
  geom_sf(data = pacea::bc_coast |> rotate_sf(), fill = "grey90", linewidth = 0.1) +
  geom_sf(aes(fill = exp(est) * 100, colour = exp(est) * 100)) +
  scale_fill_viridis_c(trans = "fourth_root_power",
    breaks = c(1, 5, 20, 80)) +
  scale_colour_viridis_c(trans = "fourth_root_power",
    breaks = c(1, 5, 20, 80)) +
  gfplot::coord_sf_auto(sf_obj = hbll_sf |> rotate_sf(), ylim = ylim) +
  facet_wrap(~ species, ncol = 6,
    labeller = labeller(species = function(x) tools::toTitleCase(tolower(x)))) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8)) +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.vjust = 0.8        # 0.5 centres it vertically alongside the bar
  )) +
  labs(
    fill = "Mean expected count per 100 hooks",
    colour = "Mean expected count per 100 hooks"
  )

ggsave(file.path(fig_dir, "predicted-distributions.png"), width = 9, height = 5.5)
