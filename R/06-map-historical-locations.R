
source(here::here("R", "00-setup.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# library(ggsidekick)
library(tidyr)
library(patchwork)
library(dplyr)

sample_dir <- here::here("data-generated", "sampled-data")

angle <- -40
fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

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

sp <- "yelloweye rockfish"
survey_cols <- c("HBLL OUT N" = "#0072B2", "HBLL OUT S" = "#D55E00", "HBLL INS N" = "#009E73")
display_mpa <- readRDS(file.path("data-generated", "spatial", "simple-mpa-500m.rds"))

ye_files <- list.files(file.path("data-generated", "cleaned-species-data"), pattern = "yelloweye-rockfish.*", full.names = TRUE)
ye_data <- purrr::map_dfr(ye_files, readRDS)


# historical_locations <- readRDS(file.path("data-generated", "historical-locations.rds")) |>
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# historical_n_sets <- historical_locations |>
#   group_by(restricted) |>
#   st_drop_geometry() |>
#   summarise(n = n())

historical_n_sets <- ye_data |>
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