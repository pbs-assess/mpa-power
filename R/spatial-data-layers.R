# ---
# title: "Rmarkdown HTML including Leaflet"
# output: html_document
# ---

# ```{r}
# knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

library(flextable)
library(tidyverse)
library(plotly)
library(sf)
library(patchwork)

source(here::here("R", "00-utils.R"))
# source(here::here("R", "00-fit-sim-functions.R"))

# Plot settings
theme_set(gfplot::theme_pbs())
# plot_limits <- coord_sf(xlim = c(-134, -124), ylim = c(48.5, 54.5), crs = 4326)
plot_limits <- coord_sf(xlim = c(-134, -124), ylim = c(50, 54.5), crs = 4326)

# Prepare spatial data layers
cache_dir <- "~/R_DFO/gfsynopsis-2024-data/report/data-cache-2025-03"

# HBLL survey blocks
sb <- load_survey_blocks(type = "polygon")
hbll_blocks <- sb |>
  filter(str_detect(survey_abbrev, "HBLL"))
syn_blocks <- sb |>
  filter(str_detect(survey_abbrev, "SYN"))

survey_locations_sf <- readRDS(file.path(cache_dir, "yelloweye-rockfish.rds"))$survey_sets |>
  distinct(fishing_event_id, latitude, longitude) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(st_crs(hbll_blocks))

sampled_blocks <- st_join(hbll_blocks, survey_locations_sf) |>
  filter(!is.na(fishing_event_id))
unsampled_blocks <- st_join(hbll_blocks, survey_locations_sf) |>
  filter(is.na(fishing_event_id))

# MPAs
shape <- st_read(here::here("data-raw", "spatial", "Spatial_Q.gdb"),
  layer = "Q1_FULL_March2023") |>
  janitor::clean_names()
# Commercial longline
comm_ll_allowance <- readRDS(here::here("data-generated", "spatial", "comm_ll_allowance.rds"))

comm_ll_allowance_simplified <- comm_ll_allowance |>
  st_simplify(dTolerance = 200)

# old layer
old_hu_ll <- st_read(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"),
  layer = "hu_co_demersalfishing_bottomlongline_d")

# Newest layer from Kevin
hu_ll <- st_read(here::here("data-raw", "spatial", "hu.gdb"),
  layer = "hu_co_demersalfishing_bottomlongline_d")

# Yelloweye commercial longline CPUE
ye_dat <- readRDS(file.path(cache_dir, "yelloweye-rockfish.rds"))
ye_ll_cpue <- ye_dat$cpue_spatial_ll
ye_trawl_cpue <- ye_dat$cpue_spatial

#
ye_ll_sf <- st_as_sf(ye_ll_cpue, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(st_crs(hbll_blocks))

ye_trawl_sf <- st_as_sf(ye_trawl_cpue, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(st_crs(hbll_blocks))

# Make grid that matches the hbll blocks
summarise_cpue_by_grid <- function(cpue_sf, grid_sf, source) {
  ref_bbox <- st_bbox(grid_sf)
  cell_size <- 2000 # Width of one cell

  matching_grid_sf <- st_make_grid(
    x = cpue_sf,      # Use points to define the full extent
    cellsize = cell_size,   # Use the same cell size as the reference grid
    offset = c(ref_bbox$xmin, ref_bbox$ymin), # Use the same origin
    square = TRUE           # Ensure it's a square grid
  ) %>%
    st_sf(grid_id = 1:length(.), geometry = .)

  ll_in_grid <- st_join(cpue_sf, matching_grid_sf)

  # Group by grid cell, summarize data, and add the privacy flag
  aggregated_grid_sf <- ll_in_grid |>
    filter(!is.na(cpue)) |>
    st_drop_geometry() |> # faster to calculate without geometry
    group_by(grid_id) |>
    summarise(
      n_points = n(),
      total_cpue = sum(cpue, na.rm = TRUE),
      geo_mean_cpue = exp(mean(log(cpue), na.rm = FALSE))  # geometric mean
    ) |>
    ungroup() |>
    mutate(is_private = n_points < 3,
      source = source) |>
    right_join(matching_grid_sf, by = "grid_id") |>
    filter(!is.na(geo_mean_cpue)) |>
    st_as_sf()
}

ye_ll_cpue_grid <- summarise_cpue_by_grid(ye_ll_sf, hbll_blocks, "longline")
ye_trawl_cpue_grid <- summarise_cpue_by_grid(ye_trawl_sf, syn_blocks, "trawl")

ggplot() +
  geom_sf(data = comm_ll_allowance, aes(fill = activity_allowance_label), alpha = 0.7) +
  geom_sf(data = pacea::bc_coast, fill = "grey99") +
  scale_fill_manual(name = "activity_allowance_label", values = c("allowed" = "#0072B2",
    "not allowed" = "#D55E00", "conditional" = "#F0E442",
    "not applicable" = "#999999"), na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = ye_ll_cpue_grid |> filter(is_private), aes(fill = geo_mean_cpue), colour = NA) +
  scale_fill_viridis_c(trans = "log10", na.value = "grey90") +
  geom_sf(data = hu_ll, fill = "grey90", colour = "grey80", alpha = 0.5) +
  geom_sf(data = hbll_blocks, fill = "NA", colour = "red", linewidth = 0.05) +
  geom_sf(data = sampled_blocks, colour = "black", fill = NA, linewidth = 0.05) +
  geom_sf(data = comm_ll_allowance, aes(colour = activity_allowance_label), fill = NA) +
  scale_colour_manual(name = "activity_allowance_label", values = c("allowed" = "#0072B2",
    "not allowed" = "#D55E00", "conditional" = "#F0E442",
    "not applicable" = "#999999"), na.value = "grey90") +
  plot_limits
ggsave(here::here("draft-figures", "static-plot3.png"), width = 19, height = 11.5)

ggplot() +
  geom_sf(data = comm_ll_allowance, aes(fill = activity_allowance_label), colour = NA,alpha = 0.3) +
  scale_fill_manual(name = "activity_allowance_label", values = c("allowed" = "#0072B2",
    "not allowed" = "#D55E00", "conditional" = "#F0E442",
    "not applicable" = "#999999"), na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = hbll_blocks, aes(fill = survey_abbrev), colour = NA, linewidth = 0.05, alpha = 0.4) +
  scale_fill_brewer(palette = "Dark2") +
  plot_limits

# -----------------------------------------------------------------------------
# Proportion of MPA zone overlap with various layers

# Existing HBLL blocks
if (!file.exists(here::here("data-generated", "spatial", "hbll_mpa_overlap.rds"))) {
hbll_mpa_overlap <- st_intersection(
    hbll_blocks |> st_transform(st_crs(comm_ll_allowance)),
    comm_ll_allowance) %>%
  mutate(hbll_mpa_area = st_area(.))

sampled_mpa_overlap <- st_intersection(
    sampled_blocks |> st_transform(st_crs(comm_ll_allowance)),
    comm_ll_allowance) %>%
  mutate(sampled_mpa_area = st_area(.))

ll_mpa_overlap <- st_intersection(
    ye_ll_cpue_grid |> st_transform(st_crs(comm_ll_allowance)),
    comm_ll_allowance) %>%
  mutate(cpue_mpa_area = st_area(.))
saveRDS(ll_mpa_overlap, here::here("data-generated", "spatial", "ll-mpa-overlap-yelloweye.rds"))
} else {
  ll_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "ll-mpa-overlap-yelloweye.rds"))
}

if (!file.exists(here::here("data-generated", "spatial", "sampled-cpue-mpa-overlap-yelloweye.rds"))) {
  sampled_cpue_mpa_overlap <- st_intersection(
      ye_ll_cpue_grid |> st_transform(st_crs(comm_ll_allowance)),
      comm_ll_allowance) |>
      st_intersection(sampled_blocks |> st_transform(st_crs(comm_ll_allowance))) %>%
    mutate(sampled_cpue_area = st_area(.))
  saveRDS(sampled_cpue_mpa_overlap, here::here("data-generated", "spatial", "sampled-cpue-mpa-overlap-yelloweye.rds"))
} else {
  sampled_cpue_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "sampled-cpue-mpa-overlap-yelloweye.rds"))
}

# Question - I think it isn't too late to add more baseline sampling so maybe
# this is the best combination to look at for a start?
if (!file.exists(here::here("data-generated", "spatial", "hbll-cpue-mpa-overlap-yelloweye.rds"))) {
  hbll_cpue_mpa_overlap <- st_intersection(
      ye_ll_cpue_grid |> st_transform(st_crs(comm_ll_allowance)),
      comm_ll_allowance) |>
      st_intersection(hbll_blocks |> st_transform(st_crs(comm_ll_allowance))) %>%
    mutate(hbll_cpue_area = st_area(.))

  saveRDS(hbll_cpue_mpa_overlap, here::here("data-generated", "spatial", "hbll-cpue-mpa-overlap-yelloweye.rds"))
} else {
  hbll_cpue_mpa_overlap <- readRDS(here::here("data-generated", "spatial", "hbll-cpue-mpa-overlap-yelloweye.rds"))
}

# Plotting overlaps
p1 <- ggplot(data = hbll_mpa_overlap) +
  geom_sf(data = comm_ll_allowance, aes(fill = activity_allowance_label), alpha = 0.4, colour = NA) +
  geom_sf(aes(fill = activity_allowance_label)) +
  geom_sf(data = pacea::bc_coast, fill = "grey99") +
  scale_fill_manual(name = "Proposed Longline Allowance", values = c("allowed" = "#0072B2",
    "not allowed" = "#D55E00", "conditional" = "#F0E442",
    "not applicable" = "#999999"), na.value = "grey90") +
  plot_limits +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8)) +
  ggtitle("MPA zone overlap with existing HBLL blocks")
p1

p2 <- p1 %+% sampled_mpa_overlap +
  ggtitle("MPA zone overlap with sampled HBLL blocks")

p3 <- p1 %+% ll_mpa_overlap +
  ggtitle("MPA zone overlap with longline CPUE (NOT PRIVACY COMPLIANT)")

p1 / p2 / p3
ggsave(here::here("draft-figures", "mpa-overlap-hbll-cpue.pdf"), width = 9.5, height = 16.5)

p4 <- p1 %+% sampled_cpue_mpa_overlap +
  ggtitle("MPA - CPUE - sampled overlap ")
p4

p5 <- p1 %+% hbll_cpue_mpa_overlap +
  ggtitle("MPA - CPUE - HBLL overlap ")
p5

p4 / p5
ggsave(here::here("draft-figures", "combined-overlap.pdf"), width = 9.5, height = 12.5)


# p5 %+% geom_sf(data = ye_trawl_cpue_grid, fill = "grey80", colour = NA, linewidth = 0.05) +
#   plot_limits

p6 <- p5 %+% geom_sf(data = ye_ll_cpue_grid, fill = "pink", alpha = 0.4, colour = NA, linewidth = 0.05) +
  plot_limits
p6$layers <- p6$layers[c(4, 1, 2, 3)]
p6





# # Leaflet shenanigans
# # -----------------------------------------------------------------------------
# # Load libraries
# library(leaflet)
# library(leafgl) # not sure this is doign anything
# library(sf)
# library(tidyverse)
# library(pacea) # For the real bc_coast if you have it


# # Create color palettes
# pal_comm <- colorFactor(
#   palette = c("#0072B2", "#D55E00", "#F0E442", "#999999"),
#   domain = comm_ll_allowance$activity_allowance_label
# )
# pal_cpue <- colorNumeric(
#   palette = "viridis",
#   domain = log10(aggregated_grid_sf$geo_mean_cpue),
#   na.color = "grey90"
# )

# # Build the leaflet map
# leaflet_map <- leaflet() %>%
#   # Add a basemap (optional, can be toggled)
#   # addProviderTiles(providers$CartoDB.Positron, group = "Basemap") %>%

#   # Set the initial view
#   setView(lng = -129, lat = 52.5, zoom = 6) %>%

#   # Layer 0: BC Coast (using the proper sf object from pacea)
#   # This layer will be our base, so we don't put it in the toggle control
#   addPolygons(
#     data = pacea::bc_coast,
#     fillColor = "grey99",
#     color = "grey60",
#     weight = 1
#   ) %>%
#   htmlwidgets::onRender("
#     function(el, x) {
#       el.style.backgroundColor = 'white';
#     }
#   ") |>

#   # Layer 1: Commercial Allowance
#   addPolygons(
#     data = comm_ll_allowance_simplified |> st_transform(st_crs(4326)) |> st_zm(),
#     fillColor = ~pal_comm(activity_allowance_label),
#     fillOpacity = 0.7, weight = 1,
#     color = ~pal_comm(activity_allowance_label),
#     popup = ~paste("Status:", activity_allowance_label),
#     group = "Commercial Allowance"
#   ) %>%

#   # Layer 2: CPUE Grid
#   addPolygons(
#     data = ye_ll_cpue_grid |> st_transform(st_crs(4326)) |> st_zm(),
#     fillColor = ~pal_cpue(log10(geo_mean_cpue)),
#     fillOpacity = 0.8,
#     weight = 0, # No border
#     popup = ~paste("Mean CPUE:", round(geo_mean_cpue, 2)),
#     group = "Commercial CPUE (longline)"
#   ) %>%

#   # Layer 3: Human use layer
#   addPolygons(
#     data = hu_ll |> st_transform(st_crs(4326)) |> st_zm(),
#     color = "grey80",
#     fillColor = "grey80",
#     fillOpacity = 0.3,
#     group = "Human use layer"
#   ) %>%

#   # Layer 4: HBLL Survey Blocks (using addPolylines for outlines)
#   addPolylines(
#     data = hbll_blocks |> st_transform(st_crs(4326)) |> st_zm(),
#     color = "red",
#     fill = FALSE,
#     weight = 1,
#     group = "HBLL Blocks"
#   ) %>%

#   # Layer 5: Sampled Blocks
#   addPolylines(
#     data = sampled_blocks |> st_transform(st_crs(4326)) |> st_zm(),
#     # radius = 3,
#     color = "black",
#     weight = 0.5,
#     stroke = TRUE,
#     fill = FALSE,
#     group = "Sampled Blocks"
#   ) %>%

#   # Add Legends
#   addLegend(
#     pal = pal_comm,
#     values = comm_ll_allowance$activity_allowance_label,
#     title = "Commercial Allowance",
#     position = "bottomright"
#   ) %>%
#   addLegend(
#     pal = pal_cpue,
#     values = log10(aggregated_grid_sf$geo_mean_cpue),
#     title = "Geo. Mean CPUE (log10)",
#     position = "bottomright",
#     labFormat = labelFormat(transform = function(x) round(10^x, 2))
#   ) %>%

#   # Add the layer control widget
#   addLayersControl(
#     baseGroups = "Basemap",
#     overlayGroups = c("Commercial Allowance", "Commercial CPUE", "Historical Use", "HBLL Blocks", "Sampled Blocks"),
#     options = layersControlOptions(collapsed = FALSE)
#   )

# leaflet_map
