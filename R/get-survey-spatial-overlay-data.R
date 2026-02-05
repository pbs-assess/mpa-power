# Appendix figures for the survey spatial overlay analysis
# 1) Large table summarising number of transects in each MPA site (zone) by subregion
# 2) Map of survey transects by survey type and subregion/MPA site (zone) (make heat maps)

source(here::here("R", "00-setup.R"))
source(here::here("R", "00-utils.R"))

library(dplyr)
library(purrr)
library(tidyr)
library(sf)
library(officer)
library(flextable)
library(patchwork)
library(ggrepel)

create_linestring <- function(lon_start, lat_start, lon_end, lat_end) {
  st_linestring(matrix(c(lon_start, lat_start, lon_end, lat_end),
                        ncol = 2, byrow = TRUE))
}


create_site_label_ranges <- function(data, site_col = site,
                                    label_col = map_label) {

# Format a vector of numbers as ranges
format_label_ranges <- function(labels) {
  # Sort labels
  labels_sorted <- sort(labels)

  # Find breaks (where diff > 1)
  breaks <- c(0, which(diff(labels_sorted) > 1), length(labels_sorted))

  # Create range strings
  ranges <- map_chr(seq_len(length(breaks) - 1), function(i) {
    start_idx <- breaks[i] + 1
    end_idx <- breaks[i + 1]

    start_val <- labels_sorted[start_idx]
    end_val <- labels_sorted[end_idx]

    if (start_val == end_val) {
      as.character(start_val)
    } else {
      paste0(start_val, "-", end_val)
    }
  })

  # Join with comma-space
  paste(ranges, collapse = ",")
}

  data |>
    arrange(subregion, {{ site_col }}, {{ label_col }}) |>
    group_by({{ site_col }}, subregion) |>
    summarise(
      zone_label = format_label_ranges({{ label_col }}),
      .groups = "drop"
    )
}


# Prepare directories
# -------------------
table_dir <- here::here("data-generated", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

overlay_dir <- here::here("data-generated", "spatial", "overlays")
dir.create(overlay_dir, recursive = TRUE, showWarnings = FALSE)

overlay_fig_dir <- here::here("figures", "spatial-overlays")
dir.create(overlay_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Load spatial data for plotting and overlay analysis
simple_coast <- pacea::bc_coast |>
  st_transform(crs = 3005) |>
  st_simplify(dTolerance = 100)

mpa_sf <- readRDS(here::here("data-generated", "spatial", "simple-analytical-mpa.rds")) |>
  mutate(site = gsub("_", " ", common_site_name_site_profile), map_id = map_label)

zone_labels <- mpa_sf |>
  st_drop_geometry() |>
  create_site_label_ranges(
    site_col = site,
    label_col = map_label
  )

all_subregion_sites <- mpa_sf |>
  st_drop_geometry() |>
  distinct(site, subregion) |>
  left_join(zone_labels, by = c("site", "subregion"))

display_mpa <- mpa_sf |> st_simplify(dTolerance = 100) |>
  left_join(zone_labels, by = c("site", "subregion"))
mpa_count <- unique(mpa_sf$uid) |> length()

subregions <- readRDS(here::here("data-generated", "spatial", "subregion-masks.rds"))
display_subregions <- subregions |> st_simplify(dTolerance = 100)

sp <- sp_to_hyphens("yelloweye rockfish")
sp_dat0 <- readRDS(file.path(synopsis_cache, paste0(sp, ".rds")))$survey_sets

msd0 <- readRDS(here::here("data-generated", "msd-catch.rds"))
msd <- distinct(msd0, trip_id, transect_site, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# Groundfish surveys
# -----------------------------------------------------------------------------
sp_dat <- sp_dat0 |>
  filter(survey_series_id.x %in% c(hbll_ssids, syn_ssids)) |>
  mutate(fo_id = case_when(
    survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S") ~ "StARGF_01",
    survey_abbrev %in% c("HBLL INS N") ~ "StARGF_03",
    survey_abbrev %in% c("SYN HS") ~ "StARGF_02",
    survey_abbrev %in% c("SYN WCVI") ~ "StARGF_04",
    survey_abbrev %in% c("SYN QCS") ~ "StARGF_07",
    survey_abbrev %in% c("SYN WCHG") ~ "StARGF_08"
  )) |>
  select(fo_id, ssid = survey_series_id.x, fe = fishing_event_id, trip_id,
    survey_abbrev, year, species = species_common_name, doorspread_m,
    latitude_start = latitude, longitude_start = longitude,
    latitude_end, longitude_end)

gf_line_geoms <- mapply(
  create_linestring,
  sp_dat$longitude_start, sp_dat$latitude_start,
  sp_dat$longitude_end, sp_dat$latitude_end,
  SIMPLIFY = FALSE
)

gf_line_sf <- st_sfc(gf_line_geoms, crs = 4326) %>%
  st_sf(sp_dat, geometry = .)

# HBLL transects
# -----------------------------------------------------------------------------
hbll_transects <- gf_line_sf |>
  filter(ssid %in% hbll_ssids) |>
  filter(survey_abbrev != "HBLL INS S") |>
  st_transform(crs = st_crs(mpa_sf))

if (!file.exists(file.path(overlay_dir, "hbll-mpa-sf.rds"))) {
hbll_mpa_sf <- st_join(hbll_transects, mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
saveRDS(hbll_mpa_sf, file.path(overlay_dir, "hbll-mpa-sf.rds"))
} else {
  hbll_mpa_sf <- readRDS(file.path(overlay_dir, "hbll-mpa-sf.rds"))
}

hbll_mpa_df <- st_drop_geometry(hbll_mpa_sf)

# Synoptic bottom trawl polygons (swath transects)
# -----------------------------------------------------------------------------
if (!file.exists(file.path(overlay_dir, "syn-mpa-sf.rds"))) {
  syn_transects <- gf_line_sf |> filter(ssid %in% syn_ssids)

  mean_doorspread <- mean(syn_transects$doorspread_m, na.rm = TRUE)
  doorspreads <- syn_transects |>
    mutate(doorspread_m = ifelse(is.na(doorspread_m), mean_doorspread, doorspread_m)) |>
    pull(doorspread_m)

  syn_polygons <- syn_transects |>
    st_transform(crs = 32609) |>
    st_buffer(dist = doorspreads) # Buffer to each doorspread (use mean if missing)

  # Create synoptic MPA overlay
  syn_mpa_sf <- syn_polygons |>
    st_transform(crs = st_crs(mpa_sf)) |>
    st_join(mpa_sf, join = st_intersects) |>
    mutate(in_mpa = ifelse(is.na(uid), 0, 1))

  syn_mpa_df <- st_drop_geometry(syn_mpa_sf) |>
    filter(survey_abbrev != "SYN WCVI")

  saveRDS(syn_mpa_sf, file.path(overlay_dir, "syn-mpa-sf.rds"))
} else {
  syn_mpa_sf <- readRDS(file.path(overlay_dir, "syn-mpa-sf.rds"))
}

# How much fishing pressure in how many MPAs? - e.g., low medium high;
# Doesn't want to mix trawl and longline fishing pressure
# How many MPAs
# Tomorrow we will talk about the different zones situation to prep the simulation scenarios


# Multi-species dive data
# ------------------------------------------------------------------------------
msd_lines <- mapply(
  create_linestring,
  msd$lon_start, msd$lat_start,
  msd$lon_end, msd$lat_end,
  SIMPLIFY = FALSE
)

msd_sf <- st_sfc(msd_lines, crs = 4326) %>%
  st_sf(msd, geometry = .)

msd_mpa_sf <- msd_sf |>
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
msd_mpa_df <- st_drop_geometry(msd_mpa_sf)


# hbll_overlay_sf <- left_join(display_mpa, hbll_overlay_df, by = "uid") |> mutate(survey = "HBLL")
# syn_overlay_sf <- left_join(display_mpa, syn_overlay_df, by = "uid") |> mutate(survey = "SYN")
# msd_overlay_sf <- left_join(display_mpa, msd_overlay_df, by = "uid") |> mutate(survey = "MSD")

# hbll_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = hbll_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = hbll_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of survey sets") +
#   theme(legend.position = "top") +
#   guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0) +
#   ggtitle("HBLL Survey Sets")

# syn_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = syn_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = syn_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of survey sets") +
#   # theme(legend.position = "inside",
#   #       legend.position.inside = c(0.1, 0.1)) +
#   theme(legend.position = "top",
#         axis.text.y = element_blank()) +
#   guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a()) +
#   ggtitle("Synoptic Survey Sets")

# msd_overlay_plot <-
# ggplot() +
#   geom_sf(data = simple_coast |> rotate_a(), fill = "grey95") +
#   geom_sf(data = msd_overlay_sf |> rotate_a(), aes(fill = n_transects)) +
#   # geom_sf(data = msd_mpa_sf |> rotate_a() |> st_centroid(),
#   #   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   scale_fill_viridis_c(na.value = "grey60", name = "Number of transects") +
#   theme(legend.position = "top",
#         axis.text.y = element_blank()) +
#   guides(fill = guide_coloursteps(title.position = "top", title.hjust = 0.5, barwidth = 10,
#     display = "rectangles")) +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0) +
#   ggtitle("Multispecies Dive Survey Sets")

# hbll_overlay_plot + syn_overlay_plot + msd_overlay_plot

# (hbll_overlay_plot +
#   geom_sf(data = hbll_mpa_sf |> rotate_a() |> st_centroid(),
#     shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# ) +
# (syn_overlay_plot +
#   geom_sf(data = syn_mpa_sf |> rotate_a() |> st_centroid(),
#    shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# ) +
# (msd_overlay_plot +
#   geom_sf(data = msd_mpa_sf |> rotate_a() |> st_centroid(),
#   shape = 19, size = 0.4, alpha = 0.4, colour = "orange") +
#   gfplot::coord_sf_auto(display_mpa |> rotate_a(), buffer = 0)
# )

# IPHC spatial overlay
# ------------------------------------------------------------------------------
iphc <- gfdata::load_iphc_dat(species = "yelloweye rockfish")
iphc_sf <- st_as_sf(iphc, coords = c("longitude", "latitude"), crs = 4326)
iphc_mpa_sf <- iphc_sf |>
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
iphc_mpa_df <- st_drop_geometry(iphc_mpa_sf)

# eDNA data - direct from Emily
# ------------------------------------------------------------------------------
edna0 <- readxl::read_excel(here::here("data-raw", "overlay", "SEAC_master_eDNA_spatial_file_updated_Dec_2024.xlsx")) |>
  janitor::clean_names()

edna_sf <- st_as_sf(edna0, coords = c("lon_dd", "lat_dd"), crs = 4326) |>
  mutate(survey = "eDNA")
edna_mpa_sf <- edna_sf |>
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1))
# edna_mpa_df <- st_drop_geometry(edna_mpa_sf)

# ROV data - Emily requested from Rob Skelly from ROV database
# ------------------------------------------------------------------------------
rov_sf <- st_read(here::here("data-raw", "overlay", "rov_dives_0", "c4e42f2a4df449ec9b70897ee43e5415.shp")) |>
  janitor::clean_names()

rov_mpa_sf <- rov_sf |>
  st_transform(crs = st_crs(mpa_sf)) |>
  st_join(mpa_sf, join = st_intersects) |>
  mutate(in_mpa = ifelse(is.na(uid), 0, 1)) |>
  mutate(year = lubridate::year(lubridate::ymd_hms(start_time)))

rov_mpa_df <- st_drop_geometry(rov_mpa_sf)

# ------------------------------------------------------------------------------
# Overlays
# ------------------------------------------------------------------------------
plot_subregion_sites <- function(subregion) {
  subregion_sf <- display_subregions |> filter(subregion == .env$subregion)
  sites_sf <- display_mpa |> filter(subregion == .env$subregion)
  mpa_centroids <- sites_sf |> rename(geometry = "Shape") |>
    arrange(subregion, site, map_label) |>
    distinct(subregion, site, .keep_all = TRUE) |>
    mutate(site = gsub("_", " ", site))
  ggplot() +
    geom_sf(data = simple_coast |> rotate_sf(), fill = "grey98") +
    geom_sf(data = subregion_sf |> st_simplify(dTolerance = 100) |> rotate_sf(), fill = "grey90", colour = "grey80",
      linewidth = 0.8) +
    geom_sf(data = sites_sf |> rotate_sf(),
      aes(fill = site)) +
    guides(fill = "none") +
    gfplot::coord_sf_auto(subregion_sf |> rotate_sf(), buffer = 5000) +
    geom_text_repel(
      data = mpa_centroids |> rotate_sf() |> st_centroid(),
      mapping = aes(label = zone_label, geometry = geometry),
      stat = "sf_coordinates",
      size = 3,
      bg.color = "white",
      bg.r = 0.2,
      direction = "both",
      max.overlaps = Inf,
      force = 20,
      force_pull = 2,
      point.padding = 0.5,
      box.padding = 0.5,
      xlim = c(-Inf, Inf),
      ylim = c(-Inf, Inf),
      segment.size = 0.25,
      min.segment.length = 0,
      seed = 1
    ) +
    theme(axis.title = element_blank(), axis.text = element_blank())
}

# https://stackoverflow.com/questions/78529290/position-dodge-doesnt-dodge-geom-sf-label
plot_subregion_sites <- function(subregion,
                                 manual_adjustments = NULL,
                                 buffer_dist = c(8000, 8000),
                                 coord_buffer = c(-20000, 20000, -20000, 20000)) {

    # 1. Prepare spatial data ----
    # Process buffer_dist: c(x, y) or single value
    if (length(buffer_dist) == 1) {
      buffer_dist <- c(buffer_dist, buffer_dist)
    }
    buffer_dist_x <- buffer_dist[1]
    buffer_dist_y <- buffer_dist[2]

    # coord_buffer: c(xleft, xright, ybottom, ytop) or single value for all sides
    if (length(coord_buffer) == 1) {
      coord_buffer <- c(-coord_buffer, coord_buffer, -coord_buffer, coord_buffer)
    } else if (length(coord_buffer) == 2) {
      # c(x, y) -> c(-x, x, -y, y)
      coord_buffer <- c(-coord_buffer[1], coord_buffer[1], -coord_buffer[2], coord_buffer[2])
    }

    subregion_sf <- display_subregions |> filter(subregion == .env$subregion)
    sites_sf <- display_mpa |> filter(subregion == .env$subregion)

    # Extract centroids and coordinates
    mpa_centroids <- sites_sf |>
      rename(geometry = "Shape") |>
      distinct(subregion, site, .keep_all = TRUE) |>
      mutate(site = gsub("_", " ", site)) |>
      rotate_sf() |>
      st_centroid() |>
      mutate(
        lon = st_coordinates(geometry)[, "X"],
        lat = st_coordinates(geometry)[, "Y"]
      )

    # 2. Calculate edge positions ----
    bbox <- st_bbox(subregion_sf |> rotate_sf())

    # Assign each label to its nearest edge and push straight out
    mpa_labels <- mpa_centroids |>
      mutate(
        # Distance from centroid to each potential edge position
        dist_left = abs(lon - (bbox["xmin"] - buffer_dist_x)),
        dist_right = abs(lon - (bbox["xmax"] + buffer_dist_x)),
        dist_top = abs(lat - (bbox["ymax"] + buffer_dist_y)),
        dist_bottom = abs(lat - (bbox["ymin"] - buffer_dist_y)),

        # Assign to whichever edge is closest
        min_dist = pmin(dist_left, dist_right, dist_top, dist_bottom),
        edge = case_when(
          min_dist == dist_left ~ "left",
          min_dist == dist_right ~ "right",
          min_dist == dist_top ~ "top",
          min_dist == dist_bottom ~ "bottom"
        ),

        # Push straight out to nearest edge (no spacing/sorting)
        xend = case_when(
          edge == "left" ~ bbox["xmin"] - buffer_dist_x,
          edge == "right" ~ bbox["xmax"] + buffer_dist_x,
          TRUE ~ lon  # Keep same x for top/bottom
        ),
        yend = case_when(
          edge == "top" ~ bbox["ymax"] + buffer_dist_y,
          edge == "bottom" ~ bbox["ymin"] - buffer_dist_y,
          TRUE ~ lat  # Keep same y for left/right
        ),

        # Text justification
        hjust = case_when(
          edge == "left" ~ 0,    # Left-justified
          edge == "right" ~ 1,   # Right-justified
          TRUE ~ 0.5             # Centered for top/bottom
        )
      )

    # 3. Apply manual overrides ----
    if (!is.null(manual_adjustments)) {
      for (zone_label in names(manual_adjustments)) {
        adj <- manual_adjustments[[zone_label]]
        idx <- which(mpa_labels$zone_label == zone_label)

        if (length(idx) > 0) {
          if (!is.null(adj$edge)) {
            mpa_labels$edge[idx] <- adj$edge
            # Recalculate xend/yend based on new edge
            mpa_labels$xend[idx] <- case_when(
              adj$edge == "left" ~ bbox["xmin"] - buffer_dist_x,
              adj$edge == "right" ~ bbox["xmax"] + buffer_dist_x,
              TRUE ~ mpa_labels$lon[idx]
            )
            mpa_labels$yend[idx] <- case_when(
              adj$edge == "top" ~ bbox["ymax"] + buffer_dist_y,
              adj$edge == "bottom" ~ bbox["ymin"] - buffer_dist_y,
              TRUE ~ mpa_labels$lat[idx]
            )
          }
          # Apply nudges (relative shifts)
          if (!is.null(adj$nudge_x)) mpa_labels$xend[idx] <- mpa_labels$xend[idx] + adj$nudge_x
          if (!is.null(adj$nudge_y)) mpa_labels$yend[idx] <- mpa_labels$yend[idx] + adj$nudge_y
          # Apply absolute positions (override)
          if (!is.null(adj$xend)) mpa_labels$xend[idx] <- adj$xend
          if (!is.null(adj$yend)) mpa_labels$yend[idx] <- adj$yend
          if (!is.null(adj$hjust)) mpa_labels$hjust[idx] <- adj$hjust
        }
      }
    }

    # 4. Build plot ----
    ggplot() +
      # Base layers
      geom_sf(data = simple_coast |> rotate_sf(), fill = "grey97", colour = "grey92") +
      geom_sf(
        data = subregion_sf |> st_simplify(dTolerance = 100) |> rotate_sf(),
        fill = "grey90", colour = "grey80", linewidth = 0.8
      ) +
      geom_sf(data = sites_sf |> rotate_sf(), aes(fill = site)) +
      guides(fill = "none") +
      # Straight connector lines from centroids to labels
      geom_segment(
        data = mpa_labels,
        aes(x = lon, y = lat, xend = xend, yend = yend),
        colour = "grey60",
        linewidth = 0.25
      ) +
      # Labels at edges
      geom_label(
        data = mpa_labels,
        aes(x = xend, y = yend, label = zone_label, hjust = hjust),
        size = 3,
        fill = "white",
        alpha = 0.9,
        label.padding = unit(0.15, "lines"),
        label.r = unit(0.15, "lines"),
        linewidth = 0
      ) +
      coord_sf(
        xlim = st_bbox(subregion_sf |> rotate_sf())[c("xmin", "xmax")] + coord_buffer[1:2],
        ylim = st_bbox(subregion_sf |> rotate_sf())[c("ymin", "ymax")] + coord_buffer[3:4],
        expand = FALSE,
        clip = "on"
      ) +
      theme(axis.title = element_blank(), axis.text = element_blank())
  }
#' @examples
#' # Switch sites to different edges
#' adjust_label(
#'   "Scott Islands" = list(edge = "top"),
#'   "Triangle Island" = list(edge = "right")
#' )
#' # Manual position override
#' adjust_label(
#'   "Some Site" = list(xend = 1500000, yend = 800000, hjust = 0)
#' )
adjust_label <- function(...) {
  adjustments <- list(...)

  # Auto-fill hjust based on edge if not specified
  for (site_name in names(adjustments)) {
    adj <- adjustments[[site_name]]

    if (!is.null(adj$edge) && is.null(adj$hjust)) {
      adj$hjust <- case_when(
        adj$edge == "left" ~ 0,
        adj$edge == "right" ~ 1,
        TRUE ~ 0.5
      )
      adjustments[[site_name]] <- adj
    }
  }

  adjustments
}

p1 <- plot_subregion_sites("HG",
  buffer_dist = c(38000, 20000),
  coord_buffer = c(-45000, 45000, -30000, 30000),
  manual_adjustments = adjust_label(
  "501" = list(nudge_y = -20000),
  "450-456" = list(nudge_x = 40000, nudge_y = -30000),
  "430-433" = list(edge = "right"),
  "400-404,410,412" = list(nudge_y = 35000),
  "438" = list(edge = "right"),
  "447" = list(edge = "top", nudge_x = 50000, nudge_y = -50000),
  "503" = list(nudge_y = -10000),
  "411,413-416" = list(nudge_y = -20000),
  "420-421,434,489-490" = list(edge = "right", nudge_y = -20000),
  # "422-424" = list(edge = "left"),
  "422-424" = list(edge = "right"),
  "491-496" = list(nudge_y = -10000),
  # "491-496" = list(edge = "right"),
  "506" = list(edge = "left"),
  "511" = list(edge = "right")
))
p1
# ggsave(file.path(overlay_fig_dir, "hg-sites.pdf"), width = 7.2, height = 7)
p2 <- plot_subregion_sites("NC",
  buffer_dist = c(25000, 23000),
  coord_buffer = c(-30000, 30000, -50000, 30000),
  manual_adjustments = adjust_label(
  "705-709" = list(nudge_y = 8000),
  "308" = list(nudge_y = 8000),
  "350-356" = list(nudge_y = 6000),
  "316" = list(edge = "right"),
  "361" = list(nudge_y = -10000),
  "305-307" = list(nudge_y = -30000),
  "700" = list(edge = "right",nudge_y = 60000),
  "317" = list(edge = "right", nudge_y = 80000),
  "315" = list(edge = "bottom", nudge_x = 40000),
  "348" = list(nudge_x = 20000, nudge_y = 20000),
  "340-342,344-347" = list(edge = "right"),
  "99" = list(nudge_x = -40000)
))
p2
# ggsave(file.path(overlay_fig_dir, "nc-sites.pdf"), width = 7.2, height = 7)
p3 <- plot_subregion_sites("CC",
  manual_adjustments = adjust_label(
  "830-831" = list(edge = "left")
))
p3
# ggsave(file.path(overlay_fig_dir, "cc-sites.pdf"), width = 7.2, height = 7)
p4 <- plot_subregion_sites("NVI",
  buffer_dist = c(8000, 8000),
  coord_buffer = c(-20000, 20000, -20000, 20000),
  manual_adjustments = adjust_label(
  "604" = list(edge = "right"),
  "605" = list(edge = "right"),
  "666,670-672" = list(edge = "bottom"),
  "686" = list(edge = "left"),
  "687" = list(edge = "left")
))
p4
# ggsave(file.path(overlay_fig_dir, "nvi-sites.pdf"), width = 7.2, height = 7)


overlay_sf  <- bind_rows(
  hbll_mpa_sf |>
    mutate(survey = "HBLL") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry),
  syn_mpa_sf |>
    mutate(survey = "SYN") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry),
  iphc_mpa_sf |>
    mutate(survey = "IPHC") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry),
  msd_mpa_sf |>
    mutate(survey = "MSD") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry),
  edna_mpa_sf |>
    mutate(survey = "eDNA") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry),
  rov_mpa_sf |>
    mutate(survey = "ROV") |>
    select(survey, uid, subregion, site, in_mpa, year, geometry)
) |>
  mutate(site = gsub("_", " ", site))

overlay_counts <- overlay_sf |>
  st_drop_geometry() |>
  group_by(survey, subregion) |>
  summarise(n_transects = n(), .groups = "drop") |>
  mutate(subregion = ifelse(is.na(subregion), "outside", subregion))

overlay_counts

plot_hg_nc <- function(survey) {
  st1 <-
  sf_dat <- overlay_sf |> st_centroid() |> filter(survey == .env$survey)
  (p1 + geom_sf(data = sf_dat |> filter(subregion == "HG") |> rotate_sf(), shape = 21) +
    gfplot::coord_sf_auto(display_mpa |> filter(subregion == "HG") |> rotate_sf(), buffer = 5000) +
    labs(subtitle = paste0("HG: ", overlay_counts |> filter(survey == .env$survey, subregion == "HG") |> pull(n_transects)))) +
  (p2 + geom_sf(data = sf_dat |> filter(subregion == "NC") |> rotate_sf(), shape = 21) +
    gfplot::coord_sf_auto(display_mpa |> filter(subregion == "NC") |> rotate_sf(), buffer = 5000) +
    labs(subtitle = paste0("NC: ", overlay_counts |> filter(survey == .env$survey, subregion == "NC") |> pull(n_transects))) ) +
  plot_annotation(title = survey)
}

plot_hg_nc("ROV")

plot_cc_nvi <- function(survey) {
  sf_dat <- overlay_sf |> st_centroid() |> filter(survey == .env$survey)
  (p3 + geom_sf(data = sf_dat |> filter(subregion == "CC") |> rotate_sf(), shape = 21) +
    gfplot::coord_sf_auto(display_mpa |> filter(subregion == "CC") |> rotate_sf(), buffer = 5000)) +
  (p4 + geom_sf(data = sf_dat |> filter(subregion == "NVI") |> rotate_sf(), shape = 21) +
    gfplot::coord_sf_auto(display_mpa |> filter(subregion == "NVI") |> rotate_sf(), buffer = 5000)) +
  plot_annotation(title = survey)
}

surveys <- c("HBLL", "SYN", "IPHC", "MSD", "eDNA")
cairo_pdf(file.path(overlay_fig_dir, "survey-plots.pdf"), width = 12, height = 7)
  walk(surveys, \(survey) {
    print(plot_hg_nc(survey))
    print(plot_cc_nvi(survey))
  })
dev.off()

overlay_df <- overlay_sf |> st_drop_geometry() |>
  mutate(subregion = factor(subregion, levels = c("HG", "NC", "CC", "NVI"))) |>
  mutate(survey = factor(survey, levels = c("HBLL", "SYN", "IPHC", "MSD", "eDNA")))

summary_df <- overlay_df |>
  # group_by(survey, subregion, site, uid) |>
  group_by(survey, subregion, site) |>
  summarise(n_transects = n(), .groups = "drop") |>
  right_join(all_subregion_sites, by = c("site", "subregion")) |>
  pivot_wider(names_from = survey, values_from = n_transects, names_vary = "slowest") |>
  select(-`NA`) |>
  mutate(subregion = factor(subregion, levels = c("HG", "NC", "CC", "NVI"))) |>
  arrange(subregion, site)

summary_df

mk_spatial_overlay_table <- function(summary_df) {
  summary_df |>
    flextable() |>
    width(width = c(subregion = 2.3, site = 6, zone_label = 3.3, HBLL = 1.5, SYN = 1.5, IPHC = 1.5, MSD = 1.5, eDNA = 1.5), unit = "cm") |>
    set_header_labels(
      subregion = "Subregion",
      site = "Site",
      zone_label = "Zones"
    ) |>
    theme_zebra(odd_header ="#fcfcfc", odd_body = "#EFEFEF", even_header = "transparent", even_body = "transparent") |>
    bg(bg = "transparent", part = "header")
}

ft1 <- summary_df |> filter(subregion == "HG") |> mk_spatial_overlay_table() |> set_caption("Haida Gwaii")
ft2 <- summary_df |> filter(subregion == "NC") |> mk_spatial_overlay_table() |> set_caption("North Coast")
ft3 <- summary_df |> filter(subregion == "CC") |> mk_spatial_overlay_table() |> set_caption("Central Coast")
ft4 <- summary_df |> filter(subregion == "NVI") |> mk_spatial_overlay_table() |> set_caption("North Vancouver Isl")

sect_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)


doc <- read_docx()
doc <- body_add_flextable(doc, value = ft1)
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, value = ft2)
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, value = ft3)
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, value = ft4)
doc <- body_set_default_section(doc, value = sect_properties)

print(doc, target = file.path(table_dir, "survey-summary-tables.docx"))



syn_flextable <- syn_mpa_df |>
  mutate(survey_name = case_when(
    survey_abbrev == "SYN HS" ~ "Hecate Strait",
    survey_abbrev == "SYN QCS" ~ "Queen Charlotte Sound",
    survey_abbrev == "SYN WCHG" ~ "West Coast Haida Gwaii",
    TRUE ~ survey_abbrev
  )) |>
  select(-fo_id) |>
  group_by(survey_name) |>
  summarise(
    sampling_frequency = "biannual",
    year_range = paste(min(year), max(year), sep = "-"),
    n_years = n_distinct(year),
    n_years_in_mpa = n_distinct(year[in_mpa == 1]),
    n_transects = n(),
    n_transects_in_mpa = sum(in_mpa),
    prop_in_mpa = round(100 * n_transects_in_mpa / n_transects, 1),
    n_mpas = n_distinct(uid, na.rm = TRUE)
  ) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(cols = -c(survey_name), names_to = "metric", values_to = "value") |>
  pivot_wider(
    names_from = survey_name,
    values_from = value,
    names_vary = "slowest"
  ) |>
  select(metric, `West Coast Haida Gwaii`, `Hecate Strait`, `Queen Charlotte Sound`) |>
  mutate(metric = metric_labels[metric]) |>
  flextable() |>
  set_header_labels(metric = "") |>
  # add_header_row(
  #   values = "DFO Synoptic Bottom Trawl Survey Coverage",
  #   colwidths = 4
  # ) |>
  # theme_vanilla() |>
  align(align = "right", part = "all") |>
  # bold(part = "header") |>
  autofit()

# Display the table
syn_flextable

syn_table_cap <-
"Table X. Summary of DFO Synoptic Bottom Trawl Survey (West Coast Haida Gwaii, Hecate Strait, Queen Charlotte Sound) coverage and overlap within MPA boundaries.
Open Canada reference: https://open.canada.ca/data/en/dataset/945e0f13-119b-451b-9038-50c6eb641aef"

# Save both tables to single Word document
read_docx() |>
  body_add_par(hbll_table_cap) |>
  body_add_flextable(hbll_flextable) |>
  body_add_break() |>
  body_add_par(syn_table_cap) |>
  body_add_flextable(syn_flextable) |>
  body_add_par("") |>
  print(target = file.path(table_dir, "survey-summary-tables.docx"))




msd_mpa_df |>
  filter(subregion == "CC") |>
  group_by(site) |>
  summarise(n_transects = n(), .groups = "drop") |>
  right_join(all_sites, by = "site")


