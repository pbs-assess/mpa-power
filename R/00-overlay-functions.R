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

plot_overlays <- function(survey_mpa_df, survey_mpa_sf, title, legend_title,
                          geom_size = waiver(),
                          return_data = FALSE, breaks = NULL, labels = NULL) {

  subregions <- unique(survey_mpa_df$subregion)
  factor_subregions <- factor(subregions, levels = c("HG", "NC", "CC", "NVI"))
  num_subregions <- length(sort(factor_subregions))
  message(paste0("Number of subregions: ", num_subregions, " (", paste(sort(factor_subregions), collapse = ", "), ")"))

  mpa_summary_sf <- survey_mpa_df |>
    group_by(uid) |>
    summarise(n = n(), .groups = "drop") |>
    left_join(x = display_mpa, y = _, by = "uid")

  if (return_data) {
    return(mpa_summary_sf)
  }

  north_arrow <- function() {
    ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(0.4, "cm"),
    width = unit(0.4, "cm"),
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.2, "cm"),
    style = ggspatial::north_arrow_orienteering(fill = c("black", "black"),
    text_size = 5)
  )}

  map_theme <- function() {
    theme(panel.background = element_rect(fill = ocean_colour),
        legend.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0.1, l = 0),
        legend.box.spacing = unit(0.1, "cm"),
        legend.position = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.text.align = 0.5,
        legend.text = element_text(size = 7))
  }

  # Bin data if breaks provided
  use_bins <- !is.null(breaks)

  if (use_bins) {
    mpa_summary_sf <- mpa_summary_sf |>
      mutate(n = ifelse(is.na(n), 0, n)) |>
      mutate(n_binned = cut(n, breaks = breaks, labels = labels,
                            include.lowest = TRUE, right = TRUE))
  }

  #zero_colour <- "white"
  zero_colour <- "#feedde"
  ocean_colour <- "#b8d4ed"
  land_colour <- "grey70"

    # Grey for "0", viridis for sample counts
    colours <- if (labels[1] == "0") {
      c(zero_colour, viridis::viridis(length(labels) - 1))
    } else {
      colours <- c(zero_colour, "#2E7D32", "#FDD835", "#FB8C00", "#D32F2F")
    }

    # colours <- c(zero_colour, RColorBrewer::brewer.pal(length(labels) - 1, "YlOrRd"))
    # colours <- c(zero_colour, "#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3")
    colours <- c(zero_colour, "#fdbe85", "#fd8d3c", "#d94701", "#a63603")

  tag_theme <- theme(
    plot.tag = element_text(size = 9),
    plot.tag.position = c(0.02, 0.95)
  )

  p1 <- ggplot(data = survey_mpa_sf |> rotate_sf()) +
      geom_sf(data = ne_coast |> rotate_sf(), fill = land_colour, linewidth = 0.08) +
      geom_sf(data = mpa_summary_sf |> rotate_sf(), fill = "white",
              colour = "grey30", linewidth = 0.1) +
      geom_sf(size = geom_size) +
      north_arrow() +
      guides(fill = "none") +
      map_theme() +
      gfplot::coord_sf_auto(display_mpa |> rotate_sf(), buffer = 0) +
      ggtitle(title) +
      theme(plot.title = element_text(vjust = -7, size = 9)) +
      tag_theme

  p2 <- ggplot() +
    # geom_sf(data = ne_coast |> rotate_sf(), fill = "grey90") +
    geom_sf(data = ne_coast |> rotate_sf(), fill = land_colour, linewidth = 0.08) +
    geom_sf(data = mpa_summary_sf |> rotate_sf(),
            aes(fill = if (use_bins) n_binned else n),
            colour = "grey30", linewidth = 0.1) +
    north_arrow() +
    map_theme() +
    theme(axis.text.y = element_blank()) +
    gfplot::coord_sf_auto(display_mpa |> rotate_sf(), buffer = 0) +
    guides(fill = guide_legend(
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE                             # Keep keys in order
    )) +
    tag_theme

  if (use_bins) {
    p2 <- p2 + scale_fill_manual(values = colours, na.value = "grey60", name = legend_title)
  } else {
    p2 <- p2 + scale_fill_viridis_c(na.value = "grey60", name = legend_title)
  }
  (p1 + p2) +
    plot_annotation(tag_levels = "a", tag_prefix = "", tag_suffix = ")")
}
