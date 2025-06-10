library(dplyr)
library(ggplot2)
library(sf)

theme_set(ggsidekick::theme_sleek())

draft_figs <- here::here("draft-figures")
dir.create(draft_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("data-generated", "spatial"), recursive = TRUE, showWarnings = FALSE)

# Look up table for activity allowance (from carrie but will double check for analysis)
activity_allowance_lu <- tibble::enframe(c(O = "allowed",
  X = "not allowed",
  AC = "conditional",
  na = "not applicable")) |>
  rename(activity_allowance = name, activity_allowance_label = value)

# MPA polygons
shape <- st_read(here::here("data-raw", "spatial", "Spatial_Q.gdb")) |>
  janitor::clean_names()

comm_ll <- shape |>
  select(hu_commercial_harvest_bottom_longline_demersal_hookand_line) |>
  left_join(activity_allowance_lu,
    by = c("hu_commercial_harvest_bottom_longline_demersal_hookand_line" = "activity_allowance"))

saveRDS(comm_ll, file.path("data-generated", "spatial", "comm_ll.rds"))

# Human use layers
human_layers <- st_layers(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"))
human_layers$name

hu_ll <- st_read(here::here("data-raw", "spatial", "mpatt_hu_10.gdb"),
  layer = "hu_co_demersalfishing_bottomlongline_d")

# Look at the context of the MPA network with commercial longline
ggplot(data = comm_ll |> rotate_a(a = 0)) +
  geom_sf(data = hu_ll |> rotate_a(a = 0), fill = "grey90") +
  geom_sf(aes(fill = activity_allowance_label), alpha = 0.8) +
  coord_sf(xlim = c(-134, -124), ylim = c(48.5, 54.5), crs = 4326) +
  labs(fill = "Commercial longline") +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("allowed" = "#0072B2", "not allowed" = "#D55E00", "conditional" = "#F0E442", "not applicable" = "#999999"))
ggsave(file.path(draft_figs, "human-ll-over-mpa.pdf"), width = 6.4, height = 7.5)
ggsave(file.path(draft_figs, "human-ll-over-mpa.png"), width = 6.4, height = 7.5)
