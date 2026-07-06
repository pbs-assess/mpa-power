# =============================================================================
# Link a new run_tag to an existing one
# =============================================================================
# Run this before sourcing 00-setup.R with a new run_tag, to avoid regenerating
# expensive geography/mesh prep that doesn't depend on whatever you're testing
# (e.g. phi, mpa_trend, conditioning formula).
#
# Cleaned species data, conditioning fits, AR1 parameters, and fit
# characteristics are NOT linked here -- they're regenerated fresh for every
# run_tag (cheap enough to just re-run). Only truly universal, non-formula-
# dependent prep (mesh, geography/allocations) is shared. The shared item list
# is an explicit allowlist rather than "everything except a few exclusions" --
# anything new or unexpected in base_run_tag is left out by default instead of
# silently being linked.

#' Symlink shared, non-run-specific pipeline outputs into a new run_tag
#'
#' @param base_run_tag run_tag to link from (holds the shared mesh/geography prep)
#' @param new_run_tag run_tag to link into (created if it doesn't exist)
link_run_tag <- function(base_run_tag, new_run_tag) {
  base_dir <- here::here("data-generated", base_run_tag)
  new_dir  <- here::here("data-generated", new_run_tag)
  stopifnot(dir.exists(base_dir))
  dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)

  shared_items <- c(
    "00-mesh-cache",
    "historical-locations.rds",
    "hbll-last-sampled-year.rds",
    "hbll-restricted-sf.rds",
    "hbll-allocations.rds"
  )
  shared_items <- shared_items[file.exists(file.path(base_dir, shared_items))]

  for (item in shared_items) {
    target <- file.path(new_dir, item)
    if (!file.exists(target)) {
      file.symlink(file.path("..", base_run_tag, item), target)
    }
  }

  message("Linked ", length(shared_items), " item(s) from '", base_run_tag,
          "' into '", new_run_tag, "': ", paste(shared_items, collapse = ", "))
}

# Example usage:
# link_run_tag("ms", "0-phi-crank")
# Then set run_tag <- "0-phi-crank" in sample-fit-config.R and proceed as usual.
# 01-cleaned-species-data, 01-fits, ar1-parameters.rds, and fit-characteristics.rds
# will regenerate fresh for the new run_tag (01-fit-conditioning-models.R must be
# run locally, per its server guard), as will 02-sim-data, 03-sampled-data, and
# 04-power-results.
