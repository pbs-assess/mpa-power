# =============================================================================
# Link a new run_tag to an existing one
# =============================================================================
# Run this before sourcing 00-setup.R with a new run_tag, to avoid regenerating
# expensive geography/depth/allocation prep and Stage 1 conditioning fits that
# don't depend on whatever you're changing (e.g. phi, mpa_trend). Only the
# stages you actually want fresh (usually 02-sim-data onward) get regenerated.

#' Symlink shared pipeline outputs from an existing run_tag into a new one
#'
#' @param base_run_tag run_tag to link from (already has 00-mesh-cache, 01-fits, etc.)
#' @param new_run_tag run_tag to link into (created if it doesn't exist)
link_run_tag <- function(base_run_tag, new_run_tag) {
  base_dir <- here::here("data-generated", base_run_tag)
  new_dir  <- here::here("data-generated", new_run_tag)
  stopifnot(dir.exists(base_dir))
  dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)

  # Directories independent of phi/formula: reuse wholesale
  shared_dirs <- c("00-mesh-cache", "01-fits", "01-cleaned-species-data")

  # Core files: geography/depth/allocation prep, independent of phi/formula
  shared_files <- c(
    "historical-locations.rds",
    "hbll-last-sampled-year.rds",
    "hbll-restricted-sf.rds",
    "hbll-allocations.rds",
    "ar1-parameters.rds",
    "fit-characteristics.rds"
  )

  for (item in c(shared_dirs, shared_files)) {
    target <- file.path(new_dir, item)
    if (!file.exists(target)) {
      file.symlink(file.path("..", base_run_tag, item), target)
    }
  }

  message("Linked ", length(shared_dirs), " directories and ", length(shared_files),
          " files from '", base_run_tag, "' into '", new_run_tag, "'")
}

# Example usage:
# link_run_tag("ms", "0-phi-crank")
# Then set run_tag <- "0-phi-crank" in sample-fit-config.R and proceed as usual;
# 02-sim-data, 03-sampled-data, and 04-power-results will regenerate fresh.
