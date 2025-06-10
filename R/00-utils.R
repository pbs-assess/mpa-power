#' Get Model Years
#'
#' Extracts years from an sdmTMB model object.
#'
#' @param fit A sdmTMB model object containing spatiotemporal information and time data.
#' @return A sorted numeric vector of years. For spatiotemporal models, includes all years.
#'   For non-spatiotemporal models, excludes extra time years.
#'
#' @examples
#' years <- get_model_years(fit)
#'
get_model_years <- function(fit) {
  if (any(fit$spatiotemporal != "off")) {
    sort(fit$time_lu$time_from_data)
  } else {
    sort(fit$time_lu$time_from_data[!fit$time_lu$extra_time])
  }
}

#' Get model parameters
#' @param fit Fitted model object
#' @return Combined fixed and random parameters
get_model_pars <- function(fit) {
  bind_rows(tidy(fit), tidy(fit, "ran_pars"))
}

prepare_sim_data <- function(dat, pred_grid, survey = c("HBLL"), cutoff = 10) {
  survey <- match.arg(survey)

  # Filter prediction grid to match survey and years
  pred_grid <- dplyr::filter(pred_grid, survey_abbrev %in% unique(dat$survey_abbrev))

  # Create mesh
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  # Return formatted data and mesh
  list(
    data = dat,
    mesh = mesh,
    pred_grid = pred_grid
  )
}
