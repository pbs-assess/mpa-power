source(here::here("R", "01-fit-conditioning-models.R"))
source(here::here("R", "00-fit-sim-functions.R"))

# Check AR1 deviations function working as expected
# ------------------------------------------------
b <- get_model_pars(fit_IN)
year_effects <- b$estimate[grepl("fyear", b$term)]
last_year_mean <- year_effects[length(year_effects)]
last_year <- as.numeric(gsub("fyear", "", b$term[grepl("fyear", b$term)]))[length(year_effects)]
sigma_V <- sd(year_effects)
rho_V <- 0.2

# Create 50 replicates of AR1 simulation
n_replicates <- 50
sim_ar1_replicates <- purrr::map_dfr(1:n_replicates, function(rep_id) {
  tibble(
    replicate = rep_id,
    year = last_year + 1:20,
    estimate = sim_ar1_deviations(
      rho = rho_V, sigma = sigma_V, years = 1:20
    )
  )
}, .id = "replicate") |>
 mutate(estimate = estimate + last_year_mean)

# Plot historical data with all 50 replicates
historical_data <- b |>
  filter(grepl("fyear", term)) |>
  mutate(year = as.numeric(gsub("fyear", "", term)),
         data_type = "Historical") |>
  select(year, estimate, data_type)

ggplot() +
  geom_line(data = sim_ar1_replicates,
            aes(x = year - 1, y = estimate, group = replicate),
            alpha = 0.15, color = "steelblue") +
  geom_line(data = historical_data, aes(x = year, y = estimate),
            color = "black") +
  geom_vline(xintercept = last_year, linetype = "dotted", alpha = 0.7) +
  labs(title = "Historical Year Effects + 50 AR1 Realizations",
        subtitle = paste0("σ_V = ", round(sigma_V, 3), ", ρ = ", rho_V),
        x = "Year", y = "Log abundance (year effect)")

# Formal tests for AR(1) deviations function
# --------------------------------------------
# Mean should be approximately 0
test_mean_zero <- function() {
  set.seed(123)
  devs <- sim_ar1_deviations(rho = 0.3, sigma = 1, years = 1:100000)
  mean_dev <- mean(devs)
  stopifnot(abs(mean_dev) < 0.001)  # Should be close to 0 for large n
  cat("✓ Mean ≈ 0 test passed (mean =", round(mean_dev, 3), ")\n")
}

# Marginal variance should be approximately sigma^2
test_marginal_variance <- function() {
  set.seed(123)
  sigma <- 1.5
  devs <- sim_ar1_deviations(rho = 0.4, sigma = sigma, years = 1:10000)
  empirical_sd <- sd(devs)
  relative_error <- abs(empirical_sd - sigma) / sigma
  stopifnot(relative_error < 0.001)
  cat("✓ Marginal variance test passed (target σ =", sigma, ", empirical σ =", round(empirical_sd, 3), ")\n")
}

# Test 4: Autocorrelation should match rho
test_autocorrelation <- function() {
  set.seed(123)
  rho <- 0.7
  devs <- sim_ar1_deviations(rho = rho, sigma = 1, years = 1:100000)
  empirical_rho <- cor(devs[-length(devs)], devs[-1])
  relative_error <- abs(empirical_rho - rho) / rho
  stopifnot(relative_error < 0.001)
  cat("✓ Autocorrelation test passed (target ρ =", rho, ", empirical ρ =", round(empirical_rho, 3), ")\n")
}

test_mean_zero()
test_marginal_variance()
test_autocorrelation()
