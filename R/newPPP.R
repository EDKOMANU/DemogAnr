#' @name project_population
#' @title Basic Probabilistic Population Projection
#'
#' @description This function projects population trajectories for multiple
#' administrative regions using a simple stochastic (simulation-based)
#' exponential growth model, inspired by the probabilistic projection
#' literature (e.g. Raftery et al. 2012), but without a full Bayesian model.
#'
#' For each region, `num_samples` growth-rate trajectories are simulated. The
#' annual growth rate of each trajectory is the sum of a random baseline, a
#' fertility effect (increasing in TFR), a mortality effect (decreasing in the
#' death rate), a migration effect (with the sign of net migration), and
#' random noise. Each trajectory is propagated independently through all
#' years, so that uncertainty accumulates over the projection horizon, and the
#' results are summarised per region and year.
#'
#' @param data A data frame containing projection parameters.
#' @param future_year Numeric. The target year for the projection.
#' @param base_year Numeric. The base year for the projection.
#' @param region_var Character. Column name for the region.
#' @param subregion_var Character. Column name for the subregion.
#' @param base_pop_var Character. Column name for the base population.
#' @param TFR_var Character. Column name for the Total Fertility Rate (TFR).
#' @param death_rate_var Character. Column name for the death rate.
#' @param net_migration_var Character. Column name for net migration.
#' @param num_samples Integer. Number of simulation samples (default: 2000).
#' @param random_seed Integer. Random seed for reproducibility (default: 42).
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during projection.
#'
#' @return A data frame with one row per region, subregion, and year,
#' containing the summary statistics of the simulated population trajectories:
#' `lower` (25th percentile), `median`, `mean`, and `upper` (75th percentile).
#'
#' @examples
#' data(region2000)
#' proj <- project_population(region2000,
#'   base_year = 2000, future_year = 2009,
#'   TFR_var = "TFR", base_pop_var = "base_pop",
#'   region_var = "Country", death_rate_var = "death_rate",
#'   net_migration_var = "net_migration",
#'   subregion_var = "Region",
#'   num_samples = 1000
#' )
#' head(proj)
#'
#' @references
#' Raftery, A. E., Li, N., Sevcikova, H., Gerland, P., & Heilig, G. K. (2012). Bayesian probabilistic population projections for all countries. \emph{Proceedings of the National Academy of Sciences}, 109(35), 13915-13921. \doi{10.1073/pnas.1211452109} (Inspiration for the stochastic/probabilistic projection logic.)
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 6: Population Projections.)
#'
#' Rowland, D. T. (2003). \emph{Demographic Methods and Concepts}. Oxford: Oxford University Press. ISBN 978-0198752639. (Chapter on population projections.)
#'
#' @export
#' @importFrom stats rnorm median quantile

project_population <- function(
    data, future_year, base_year,
    region_var = "region", subregion_var = "subregion",
    base_pop_var = "base_pop", TFR_var = "TFR",
    death_rate_var = "death_rate", net_migration_var = "net_migration",
    num_samples = 2000, random_seed = 42, verbose = FALSE
) {
  set.seed(random_seed)

  if (future_year <= base_year) stop("future_year must be greater than base_year.")
  num_years <- future_year - base_year

  if (verbose) {
    message(sprintf("Starting population projection for %d regions from %d to %d...", nrow(data), base_year, future_year))
  }

  results_list <- list()

  for (i in 1:nrow(data)) {
    row_data <- data[i, ]
    region <- as.character(row_data[[region_var]])
    subregion <- as.character(row_data[[subregion_var]])

    base_pop <- as.numeric(row_data[[base_pop_var]])

    if (verbose) {
      message(sprintf("  Projecting region: %s (subregion: %s) with base pop: %.0f", region, subregion, base_pop))
    }

    # Parameters assumed constant across years in this model
    TFR <- as.numeric(row_data[[TFR_var]])
    death_rate <- as.numeric(row_data[[death_rate_var]])
    net_migration <- as.numeric(row_data[[net_migration_var]])

    # One population value per simulated trajectory; uncertainty accumulates
    # across years because each trajectory is propagated independently.
    sim_pop <- rep(base_pop, num_samples)

    # Base year row
    region_results <- data.frame(
      region = region,
      subregion = subregion,
      year = base_year,
      lower = base_pop,
      median = base_pop,
      mean = base_pop,
      upper = base_pop,
      stringsAsFactors = FALSE
    )

    # Simulate year by year
    for (j in 1:num_years) {
      beta0 <- rnorm(num_samples, mean = 0, sd = 0.01)
      beta_TFR <- rnorm(num_samples, mean = 0.004, sd = 0.016)
      beta_death <- rnorm(num_samples, mean = 0.002, sd = 0.2)
      beta_mig <- rnorm(num_samples, mean = 0.002, sd = 0.137)

      # Migration effect relative to current population, preserving its sign
      # (net out-migration reduces growth)
      migration_effect <- net_migration / sim_pop

      growth_rate <- beta0 +
        beta_TFR * log(1 + TFR) -
        beta_death * log(1 + death_rate) +
        beta_mig * sign(migration_effect) * log(1 + abs(migration_effect))

      # Additional noise to simulate variability
      noise <- rnorm(num_samples, mean = 0, sd = 0.02)
      annual_growth <- growth_rate + noise

      # Update every trajectory using exponential growth over one year
      sim_pop <- sim_pop * exp(annual_growth)

      region_results <- rbind(region_results, data.frame(
        region = region,
        subregion = subregion,
        year = base_year + j,
        lower = as.numeric(quantile(sim_pop, 0.25)),
        median = median(sim_pop),
        mean = mean(sim_pop),
        upper = as.numeric(quantile(sim_pop, 0.75)),
        stringsAsFactors = FALSE
      ))
    }

    results_list[[i]] <- region_results
  }

  # Combine all region results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  return(results_df)
}
