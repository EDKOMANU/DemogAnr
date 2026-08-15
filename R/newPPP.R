#' @name project_population
#' @title Probabilistic Population Projection (Balancing Equation)
#'
#' @description Projects population trajectories for several administrative
#' regions by stochastic simulation of the demographic balancing equation
#' (Preston et al. 2001, Ch. 1). The annual continuous growth rate of a
#' population is the crude birth rate minus the crude death rate plus the
#' net-migration rate, `g = b - d + m`, and the population is carried forward
#' one year at a time as `P(t+1) = P(t) * exp(g)`.
#'
#' Uncertainty is treated in the spirit of probabilistic projection (Raftery et
#' al. 2012): rather than fixing \eqn{b}, \eqn{d} and \eqn{m}, the function
#' draws each of them from a distribution. For every one of `num_samples`
#' trajectories the three rates are drawn **once** (parameter uncertainty) and
#' then held across the horizon, so that a "high-fertility" or "high-mortality"
#' world persists over time and the projection fan widens as the horizon
#' lengthens. An optional per-year innovation (`annual_noise_sd`) adds temporal
#' variability on top.
#'
#' The distribution of each component is **user-definable**. Supply a sampler
#' `function(mean, n)` returning `n` draws centred on the region's rate to
#' `birth_dist`, `death_dist`, or `migration_dist`; leave them `NULL` to use the
#' defaults (a lognormal for the non-negative birth and death rates, a normal
#' for net migration), whose spread is set by `cv`.
#'
#' @param data A data frame with one row per region.
#' @param future_year Numeric. Target (final) year of the projection.
#' @param base_year Numeric. Base (launch) year.
#' @param region_var,subregion_var Character. Column names identifying the
#'   region and subregion.
#' @param base_pop_var Character. Column name for the base-year population.
#' @param birth_rate_var Character. Column name for the crude birth rate
#'   \eqn{b} (a per-capita annual rate, e.g. 0.032 for 32 per 1000).
#' @param death_rate_var Character. Column name for the crude death rate
#'   \eqn{d} (per-capita annual rate).
#' @param net_migration_var Character. Column name for the net-migration rate
#'   \eqn{m} (per-capita annual rate; positive for net in-migration).
#' @param birth_dist,death_dist,migration_dist Optional samplers. Each is a
#'   `function(mean, n)` returning `n` random draws of the rate for a region
#'   whose central value is `mean`. When `NULL` (default), a lognormal
#'   (`birth_dist`, `death_dist`) or normal (`migration_dist`) sampler with
#'   coefficient of variation `cv` is used.
#' @param cv Numeric. Coefficient of variation of the default samplers
#'   (default 0.10). Ignored for any component given an explicit sampler.
#' @param annual_noise_sd Numeric. Standard deviation of an optional Normal
#'   innovation added to the growth rate each year (default 0, i.e. none).
#' @param num_samples Integer. Number of simulated trajectories (default 2000).
#' @param probs Numeric length-2. Lower and upper quantiles summarising the
#'   trajectories (default `c(0.1, 0.9)`, the 10th to 90th percentile band).
#' @param random_seed Integer. Seed for reproducibility (default 42).
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   projected trajectories (median and the `probs` band, faceted by subregion)
#'   is attached to the result and shown when it is printed.
#' @param verbose Logical. If `TRUE`, prints status messages during projection.
#'
#' @return A data frame with one row per region, subregion, and year, carrying
#'   the summary of the simulated trajectories: `lower` (the `probs[1]`
#'   quantile), `median`, `mean`, and `upper` (the `probs[2]` quantile). When
#'   `graph = TRUE` the result also carries an attached plot and prints through
#'   its `print` method (`print.dem_projection`).
#'
#' @examples
#' data(region2000)
#' # Default uncertainty (lognormal births/deaths, normal migration, cv = 0.10)
#' proj <- project_population(region2000,
#'   base_year = 2000, future_year = 2010,
#'   region_var = "Country", subregion_var = "Region",
#'   base_pop_var = "base_pop",
#'   birth_rate_var = "cbr", death_rate_var = "cdr",
#'   net_migration_var = "nmr",
#'   num_samples = 500, graph = FALSE
#' )
#' head(proj)
#'
#' # User-defined distributions: wider, skewed mortality; tight fertility
#' proj2 <- project_population(region2000[1:4, ],
#'   base_year = 2000, future_year = 2010,
#'   region_var = "Country", subregion_var = "Region",
#'   base_pop_var = "base_pop",
#'   birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
#'   birth_dist     = function(mean, n) rnorm(n, mean, 0.001),
#'   death_dist     = function(mean, n) rlnorm(n, log(mean), 0.25),
#'   migration_dist = function(mean, n) rnorm(n, mean, 0.004),
#'   num_samples = 500, graph = FALSE
#' )
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' ISBN 978-0631226161. (Chapter 1: the balancing equation and growth rates.)
#'
#' Raftery, A. E., Li, N., Sevcikova, H., Gerland, P., & Heilig, G. K. (2012).
#' Bayesian probabilistic population projections for all countries.
#' \emph{Proceedings of the National Academy of Sciences}, 109(35),
#' 13915-13921. \doi{10.1073/pnas.1211452109} (Probabilistic projection through
#' component uncertainty.)
#'
#' @seealso [math_project()] for deterministic subnational projection.
#' @export
#' @importFrom stats rnorm rlnorm median quantile

project_population <- function(
    data, future_year, base_year,
    region_var = "region", subregion_var = "subregion",
    base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr",
    net_migration_var = "nmr",
    birth_dist = NULL, death_dist = NULL, migration_dist = NULL,
    cv = 0.10, annual_noise_sd = 0,
    num_samples = 2000, probs = c(0.1, 0.9),
    random_seed = 42, graph = TRUE, verbose = FALSE
) {
  if (future_year <= base_year) stop("future_year must be greater than base_year.")
  if (length(probs) != 2 || any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
    stop("'probs' must be two increasing quantiles in [0, 1].")
  }
  needed <- c(region_var, subregion_var, base_pop_var,
              birth_rate_var, death_rate_var, net_migration_var)
  miss <- needed[!needed %in% names(data)]
  if (length(miss)) stop("Column(s) not found in 'data': ", paste(miss, collapse = ", "))

  set.seed(random_seed)
  num_years <- future_year - base_year

  # Build the three component samplers. A user sampler is any function(mean, n);
  # otherwise a default is used (non-negative lognormal for rates, normal for
  # net migration), with spread controlled by 'cv'.
  birth_sampler <- .rate_sampler(birth_dist, "lognormal", cv)
  death_sampler <- .rate_sampler(death_dist, "lognormal", cv)
  mig_sampler   <- .rate_sampler(migration_dist, "normal", cv)

  if (verbose) {
    message(sprintf("Projecting %d region(s), %d to %d, %d trajectories.",
                    nrow(data), base_year, future_year, num_samples))
  }

  results_list <- vector("list", nrow(data))

  for (i in seq_len(nrow(data))) {
    row_data  <- data[i, ]
    region    <- as.character(row_data[[region_var]])
    subregion <- as.character(row_data[[subregion_var]])
    base_pop  <- as.numeric(row_data[[base_pop_var]])

    b_mean <- as.numeric(row_data[[birth_rate_var]])
    d_mean <- as.numeric(row_data[[death_rate_var]])
    m_mean <- as.numeric(row_data[[net_migration_var]])

    if (verbose) {
      message(sprintf("  %s / %s: base pop %.0f, b=%.4f d=%.4f m=%.4f",
                      region, subregion, base_pop, b_mean, d_mean, m_mean))
    }

    # Parameter uncertainty: draw each component ONCE per trajectory and hold
    # it across the horizon, so structural uncertainty persists (and the fan
    # widens with time) rather than washing out into year-to-year noise.
    b <- birth_sampler(b_mean, num_samples)
    d <- death_sampler(d_mean, num_samples)
    m <- mig_sampler(m_mean, num_samples)

    sim_pop <- rep(base_pop, num_samples)

    region_results <- data.frame(
      region = region, subregion = subregion, year = base_year,
      lower = base_pop, median = base_pop, mean = base_pop, upper = base_pop,
      stringsAsFactors = FALSE
    )

    for (j in seq_len(num_years)) {
      eps <- if (annual_noise_sd > 0) rnorm(num_samples, 0, annual_noise_sd) else 0
      g <- b - d + m + eps                       # balancing equation, per year
      sim_pop <- sim_pop * exp(g)                # P_{t+1} = P_t * exp(g)

      region_results <- rbind(region_results, data.frame(
        region = region, subregion = subregion, year = base_year + j,
        lower  = as.numeric(quantile(sim_pop, probs[1])),
        median = median(sim_pop),
        mean   = mean(sim_pop),
        upper  = as.numeric(quantile(sim_pop, probs[2])),
        stringsAsFactors = FALSE
      ))
    }
    results_list[[i]] <- region_results
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  if (graph) {
    attr(results_df, "plot") <- .plot_projection(results_df)
    class(results_df) <- c("dem_projection", "data.frame")
  }
  results_df
}

# Build a component sampler. If 'dist' is a function it is returned unchanged
# (user-defined distribution). Otherwise a default centred on 'mean' with
# coefficient of variation 'cv' is returned: a non-negative lognormal for rates
# ("lognormal"), or a normal for a quantity that may be negative ("normal").
.rate_sampler <- function(dist, family = c("lognormal", "normal"), cv) {
  family <- match.arg(family)
  if (is.function(dist)) return(dist)
  if (!is.numeric(cv) || length(cv) != 1 || cv < 0) stop("'cv' must be a single non-negative number.")
  function(mean, n) {
    mean <- mean[1]
    if (is.na(mean)) return(rep(NA_real_, n))
    if (cv == 0) return(rep(mean, n))
    if (family == "lognormal") {
      if (mean <= 0) return(rep(mean, n))
      sdlog <- sqrt(log(1 + cv^2))               # preserves E[X] = mean
      stats::rlnorm(n, meanlog = log(mean) - 0.5 * sdlog^2, sdlog = sdlog)
    } else {
      stats::rnorm(n, mean = mean, sd = abs(mean) * cv)
    }
  }
}

#' @export
print.dem_projection <- function(x, ...) {
  p <- attr(x, "plot")
  y <- x
  attr(y, "plot") <- NULL
  class(y) <- "data.frame"
  print(y, ...)
  if (!is.null(p)) print(p)
  invisible(x)
}

#' @export
plot.dem_projection <- function(x, ...) {
  p <- attr(x, "plot")
  if (is.null(p)) stop("No plot available; call project_population(..., graph = TRUE).")
  p
}
