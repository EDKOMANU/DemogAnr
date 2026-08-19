# Controlled stochastic disaggregation of an existing national projection to
# regions ("calibration"), with raking to the national control total.
# Source: Siegel & Swanson (2004), Ch. 21; Smith, Tayman & Swanson (2013),
# Ch. 6-7 (share methods and controlling); Deming & Stephan (1940) (raking).

#' @name calibrate_regions
#' @title Calibrate Regions to a National Projection (Controlled Disaggregation)
#'
#' @description Distributes an **already existing** national projection across
#' regions, letting each region deviate from the national profile and then
#' raking the regions back onto the national control total. This is the working
#' pattern used in official projection offices: the national trajectory is
#' settled first (by a cohort-component or probabilistic model), and the
#' subnational figures are then produced as a controlled disaggregation of it,
#' so that the published regional series adds up to the published national
#' series in every year.
#'
#' It is therefore **not** a projection. [project_population()] projects each
#' region independently from the balancing equation and the regional totals are
#' whatever they turn out to be; `calibrate_regions()` takes the national path
#' as given and asks only how it is split. Two consequences follow, and both are
#' the reason for the function:
#' \describe{
#'   \item{Raking}{The regions are forced to sum exactly to the national total
#'     in every year (Deming and Stephan 1940). The multiplier that does this is
#'     reported in the output, so the adjustment is visible rather than hidden.}
#'   \item{Correlated deviation}{Regional trajectories deviate *around* a shared
#'     national path instead of wandering independently. Independent regional
#'     errors cancel on aggregation and so understate national uncertainty;
#'     here the common component is carried by the national control total and
#'     only the regional shares are stochastic.}
#' }
#'
#' @details
#' Write \eqn{P_r} for the base-year population of region \eqn{r} (within an age
#' group, when `age_col` is given), \eqn{N_y} for the national control total in
#' year \eqn{y}, and \eqn{N_0 = \sum_r P_r} for the base-year total. The
#' calculation, carried out for each of `num_samples` trajectories, is
#' \deqn{U_{r,y} = P_r \cdot (N_y / N_0) \cdot \delta_{r,y}, \qquad
#'       f_y = N_y / \sum_r U_{r,y}, \qquad
#'       \hat P_{r,y} = f_y \, U_{r,y},}
#' where \eqn{U} is the unraked allocation, \eqn{f_y} the raking factor and
#' \eqn{\delta_{r,y}} the regional deviation. The deviation is drawn **once per
#' region and trajectory** and then phased in linearly over the horizon,
#' \eqn{\delta_{r,y} = 1 + (z_r - 1) \, w_y} with \eqn{w_y} running from 0 in
#' the base year to 1 at the end of the horizon, where \eqn{z_r} is the ratio of
#' a draw to the region's own base share. Holding the draw across the horizon
#' (rather than re-drawing each year) makes a region that is set to gain share
#' keep gaining it, which is how regional divergence actually behaves; the fan
#' widens smoothly with time and the base year is reproduced exactly.
#'
#' When `deviation` is `NULL` or 0 there is no deviation, every raking factor is
#' 1, and the result is the ordinary constant-share allocation
#' \eqn{\hat P_{r,y} = N_y \, P_r / N_0}. This is the deterministic special
#' case, and is worth running first as a check.
#'
#' Only the mean of the trajectories is additive. Because a sum of quantiles is
#' not the quantile of a sum, the `lower`, `median` and `upper` columns do not
#' add to the national total; the `calibrated` column (the trajectory mean) does,
#' exactly, and every individual trajectory does too.
#'
#' @param data A data frame of base-year regional population, one row per region
#'   (or per region and age group when `age_col` is used).
#' @param region_col Character. Column name identifying the region.
#' @param pop_col Character. Column name for the base-year regional population.
#' @param national A data frame holding the national projection that is being
#'   disaggregated: one row per year (or per year and age group), with the
#'   projected total already computed.
#' @param year_col Character. Column name for the year in `national`.
#' @param total_col Character. Column name for the projected national total in
#'   `national`.
#' @param age_col Character or `NULL`. Optional age-group column, present in
#'   both `data` and `national`. When given, shares, deviation and raking are
#'   all computed within age group, so the age profile of each region is
#'   calibrated separately and the regions add to the national figure in every
#'   age group as well as in every year.
#' @param deviation The size and shape of the regional deviation. One of:
#'   `NULL` (the default; no deviation, giving the deterministic constant-share
#'   allocation), a single number read as the coefficient of variation of a
#'   mean-preserving lognormal, a numeric vector of per-region coefficients of
#'   variation (named by region, or in the row order of the regions), or a
#'   sampler `function(mean, n)` returning `n` draws of a region's share centred
#'   on `mean`, in the same form as the samplers of [project_population()].
#' @param base_year Numeric or `NULL`. The year the `data` population refers to.
#'   Defaults to the earliest year in `national`.
#' @param num_samples Integer. Number of simulated trajectories (default 1000).
#'   Ignored when `deviation` is `NULL` or 0, where one deterministic trajectory
#'   is enough.
#' @param rake Logical. If `TRUE` (default) the regions are scaled onto the
#'   national control total. If `FALSE` the unraked allocation is returned
#'   instead, with the raking factors still reported so that the size of the
#'   discrepancy can be seen.
#' @param probs Numeric length-2. Lower and upper quantiles summarising the
#'   trajectories (default `c(0.1, 0.9)`).
#' @param random_seed Integer. Seed for reproducibility (default 42).
#' @param graph Logical. If `TRUE` (default) a \pkg{ggplot2} plot of the
#'   calibrated regional series is attached and shown on printing.
#'
#' @return An object of class `dem_calibration`: a list with
#'   \describe{
#'     \item{`table`}{The full working table, one row per region and year (and
#'       age group): the base population `base_pop`, its `share` of the national
#'       base, the national control total `national`, the national growth factor
#'       `growth`, the mean deviation `deviation`, the unraked allocation
#'       `unraked`, the raking factor `rake_factor`, the calibrated population
#'       `calibrated`, and the `lower`, `median` and `upper` summary of the
#'       trajectories.}
#'     \item{`rake`}{The raking factors by year (and age group), with the unraked
#'       and control totals they reconcile.}
#'     \item{`national`}{The control totals used.}
#'     \item{`check`}{Per year (and age group), the calibrated regional sum, the
#'       control total, and their difference: zero when `rake = TRUE`.}
#'     \item{`settings`}{The base year, horizon, number of trajectories and
#'       quantiles used.}
#'   }
#'
#' @examples
#' data(region2000)
#'
#' # A national projection that has already been made elsewhere: Ghana's total
#' # population growing at 2.5 percent a year from the 2000 regional base.
#' base_total <- sum(region2000$base_pop)
#' national <- data.frame(
#'   year  = seq(2000, 2020, 5),
#'   total = round(base_total * exp(0.025 * seq(0, 20, 5)))
#' )
#'
#' # Deterministic constant-share allocation (no deviation): the regions keep
#' # their 2000 shares and every raking factor is 1.
#' calibrate_regions(region2000,
#'   region_col = "Region", pop_col = "base_pop",
#'   national = national, year_col = "year", total_col = "total",
#'   graph = FALSE
#' )
#'
#' # With deviation: each region is allowed to drift from the national profile
#' # by a lognormal with a coefficient of variation of 0.08, phased in over the
#' # horizon, and the regions are raked back onto the national total.
#' cal <- calibrate_regions(region2000,
#'   region_col = "Region", pop_col = "base_pop",
#'   national = national, year_col = "year", total_col = "total",
#'   deviation = 0.08, num_samples = 500, graph = FALSE
#' )
#' cal$check
#'
#' # Per-region deviation: a fast-urbanising region is allowed to move further
#' # from the national profile than the rest.
#' sds <- stats::setNames(rep(0.05, nrow(region2000)), region2000$Region)
#' sds["Greater Accra"] <- 0.15
#' calibrate_regions(region2000,
#'   region_col = "Region", pop_col = "base_pop",
#'   national = national, year_col = "year", total_col = "total",
#'   deviation = sds, num_samples = 500, graph = FALSE
#' )
#'
#' # By age group. The national projection is by age as well as year, and each
#' # region's age profile is calibrated to it. Here the 2010 census age
#' # distribution is split between two illustrative regions and carried forward.
#' data(gphc2010)
#' nat_age <- data.frame(
#'   Age   = rep(gphc2010$Age, 3),
#'   year  = rep(c(2010, 2015, 2020), each = nrow(gphc2010)),
#'   total = round(rep(gphc2010$Pop, 3) * rep(c(1, 1.13, 1.28), each = nrow(gphc2010)))
#' )
#' reg_age <- data.frame(
#'   Region = rep(c("North", "South"), each = nrow(gphc2010)),
#'   Age    = rep(gphc2010$Age, 2),
#'   pop    = c(round(gphc2010$Pop * 0.4), gphc2010$Pop - round(gphc2010$Pop * 0.4))
#' )
#' cal_age <- calibrate_regions(reg_age,
#'   region_col = "Region", pop_col = "pop", age_col = "Age",
#'   national = nat_age, year_col = "year", total_col = "total",
#'   deviation = 0.05, num_samples = 200, graph = FALSE
#' )
#' head(cal_age$table, 10)
#'
#' @references
#' Deming, W. E., & Stephan, F. F. (1940). On a least squares adjustment of a
#' sampled frequency table when the expected marginal totals are known.
#' \emph{The Annals of Mathematical Statistics}, 11(4), 427-444.
#' \doi{10.1214/aoms/1177731829} (Iterative proportional fitting, the origin of
#' raking to known margins.)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and
#' Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press.
#' ISBN 978-0126419559. (Subnational projection and controlling to independent
#' totals.)
#'
#' Smith, S. K., Tayman, J., & Swanson, D. A. (2013). \emph{A Practitioner's
#' Guide to State and Local Population Projections}. Dordrecht: Springer.
#' \doi{10.1007/978-94-007-7551-0} (Share-of-total methods and the practice of
#' controlling subnational forecasts to a higher-level total.)
#'
#' Raftery, A. E., Li, N., Sevcikova, H., Gerland, P., & Heilig, G. K. (2012).
#' Bayesian probabilistic population projections for all countries.
#' \emph{Proceedings of the National Academy of Sciences}, 109(35),
#' 13915-13921. \doi{10.1073/pnas.1211452109} (Why correlated rather than
#' independent subnational error matters for the aggregate.)
#'
#' @seealso [project_population()] for independent probabilistic projection of
#'   each region, and [math_project()] for deterministic share-based
#'   subnational projection.
#' @export
#' @importFrom stats median quantile setNames

calibrate_regions <- function(
    data, region_col, pop_col,
    national, year_col, total_col,
    age_col = NULL,
    deviation = NULL,
    base_year = NULL,
    num_samples = 1000,
    rake = TRUE,
    probs = c(0.1, 0.9),
    random_seed = 42, graph = TRUE
) {
  # ---- validation -----------------------------------------------------------
  need_d <- c(region_col, pop_col, age_col)
  miss <- need_d[!need_d %in% names(data)]
  if (length(miss)) stop("Column(s) not found in 'data': ", paste(miss, collapse = ", "))
  need_n <- c(year_col, total_col, age_col)
  miss <- need_n[!need_n %in% names(national)]
  if (length(miss)) stop("Column(s) not found in 'national': ", paste(miss, collapse = ", "))
  if (length(probs) != 2 || any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
    stop("'probs' must be two increasing quantiles in [0, 1].")
  }
  if (!is.numeric(num_samples) || num_samples < 1) stop("'num_samples' must be a positive integer.")

  regions <- as.character(data[[region_col]])
  base_pop <- as.numeric(data[[pop_col]])
  if (any(is.na(base_pop))) stop("Missing values in '", pop_col, "'.")
  if (any(base_pop < 0)) stop("Negative population in '", pop_col, "'.")

  years <- sort(unique(as.numeric(national[[year_col]])))
  if (is.null(base_year)) base_year <- min(years)
  horizon <- max(years) - base_year
  if (horizon < 0) stop("'base_year' is later than every year in 'national'.")

  # Age is handled as a stratum: everything below is done within stratum, so the
  # age-group and no-age cases run through exactly the same code.
  if (is.null(age_col)) {
    strata_d <- rep("_all_", nrow(data))
    strata_n <- rep("_all_", nrow(national))
  } else {
    strata_d <- as.character(data[[age_col]])
    strata_n <- as.character(national[[age_col]])
    if (!all(strata_d %in% strata_n)) {
      stop("Age group(s) in 'data' with no control total in 'national': ",
           paste(unique(strata_d[!strata_d %in% strata_n]), collapse = ", "))
    }
  }

  # ---- the deviation sampler(s) --------------------------------------------
  # One sampler per region, so that a per-region spread can be given. A single
  # number, or NULL, is recycled across regions.
  n_reg <- length(unique(regions))
  samplers <- .deviation_samplers(deviation, unique(regions))
  deterministic <- is.null(deviation) ||
    (is.numeric(deviation) && all(deviation == 0))
  n_sim <- if (deterministic) 1L else as.integer(num_samples)

  set.seed(random_seed)

  # ---- calibration, stratum by stratum -------------------------------------
  out_rows <- list()
  rake_rows <- list()

  for (s in unique(strata_d)) {
    idx_d <- which(strata_d == s)
    reg_s <- regions[idx_d]
    pop_s <- base_pop[idx_d]
    n0 <- sum(pop_s)
    if (n0 <= 0) stop("Base population sums to zero in stratum '", s, "'.")
    share_s <- pop_s / n0

    # z is the ratio of a drawn share to the region's own share: 1 means no
    # deviation. Drawn ONCE per region and trajectory and held across the
    # horizon, so a region that is set to gain share keeps gaining it.
    z <- matrix(1, nrow = length(idx_d), ncol = n_sim)
    if (!deterministic) {
      for (i in seq_along(idx_d)) {
        draw <- samplers[[reg_s[i]]](share_s[i], n_sim)
        z[i, ] <- if (share_s[i] > 0) draw / share_s[i] else 1
      }
      z[!is.finite(z)] <- 1
      z[z < 0] <- 0                       # a share cannot be negative
    }

    nat_s <- national[strata_n == s, , drop = FALSE]
    nat_y <- as.numeric(nat_s[[year_col]])
    nat_t <- as.numeric(nat_s[[total_col]])
    ord <- order(nat_y); nat_y <- nat_y[ord]; nat_t <- nat_t[ord]

    for (k in seq_along(nat_y)) {
      y <- nat_y[k]; N <- nat_t[k]
      w <- if (horizon > 0) (y - base_year) / horizon else 0
      w <- max(0, min(1, w))              # base year exactly reproduced; no
                                          # extrapolation of the deviation
      delta <- 1 + (z - 1) * w            # phased-in deviation, region x sim
      growth <- N / n0                    # national growth factor for the year
      unraked <- pop_s * growth * delta   # region x sim
      colsum <- colSums(unraked)
      f <- if (rake) N / colsum else rep(1, n_sim)
      f[!is.finite(f)] <- 1
      calibrated <- sweep(unraked, 2, f, "*")

      out_rows[[length(out_rows) + 1L]] <- data.frame(
        region      = reg_s,
        age         = if (is.null(age_col)) NA_character_ else s,
        year        = y,
        base_pop    = pop_s,
        share       = share_s,
        national    = N,
        growth      = growth,
        deviation   = rowMeans(delta),
        unraked     = rowMeans(unraked),
        rake_factor = mean(f),
        calibrated  = rowMeans(calibrated),
        lower       = apply(calibrated, 1, quantile, probs = probs[1]),
        median      = apply(calibrated, 1, median),
        upper       = apply(calibrated, 1, quantile, probs = probs[2]),
        stringsAsFactors = FALSE
      )

      rake_rows[[length(rake_rows) + 1L]] <- data.frame(
        age         = if (is.null(age_col)) NA_character_ else s,
        year        = y,
        unraked     = mean(colsum),
        control     = N,
        rake_factor = mean(f),
        stringsAsFactors = FALSE
      )
    }
  }

  tab <- do.call(rbind, out_rows)
  rk  <- do.call(rbind, rake_rows)
  if (is.null(age_col)) { tab$age <- NULL; rk$age <- NULL }
  ord <- if (is.null(age_col)) order(tab$year, tab$region) else
    order(tab$age, tab$year, tab$region)
  tab <- tab[ord, , drop = FALSE]
  rownames(tab) <- NULL
  rk <- rk[if (is.null(age_col)) order(rk$year) else order(rk$age, rk$year), ,
           drop = FALSE]
  rownames(rk) <- NULL

  # The additivity check, reported rather than asserted: the calibrated column
  # is the trajectory mean, which is exactly additive, so 'difference' is zero
  # (to floating point) whenever rake = TRUE.
  grp <- if (is.null(age_col)) list(year = tab$year) else
    list(age = tab$age, year = tab$year)
  agg <- stats::aggregate(list(regional_sum = tab$calibrated), grp, sum)
  ctl <- stats::aggregate(list(control = tab$national), grp, function(v) v[1])
  check <- merge(agg, ctl, by = names(grp))
  check$difference <- check$regional_sum - check$control
  check <- check[order(check[[length(grp)]]), , drop = FALSE]
  rownames(check) <- NULL

  out <- list(
    table = tab, rake = rk, national = national, check = check,
    settings = list(base_year = base_year, horizon = horizon,
                    num_samples = n_sim, rake = rake, probs = probs,
                    deterministic = deterministic, n_regions = n_reg,
                    by_age = !is.null(age_col))
  )
  class(out) <- "dem_calibration"
  if (graph) out$plot <- .plot_calibration(tab, deterministic)
  out
}

# Build one deviation sampler per region. 'deviation' may be NULL (no
# deviation), a single coefficient of variation, a per-region vector of
# coefficients of variation (named, or in region order), or a user sampler
# function(mean, n) shared by every region. The default family is the
# mean-preserving lognormal of .rate_sampler(), so a deviated share has the
# region's own share as its expectation and cannot go negative.
.deviation_samplers <- function(deviation, region_names) {
  n <- length(region_names)
  if (is.null(deviation)) {
    return(stats::setNames(rep(list(.rate_sampler(NULL, "lognormal", 0)), n),
                           region_names))
  }
  if (is.function(deviation)) {
    return(stats::setNames(rep(list(deviation), n), region_names))
  }
  if (!is.numeric(deviation)) {
    stop("'deviation' must be NULL, a number, a numeric vector, or a function(mean, n).")
  }
  if (any(deviation < 0, na.rm = TRUE)) stop("'deviation' must be non-negative.")
  if (length(deviation) == 1) deviation <- rep(deviation, n)
  if (length(deviation) != n) {
    stop("'deviation' must be of length 1 or one value per region (", n, ").")
  }
  if (!is.null(names(deviation))) {
    miss <- region_names[!region_names %in% names(deviation)]
    if (length(miss)) stop("No 'deviation' given for region(s): ",
                           paste(miss, collapse = ", "))
    deviation <- deviation[region_names]
  }
  stats::setNames(
    lapply(seq_len(n), function(i) .rate_sampler(NULL, "lognormal", deviation[i])),
    region_names)
}

#' @export
print.dem_calibration <- function(x, ...) {
  s <- x$settings
  cat("Regional calibration to a national control total\n")
  cat(strrep("-", 78), "\n", sep = "")
  cat(sprintf("Base year %g, horizon %g year(s), %d region(s)%s\n",
              s$base_year, s$horizon, s$n_regions,
              if (s$by_age) ", by age group" else ""))
  cat(sprintf("Deviation: %s   Raking: %s   Trajectories: %d\n",
              if (s$deterministic) "none (constant share)" else "stochastic",
              if (s$rake) "on" else "off", s$num_samples))
  cat(strrep("-", 78), "\n", sep = "")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], function(v) if (max(abs(v), na.rm = TRUE) > 100)
    round(v, 0) else round(v, 4))
  print(tab, row.names = FALSE)
  cat(strrep("-", 78), "\n", sep = "")
  cat("Raking factors\n")
  rk <- x$rake
  rk$unraked <- round(rk$unraked, 0)
  rk$rake_factor <- round(rk$rake_factor, 6)
  print(rk, row.names = FALSE)
  cat(strrep("-", 78), "\n", sep = "")
  cat(sprintf("Largest regional-sum minus control-total difference: %.6g\n",
              max(abs(x$check$difference))))
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}

#' @export
plot.dem_calibration <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call calibrate_regions(..., graph = TRUE).")
  x$plot
}
