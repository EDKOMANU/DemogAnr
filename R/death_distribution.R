# Death distribution methods: estimating the completeness of death
# registration relative to an enumerated population.
#
# Sources: Hill (1987) for the generalized growth balance; Preston & Coale
# (1982) and Bennett & Horiuchi (1981) for synthetic extinct generations;
# Hill, You & Choi (2009) for the combined procedure. The presentation
# follows Moultrie et al. (2013), Tools for Demographic Estimation.

# Shared preparation: from two censuses and intercensal deaths, build the
# average (person-year) age distribution, the annual death rate schedule, and
# the cumulated quantities above each exact age that both methods work from.
.dd_prepare <- function(data, age, pop1, pop2, deaths, t1, t2, deaths_annual) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (nm in c(age, pop1, pop2, deaths)) {
    if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
  }
  ages <- as.numeric(data[[age]])
  if (is.unsorted(ages, strictly = TRUE)) {
    stop("Ages must be the strictly increasing lower bounds of the age groups.")
  }
  N1 <- as.numeric(data[[pop1]])
  N2 <- as.numeric(data[[pop2]])
  D  <- as.numeric(data[[deaths]])
  if (any(is.na(c(ages, N1, N2, D)))) stop("Missing values are not allowed.")
  if (any(N1 <= 0) || any(N2 <= 0)) stop("Census counts must be positive.")
  if (any(D < 0)) stop("Deaths must be non-negative.")

  t <- t2 - t1
  if (!is.finite(t) || t <= 0) stop("'t2' must be later than 't1'.")

  k <- length(ages)
  n <- c(diff(ages), Inf)

  # Person-years lived in each age group over the intercensal period are
  # approximated by the geometric mean of the two enumerations, which is the
  # exact average when the group grows exponentially.
  Nbar <- sqrt(N1 * N2)
  # Deaths per year in each age group
  Dann <- if (deaths_annual) D else D / t

  # Cumulated above each exact age x
  Nbar_open <- rev(cumsum(rev(Nbar)))
  D_open    <- rev(cumsum(rev(Dann)))
  N1_open   <- rev(cumsum(rev(N1)))
  N2_open   <- rev(cumsum(rev(N2)))

  # Growth rate of the population above each exact age
  r_open <- log(N2_open / N1_open) / t

  # People reaching exact age x each year. The density just below x is
  # Nbar[i-1]/n[i-1] and just above it Nbar[i]/n[i]; their mean estimates the
  # density at x itself. The open interval has no finite width, so there only
  # the group below is used.
  entries <- rep(NA_real_, k)
  for (i in seq_len(k)[-1]) {
    below <- Nbar[i - 1] / n[i - 1]
    above <- if (is.finite(n[i])) Nbar[i] / n[i] else NA_real_
    entries[i] <- if (is.na(above)) below else (below + above) / 2
  }

  list(ages = ages, n = n, k = k, t = t, N1 = N1, N2 = N2,
       Nbar = Nbar, Dann = Dann, entries = entries,
       Nbar_open = Nbar_open, D_open = D_open,
       N1_open = N1_open, N2_open = N2_open, r_open = r_open)
}

# Rows to fit over, given an age range
.dd_fit_rows <- function(ages, age_range, k) {
  keep <- which(ages >= age_range[1] & ages <= age_range[2])
  keep <- keep[keep >= 2 & keep <= k - 1]   # need a group below and above
  if (length(keep) < 3) {
    stop("Fewer than three age groups fall in 'age_range' (",
         age_range[1], " to ", age_range[2], "). Widen it, or supply more ",
         "age groups.")
  }
  keep
}


#' Generalized Growth Balance Estimate of Death Registration Completeness
#'
#' Estimates the completeness of death registration relative to an enumerated
#' population, and the change in census coverage between two censuses, by the
#' generalized growth balance method (Hill 1987). Where registration is
#' incomplete -- the normal case wherever civil registration is still
#' developing -- a life table built directly from registered deaths
#' understates mortality, and the completeness factor from this method is the
#' correction to apply first.
#'
#' @details
#' For a population open only at age 0, the rate at which people enter the age
#' segment \eqn{x+} equals the rate at which they leave it by death plus the
#' rate at which the segment grows:
#' \deqn{b(x+) = d(x+) + r(x+).}
#' If deaths are registered with completeness \eqn{c}, and the two censuses
#' enumerate with relative coverage \eqn{k_1/k_2}, the observed quantities
#' satisfy
#' \deqn{b(x+) - r(x+) = \frac{1}{t}\log\frac{k_1}{k_2}
#'   + \frac{1}{c}\,d^{obs}(x+),}
#' so a straight line fitted to \eqn{b(x+) - r(x+)} against \eqn{d^{obs}(x+)}
#' over a range of ages gives the completeness as the reciprocal of the slope
#' and the relative census coverage from the intercept.
#'
#' Person-years by age group are taken as the geometric mean of the two
#' enumerations, and the number of people reaching exact age \eqn{x} each year
#' as the mean of the two adjacent age groups divided by the group width.
#'
#' The fit is conventionally restricted to adult ages -- `age_range` defaults
#' to 5 to 65 -- because migration distorts the young adult ages and age
#' misreporting the oldest. Inspect `plot()` before trusting the result: the
#' points should lie close to a straight line, and a systematic curve is a
#' sign that the method's assumptions (a closed population, constant
#' completeness by age) do not hold.
#'
#' @param data A data frame with one row per age group.
#' @param age Column name for the lower bound of each age group.
#' @param pop1,pop2 Column names for the population enumerated at the first and
#'   second census.
#' @param deaths Column name for deaths by age. Taken as total deaths over the
#'   intercensal period unless `deaths_annual = TRUE`.
#' @param t1,t2 Decimal years of the two censuses (for example 2000.3 and
#'   2010.4).
#' @param deaths_annual Logical; `TRUE` if `deaths` are already an annual
#'   average rather than the intercensal total (default `FALSE`).
#' @param age_range Two ages giving the range of exact ages \eqn{x} to fit the
#'   line over (default `c(5, 65)`).
#' @param graph Logical. If `TRUE` (default), the diagnostic plot of the fitted
#'   line is attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_ggb`: a list with the estimated
#'   `completeness` of death registration, the `coverage_ratio` \eqn{k_1/k_2},
#'   the fitted `slope` and `intercept`, the `r_squared` of the fit, the
#'   `age_range` used, and a `table` of the quantities at every age.
#'
#' @examples
#' # A population whose deaths are registered at 70 per cent completeness.
#' # The censuses are ten years apart and enumerate equally well, so all three
#' # methods should return a completeness near 0.70.
#' d <- data.frame(
#'   age = c(0, seq(5, 85, 5)),
#'   n1  = c(41255, 36289, 31948, 28121, 24743, 21758, 19116, 16770, 14677,
#'           12798, 11094, 9527, 8060, 6659, 5298, 3970, 2706, 2614),
#'   n2  = c(52973, 46596, 41022, 36108, 31770, 27938, 24545, 21533, 18846,
#'           16433, 14245, 12233, 10349, 8550, 6802, 5098, 3475, 3356),
#'   dth = c(516, 134, 126, 123, 126, 136, 154, 184, 229, 294, 385, 508,
#'           668, 865, 1084, 1282, 1380, 2833)
#' )
#' ggb(d, age = "age", pop1 = "n1", pop2 = "n2", deaths = "dth",
#'        t1 = 2000, t2 = 2010, graph = FALSE)
#'
#' @references
#' Hill, K. (1987). Estimating census and death registration completeness. \emph{Asian and Pacific Population Forum}, 1(3), 8-13, 23-24.
#'
#' Brass, W. (1975). \emph{Methods for Estimating Fertility and Mortality from Limited and Defective Data}. Chapel Hill: Carolina Population Center.
#'
#' Hill, K., You, D., & Choi, Y. (2009). Death distribution methods for estimating adult mortality. \emph{Demographic Research}, 21, 235-254. \doi{10.4054/DemRes.2009.21.9}
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: IUSSP. \url{https://demographicestimation.iussp.org/}
#'
#' @seealso [seg()] for the synthetic extinct generations estimate, and
#'   [ggb_seg()] for the combined procedure.
#' @export
ggb <- function(data, age, pop1, pop2, deaths, t1, t2,
                deaths_annual = FALSE, age_range = c(5, 65), graph = TRUE) {
  p <- .dd_prepare(data, age, pop1, pop2, deaths, t1, t2, deaths_annual)
  k <- p$k; ages <- p$ages; Nbar <- p$Nbar

  # People reaching exact age x each year: the mean of the age groups either
  # side of x, divided by the width of a group.
  entries <- p$entries
  b_open <- entries / p$Nbar_open          # entry rate into x+
  d_open <- p$D_open / p$Nbar_open         # observed death rate above x
  y <- b_open - p$r_open                   # the left-hand side of the relation

  keep <- .dd_fit_rows(ages, age_range, k)
  fit <- stats::lm(y[keep] ~ d_open[keep])
  intercept <- unname(stats::coef(fit)[1])
  slope     <- unname(stats::coef(fit)[2])
  if (!is.finite(slope) || slope <= 0) {
    stop("The fitted slope is not positive, so no completeness can be read ",
         "from it. Check the inputs and the 'age_range'.")
  }

  completeness <- 1 / slope
  coverage <- exp(p$t * intercept)         # k1 / k2

  tab <- data.frame(
    age = ages, N1 = p$N1, N2 = p$N2, Nbar = Nbar,
    deaths_annual = p$Dann,
    N_open = p$Nbar_open, D_open = p$D_open,
    entry_rate = b_open, growth_rate = p$r_open, death_rate = d_open,
    y = y, fitted = ages %in% ages[keep]
  )

  out <- list(completeness = completeness, coverage_ratio = coverage,
              slope = slope, intercept = intercept,
              r_squared = summary(fit)$r.squared,
              age_range = age_range, t = p$t, model = fit, table = tab)
  class(out) <- "dem_ggb"
  if (graph) out$plot <- .plot_ggb(tab, intercept, slope, age_range)
  out
}

#' @export
print.dem_ggb <- function(x, ...) {
  cat("Generalized growth balance (Hill 1987)\n")
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("Fitted over exact ages %g to %g, %d points, R-squared %.4f\n",
              x$age_range[1], x$age_range[2], sum(x$table$fitted),
              x$r_squared))
  cat(sprintf("  slope %.5f, intercept %.6f\n", x$slope, x$intercept))
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("  Completeness of death registration: %.1f%%\n",
              100 * x$completeness))
  cat(sprintf("  Census coverage ratio k1/k2:        %.4f\n",
              x$coverage_ratio))
  cat(sprintf("  Deaths should be inflated by a factor of %.3f\n",
              1 / x$completeness))
  cat("\nInspect plot() before use: the points should lie on a straight line.\n")
  invisible(x)
}

#' @export
plot.dem_ggb <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call ggb(..., graph = TRUE).")
  x$plot
}


#' Synthetic Extinct Generations Estimate of Death Registration Completeness
#'
#' Estimates the completeness of death registration by the synthetic extinct
#' generations method of Preston and Coale (1982), in the form generalized to
#' non-stable populations by Bennett and Horiuchi (1981).
#'
#' @details
#' The population aged \eqn{x} can be reconstructed from the deaths that the
#' cohort will go on to experience, provided those deaths are inflated for the
#' growth of the population:
#' \deqn{\hat{N}(x) = \sum_{a \ge x} {}_nD(a)\,
#'   \exp\!\big(r(a)\,(a + n/2 - x)\big),}
#' where \eqn{r(a)} is the intercensal growth rate of age group \eqn{a}. This
#' is the number at *exact* age \eqn{x}, so it is compared with the census
#' estimate of the same quantity -- the mean of the densities in the age groups
#' either side of \eqn{x} -- and not with the population above \eqn{x}. If
#' registration were complete the two would agree; where deaths are registered
#' to a degree \eqn{c}, the ratio \eqn{\hat{N}(x)/N(x)} estimates \eqn{c}.
#' The reported completeness is the average of that ratio over `age_range`.
#'
#' For a stable population the relation is exact: substituting
#' \eqn{N(a) = B e^{-ra} l(a)} and \eqn{D(a) = N(a)\mu(a)} into the sum
#' returns \eqn{B e^{-rx} l(x) = N(x)}.
#'
#' Unlike [ggb()], this method assumes the two censuses cover the population
#' equally well. When they do not, run [ggb_seg()], which corrects the relative
#' coverage first.
#'
#' @inheritParams ggb
#' @param age_range Two ages giving the range of exact ages \eqn{x} over which
#'   the completeness ratios are averaged (default `c(5, 65)`).
#' @param open_mean_age Mean age at death in the open-ended interval. When
#'   `NULL` (default) it is taken as \eqn{x + 1/M}, assuming deaths decline
#'   exponentially above the open age. Using the interval's nominal midpoint
#'   instead biases the completeness downwards by several per cent, so this is
#'   worth setting where the open age is low or mortality there is known.
#'
#' @return An object of class `dem_seg`: a list with the estimated
#'   `completeness`, its standard deviation across the ages averaged
#'   (`completeness_sd`, a rough indication of how stable the estimate is),
#'   the `age_range` used, and a `table` of the age-specific ratios.
#'
#' @examples
#' # A population whose deaths are registered at 70 per cent completeness.
#' # The censuses are ten years apart and enumerate equally well, so all three
#' # methods should return a completeness near 0.70.
#' d <- data.frame(
#'   age = c(0, seq(5, 85, 5)),
#'   n1  = c(41255, 36289, 31948, 28121, 24743, 21758, 19116, 16770, 14677,
#'           12798, 11094, 9527, 8060, 6659, 5298, 3970, 2706, 2614),
#'   n2  = c(52973, 46596, 41022, 36108, 31770, 27938, 24545, 21533, 18846,
#'           16433, 14245, 12233, 10349, 8550, 6802, 5098, 3475, 3356),
#'   dth = c(516, 134, 126, 123, 126, 136, 154, 184, 229, 294, 385, 508,
#'           668, 865, 1084, 1282, 1380, 2833)
#' )
#' seg(d, age = "age", pop1 = "n1", pop2 = "n2", deaths = "dth",
#'        t1 = 2000, t2 = 2010, graph = FALSE)
#'
#' @references
#' Preston, S. H., & Coale, A. J. (1982). Age structure, growth, attrition and accession: a new synthesis. \emph{Population Index}, 48(2), 217-259. \doi{10.2307/2735961}
#'
#' Bennett, N. G., & Horiuchi, S. (1981). Estimating the completeness of death registration in a closed population. \emph{Population Index}, 47(2), 207-221. \doi{10.2307/2736447}
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: IUSSP. \url{https://demographicestimation.iussp.org/}
#'
#' @seealso [ggb()], and [ggb_seg()] for the combined procedure.
#' @export
seg <- function(data, age, pop1, pop2, deaths, t1, t2,
                deaths_annual = FALSE, age_range = c(5, 65),
                open_mean_age = NULL, graph = TRUE) {
  p <- .dd_prepare(data, age, pop1, pop2, deaths, t1, t2, deaths_annual)
  k <- p$k; ages <- p$ages

  # Growth rate within each age group, used to inflate that group's deaths
  r_group <- log(p$N2 / p$N1) / p$t

  # Mean age at death within each group. For closed groups the midpoint is
  # accurate to a tenth of a year. The open interval is different: deaths
  # there are spread over decades, and its nominal midpoint biases the
  # completeness downwards by several per cent. Assuming deaths decline
  # exponentially above the open age puts their mean at x + 1/M -- but the
  # observed rate there is itself deflated by the completeness being
  # estimated (M_obs = c M), so mean age and completeness determine each
  # other and the pair is solved by iteration.
  mid <- ages + p$n / 2
  keep <- .dd_fit_rows(ages, age_range, k)
  m_open <- if (p$Nbar[k] > 0) p$Dann[k] / p$Nbar[k] else NA_real_

  # N(x) is the number at EXACT age x that the later deaths of this cohort
  # imply, so it must be compared with the census estimate of the same thing
  # -- not with the population above x, which spans decades.
  estimate <- function(mid_open) {
    m <- mid; m[k] <- mid_open
    Nhat <- vapply(seq_len(k), function(i) {
      j <- i:k
      sum(p$Dann[j] * exp(r_group[j] * (m[j] - ages[i])))
    }, numeric(1))
    list(Nhat = Nhat, ratio = Nhat / p$entries)
  }

  fixed_open <- !is.null(open_mean_age)
  open_age <- if (fixed_open) {
    as.numeric(open_mean_age)[1]
  } else {
    ages[k] + (if (k >= 2) p$n[k - 1] else 5) / 2
  }
  res <- estimate(open_age)
  if (!fixed_open && is.finite(m_open) && m_open > 0) {
    cc <- 1
    for (it in 1:100) {
      open_age <- ages[k] + cc / m_open
      res <- estimate(open_age)
      new_c <- mean(res$ratio[keep])
      if (!is.finite(new_c) || new_c <= 0) break
      done <- abs(new_c - cc) < 1e-10
      cc <- new_c
      if (done) break
    }
  }

  Nhat <- res$Nhat
  ratio <- res$ratio
  completeness <- mean(ratio[keep])
  if (!is.finite(completeness) || completeness <= 0) {
    stop("The estimated completeness is not positive. Check the inputs and ",
         "the 'age_range'.")
  }

  tab <- data.frame(
    age = ages, N_exact = p$entries, D_open = p$D_open,
    N_hat = Nhat, ratio = ratio,
    fitted = ages %in% ages[keep]
  )

  out <- list(completeness = completeness,
              completeness_sd = stats::sd(ratio[keep]),
              open_mean_age = open_age,
              age_range = age_range, t = p$t, table = tab)
  class(out) <- "dem_seg"
  if (graph) out$plot <- .plot_seg(tab, completeness, age_range)
  out
}

#' @export
print.dem_seg <- function(x, ...) {
  cat("Synthetic extinct generations (Preston-Coale; Bennett-Horiuchi)\n")
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("Averaged over exact ages %g to %g (%d points)\n",
              x$age_range[1], x$age_range[2], sum(x$table$fitted)))
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("  Completeness of death registration: %.1f%% (sd %.3f)\n",
              100 * x$completeness, x$completeness_sd))
  cat(sprintf("  Deaths should be inflated by a factor of %.3f\n",
              1 / x$completeness))
  cat("\nA large sd, or a trend in plot(), means the ratios are not flat\n")
  cat("and the estimate should not be relied on.\n")
  invisible(x)
}

#' @export
plot.dem_seg <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call seg(..., graph = TRUE).")
  x$plot
}


#' Combined Growth Balance and Extinct Generations Estimate
#'
#' Runs the two death distribution methods in the sequence recommended by Hill,
#' You and Choi (2009): the generalized growth balance first, to estimate how
#' the coverage of the two censuses differs, and then synthetic extinct
#' generations on censuses adjusted for that difference. [seg()] on its own
#' assumes the two censuses enumerate equally well, which is what this
#' correction removes.
#'
#' @inheritParams ggb
#' @inheritParams seg
#' @param age_range Two ages for both stages (default `c(5, 65)`).
#'
#' @return An object of class `dem_ggb_seg`: a list with the final
#'   `completeness`, the `coverage_ratio` from the growth balance stage, and
#'   the two component fits, `ggb` and `seg`, the latter run on the adjusted
#'   censuses.
#'
#' @examples
#' # A population whose deaths are registered at 70 per cent completeness.
#' # The censuses are ten years apart and enumerate equally well, so all three
#' # methods should return a completeness near 0.70.
#' d <- data.frame(
#'   age = c(0, seq(5, 85, 5)),
#'   n1  = c(41255, 36289, 31948, 28121, 24743, 21758, 19116, 16770, 14677,
#'           12798, 11094, 9527, 8060, 6659, 5298, 3970, 2706, 2614),
#'   n2  = c(52973, 46596, 41022, 36108, 31770, 27938, 24545, 21533, 18846,
#'           16433, 14245, 12233, 10349, 8550, 6802, 5098, 3475, 3356),
#'   dth = c(516, 134, 126, 123, 126, 136, 154, 184, 229, 294, 385, 508,
#'           668, 865, 1084, 1282, 1380, 2833)
#' )
#' ggb_seg(d, age = "age", pop1 = "n1", pop2 = "n2", deaths = "dth",
#'            t1 = 2000, t2 = 2010, graph = FALSE)
#'
#' @references
#' Hill, K., You, D., & Choi, Y. (2009). Death distribution methods for estimating adult mortality. \emph{Demographic Research}, 21, 235-254. \doi{10.4054/DemRes.2009.21.9}
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: IUSSP. \url{https://demographicestimation.iussp.org/}
#'
#' @seealso [ggb()], [seg()]
#' @export
ggb_seg <- function(data, age, pop1, pop2, deaths, t1, t2,
                    deaths_annual = FALSE, age_range = c(5, 65),
                    open_mean_age = NULL, graph = TRUE) {
  g <- ggb(data, age, pop1, pop2, deaths, t1, t2, deaths_annual,
           age_range, graph = graph)

  # Bring the second census onto the coverage basis of the first: if the GGB
  # says k1/k2 = 1.02, the second census under-enumerates by that factor
  # relative to the first, so scaling it up makes the pair consistent.
  adj <- data
  adj[[pop2]] <- as.numeric(data[[pop2]]) * g$coverage_ratio

  s <- seg(adj, age, pop1, pop2, deaths, t1, t2, deaths_annual,
           age_range, open_mean_age, graph = graph)

  out <- list(completeness = s$completeness,
              coverage_ratio = g$coverage_ratio,
              age_range = age_range, ggb = g, seg = s)
  class(out) <- "dem_ggb_seg"
  out
}

#' @export
print.dem_ggb_seg <- function(x, ...) {
  cat("Combined growth balance and extinct generations (Hill, You & Choi 2009)\n")
  cat(strrep("-", 70), "\n", sep = "")
  cat(sprintf("  Stage 1, growth balance: census coverage ratio k1/k2 = %.4f\n",
              x$coverage_ratio))
  cat(sprintf("           (its own completeness estimate was %.1f%%)\n",
              100 * x$ggb$completeness))
  cat(sprintf("  Stage 2, extinct generations on the adjusted censuses:\n"))
  cat(strrep("-", 70), "\n", sep = "")
  cat(sprintf("  Completeness of death registration: %.1f%% (sd %.3f)\n",
              100 * x$completeness, x$seg$completeness_sd))
  cat(sprintf("  Deaths should be inflated by a factor of %.3f\n",
              1 / x$completeness))
  invisible(x)
}
