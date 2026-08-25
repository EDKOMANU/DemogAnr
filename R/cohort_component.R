# Cohort-component population projection (female-dominant Leslie matrix).
# Source: Preston, Heuveline & Guillot (2001), Chapter 6.

#' Cohort-Component Population Projection
#'
#' Projects a female population forward in \eqn{n}-year steps using the
#' cohort-component method: each age group is survived forward with life-table
#' survivorship ratios, and a new youngest group is generated from age-specific
#' fertility. The projection is represented by a Leslie matrix `L`, whose
#' sub-diagonal holds the survivorship ratios \eqn{{}_nL_{x+n}/{}_nL_x} (with the
#' last closed group and the open interval pooled) and whose first row holds the
#' fertility contributions
#' \eqn{\tfrac{{}_nL_0}{2\ell_0}\,(F^f_x + F^f_{x+n}\,{}_nL_{x+n}/{}_nL_x)}.
#'
#' Age groups must be of uniform width \eqn{n} (for example five-year groups
#' `0, 5, ..., 85+`). Fertility rates are converted to female births with
#' `fraction_female`.
#'
#' @param data A data frame with one row per (uniform-width) age group.
#' @param age Column name for the (lower bound of the) age group.
#' @param population Column name for the base-year female population.
#' @param nLx Column name for female life-table person-years `nLx`.
#' @param asfr Column name for the age-specific fertility rate (0 outside the
#'   reproductive ages).
#' @param srb,fraction_female Sex ratio at birth and the derived (or supplied)
#'   fraction of births that are female.
#' @param radix Radix \eqn{\ell_0} of `nLx` (default 100000).
#' @param n Age-interval width; defaults to the spacing of `age`.
#' @param steps Number of `n`-year steps to project (default 1).
#' @param graph Logical; if `TRUE` (default) a plot of the projected total
#'   population by year is attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_ccm`: a list with the Leslie `matrix`, a
#'   `projection` data frame (population by age group and year), the total
#'   population by `year`, the implied long-run growth ratio `lambda` (dominant
#'   eigenvalue), and the base interval `n`.
#'
#' @examples
#' # Preston et al. (2001) Box 6.1: Sweden, females, baseline 1993.
#' # Projected totals: 4,449,570 (1998) and 4,478,712 (2003).
#' sweden <- data.frame(
#'   Age  = c(seq(0, 80, 5), 85),
#'   N93  = c(293395, 248369, 240012, 261346, 285209, 314388, 281290, 286923,
#'            304108, 324946, 247613, 211351, 215140, 221764, 223506, 183654,
#'            141990, 112424),
#'   Lx   = c(497487, 497138, 496901, 496531, 495902, 495168, 494213, 492760,
#'            490447, 486613, 480665, 471786, 457852, 436153, 402775, 350358,
#'            271512, 291707),
#'   Fx   = c(0, 0, 0, 0.0120, 0.0908, 0.1499, 0.1125, 0.0441, 0.0074, 0.0003,
#'            rep(0, 8)))
#' cohort_component(sweden, age = "Age", population = "N93", nLx = "Lx",
#'                  asfr = "Fx", steps = 2, graph = FALSE)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' Chapter 6.
#'
#' @seealso [stable_population()], [math_project()]
#' @export
cohort_component <- function(data, age, population, nLx, asfr,
                             srb = 1.05, fraction_female = NULL, radix = 1e5,
                             n = NULL, steps = 1, graph = TRUE) {
  ages <- as.numeric(data[[age]])
  if (is.null(fraction_female)) fraction_female <- 1 / (1 + srb)
  if (is.null(n)) {
    w <- diff(ages)
    if (length(unique(round(w, 6))) != 1)
      stop("Age groups must have uniform width; supply 'n' or use even groups.")
    n <- w[1]
  }
  k  <- length(ages)
  L  <- as.numeric(data[[nLx]])
  P0 <- as.numeric(data[[population]])
  Ff <- as.numeric(data[[asfr]]) * fraction_female
  Ff[is.na(Ff)] <- 0

  # Leslie matrix
  M <- matrix(0, k, k)
  for (i in 2:(k - 1)) M[i, i - 1] <- L[i] / L[i - 1]       # survive into next group
  open <- L[k] / (L[k - 1] + L[k])                          # pooled open interval
  M[k, k - 1] <- open
  M[k, k]     <- open
  for (a in 1:(k - 1)) {                                    # fertility (first row)
    sa <- L[a + 1] / L[a]
    M[1, a] <- (L[1] / (2 * radix)) * (Ff[a] + Ff[a + 1] * sa)
  }
  M[1, k] <- (L[1] / (2 * radix)) * Ff[k]

  # project
  proj <- matrix(0, k, steps + 1); proj[, 1] <- P0
  for (t in seq_len(steps)) proj[, t + 1] <- M %*% proj[, t]
  years <- seq(0, by = n, length.out = steps + 1)
  colnames(proj) <- paste0("y", years)

  lambda <- max(Re(eigen(M, only.values = TRUE)$values))    # long-run growth ratio
  totals <- colSums(proj)

  pdf <- data.frame(age = ages, proj); names(pdf) <- c("age", paste0("y", years))
  out <- list(matrix = M, projection = pdf,
              totals = stats::setNames(totals, paste0("y", years)),
              lambda = lambda, r_annual = log(lambda) / n, n = n, years = years)
  class(out) <- "dem_ccm"
  if (graph) {
    out$plot <- .plot_series(years, totals, ylab = "Total population",
                             title = "Projected population", geom = "col")
  }
  out
}

#' @export
print.dem_ccm <- function(x, ...) {
  cat("Cohort-component projection\n")
  cat(strrep("-", 46), "\n", sep = "")
  print(round(x$projection), row.names = FALSE)
  cat(strrep("-", 46), "\n", sep = "")
  cat("Total population by year:\n"); print(round(x$totals))
  cat(sprintf("Long-run growth ratio (lambda) = %.5f  ->  r = %.5f / year\n",
              x$lambda, x$r_annual))
  invisible(x)
}

#' @export
plot.dem_ccm <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call cohort_component(..., graph = TRUE).")
  x$plot
}
