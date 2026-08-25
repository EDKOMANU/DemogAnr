# Reproduction measures and the stable population model.
# Sources: Preston, Heuveline & Guillot (2001), Chapters 5 and 7.

#' Reproduction Measures: GRR, NRR, and the Intrinsic Growth Rate
#'
#' Computes the gross and net reproduction rates from an age schedule of
#' fertility and a female life table, together with the total fertility rate,
#' the mean age of childbearing, and an approximate intrinsic growth rate.
#'
#' The gross reproduction rate is \eqn{GRR = \sum_x n \cdot F^f_x}, the average
#' number of daughters a woman would bear if she survived through her
#' reproductive span; the net reproduction rate
#' \eqn{NRR = \sum_x F^f_x \cdot {}_nL_x / \ell_0} discounts that figure for
#' mortality. \eqn{F^f_x} is the female maternity rate, obtained from the
#' age-specific fertility rate by the fraction of births that are female.
#'
#' @param data A data frame with one row per reproductive age group.
#' @param age Column name for the (lower bound of the) age group.
#' @param asfr Column name for the age-specific fertility rate. Supply either
#'   `asfr`, or both `births` and `women`.
#' @param births,women Column names for births and mid-year women, used to
#'   form the rate when `asfr` is not given.
#' @param nLx Column name for the female life-table person-years `nLx`
#'   (radix `radix`). Required for the NRR; if omitted only the GRR and TFR are
#'   returned.
#' @param radix The life-table radix \eqn{\ell_0} of `nLx` (default 100000).
#' @param srb Sex ratio at birth (default 1.05), used to derive the default
#'   `fraction_female`.
#' @param fraction_female Proportion of births that are female. Defaults to
#'   `1/(1 + srb)`. Set to 1 if the supplied rates already refer to female
#'   births only.
#' @param n Age-interval width(s); a scalar or per-row vector. Defaults to the
#'   spacing of `age`.
#' @param graph Logical; if `TRUE` (default) a plot of the net maternity
#'   function is attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_reproduction`: a list with `TFR`, `GRR`,
#'   `NRR`, `mean_age` (mean age of childbearing in the net schedule),
#'   `r_approx` (\eqn{\ln(NRR)/A}), and a `table` of the per-age contributions.
#'
#' @examples
#' # Preston et al. (2001) Box 5.5: United States, 1991 (female births)
#' us1991 <- data.frame(
#'   age   = seq(10, 45, 5),
#'   fbir  = c(5816, 253979, 532712, 596823, 431694, 162005, 25531, 829),
#'   women = c(8620, 8371, 9419, 10325, 11125, 10344, 9496, 7188) * 1000,
#'   nLx   = c(494603, 493804, 492552, 491138, 489356, 486941, 483577, 478475)
#' )
#' reproduction(us1991, age = "age", births = "fbir", women = "women",
#'              nLx = "nLx", radix = 1e5, fraction_female = 1, graph = FALSE)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' Chapter 5.
#'
#' @seealso [stable_population()]
#' @export
reproduction <- function(data, age, asfr = NULL, births = NULL, women = NULL,
                         nLx = NULL, radix = 1e5, srb = 1.05,
                         fraction_female = NULL, n = NULL, graph = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  ages <- as.numeric(data[[age]])
  if (is.null(fraction_female)) fraction_female <- 1 / (1 + srb)

  if (!is.null(asfr)) {
    fx <- as.numeric(data[[asfr]])
  } else if (!is.null(births) && !is.null(women)) {
    fx <- as.numeric(data[[births]]) / as.numeric(data[[women]])
  } else {
    stop("Supply either 'asfr', or both 'births' and 'women'.")
  }

  if (is.null(n)) n <- .interval_widths(ages)
  if (length(n) == 1) n <- rep(n, length(ages))

  matf <- fx * fraction_female            # female maternity rate
  TFR  <- sum(n * fx)
  GRR  <- sum(n * matf)

  tbl <- data.frame(age = ages, n = n, asfr = fx, maternity_f = matf)

  if (!is.null(nLx)) {
    Lx  <- as.numeric(data[[nLx]]) / radix
    net <- matf * Lx
    NRR <- sum(net)
    mid <- ages + n / 2
    A   <- sum(mid * net) / NRR
    r   <- log(NRR) / A
    tbl$nLx <- as.numeric(data[[nLx]])
    tbl$net_contribution <- net
  } else {
    NRR <- NA_real_; A <- NA_real_; r <- NA_real_
  }

  out <- list(TFR = TFR, GRR = GRR, NRR = NRR, mean_age = A,
              r_approx = r, radix = radix, fraction_female = fraction_female,
              table = tbl)
  class(out) <- "dem_reproduction"
  if (graph && !is.null(nLx)) {
    out$plot <- .plot_series(tbl$age, tbl$net_contribution,
                             ylab = "Net maternity contribution",
                             title = "Net maternity function", geom = "col")
  }
  out
}

#' Stable Population Model
#'
#' Solves for the intrinsic growth rate implied by a set of age-specific
#' fertility and mortality schedules, and derives the corresponding stable
#' population: the intrinsic birth and death rates, the stable age
#' distribution, and its mean age.
#'
#' The intrinsic growth rate \eqn{r} is the root of Lotka's equation
#' \eqn{\sum_x e^{-r(x + n/2)}\, ({}_nL_x/\ell_0)\, F^f_x = 1}, found by Coale's
#' iterative procedure. The stable birth rate is then
#' \eqn{b = 1 / \sum_x e^{-r(x + n/2)}\, {}_nL_x/\ell_0}, the death rate is
#' \eqn{d = b - r}, and the stable age distribution is
#' \eqn{{}_nc_x = b\, e^{-r(x + n/2)}\, {}_nL_x/\ell_0}.
#'
#' @param data A data frame with one row per age group, spanning all ages.
#' @param age Column name for the (lower bound of the) age group.
#' @param asfr Column name for the age-specific fertility rate (0 outside the
#'   reproductive ages).
#' @param nLx Column name for the female life-table person-years `nLx`.
#' @param radix The life-table radix of `nLx` (default 100000).
#' @param srb,fraction_female Sex ratio at birth and the derived (or supplied)
#'   fraction of births that are female; see [reproduction()].
#' @param n Age-interval width(s); defaults to the spacing of `age`.
#' @param tol,max_iter Convergence tolerance and iteration cap for the Coale
#'   procedure.
#' @param graph Logical; if `TRUE` (default) a plot of the stable age
#'   distribution is attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_stable`: a list with the intrinsic rate `r`,
#'   birth rate `b`, death rate `d`, net reproduction rate `NRR`, `mean_age` of
#'   the stable population, number of `iterations`, whether the solution
#'   `converged`, and a `table` of the stable
#'   age distribution `cx`.
#'
#' @examples
#' # Preston et al. (2001) Box 7.1: Egypt, 1997 (intrinsic r = 0.01424)
#' egypt <- data.frame(
#'   age  = seq(15, 45, 5),
#'   Lx   = c(4.66740, 4.63097, 4.58518, 4.53206, 4.46912, 4.39135, 4.28969),
#'   ma   = c(0.00567, 0.06627, 0.11204, 0.07889, 0.05075, 0.01590, 0.00610)
#' )
#' stable_population(egypt, age = "age", asfr = "ma", nLx = "Lx",
#'                   radix = 1, fraction_female = 1, graph = FALSE)$r
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' Chapter 7. Coale, A. J. (1957). A new method for calculating Lotka's r.
#' \emph{Population Studies}, 11(1), 92-94.
#'
#' @seealso [reproduction()], [lifetable()]
#' @export
stable_population <- function(data, age, asfr, nLx, radix = 1e5,
                              srb = 1.05, fraction_female = NULL, n = NULL,
                              tol = 1e-9, max_iter = 100, graph = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  ages <- as.numeric(data[[age]])
  if (is.null(fraction_female)) fraction_female <- 1 / (1 + srb)
  if (is.null(n)) n <- .interval_widths(ages)
  if (length(n) == 1) n <- rep(n, length(ages))

  Lx   <- as.numeric(data[[nLx]]) / radix
  matf <- as.numeric(data[[asfr]]) * fraction_female
  matf[is.na(matf)] <- 0
  mid  <- ages + n / 2

  NRR <- sum(Lx * matf)
  if (NRR <= 0) stop("Net reproduction rate is zero; check 'asfr' and 'nLx'.")
  A <- sum(mid * Lx * matf) / NRR          # mean age of childbearing

  # Coale's iterative solution of Lotka's equation y(r) = 1.
  r <- log(NRR) / A
  it <- 0L
  converged <- FALSE
  repeat {
    it <- it + 1L
    y  <- sum(exp(-r * mid) * Lx * matf)
    if (abs(y - 1) < tol) { converged <- TRUE; break }
    if (it >= max_iter) break
    r <- r + (y - 1) / A
  }
  if (!converged) {
    warning(sprintf(
      paste0("Lotka's equation did not converge in %d iteration(s): ",
             "|y(r) - 1| = %.3g, still above tol = %.3g. The returned 'r' ",
             "(and the b, d and c(x) derived from it) are provisional; ",
             "raise 'max_iter' or relax 'tol'."),
      max_iter, abs(y - 1), tol), call. = FALSE)
  }

  denom <- sum(exp(-r * mid) * Lx)         # over ALL ages present
  b  <- 1 / denom
  d  <- b - r
  cx <- b * exp(-r * mid) * Lx
  mean_age <- sum(mid * cx) / sum(cx)

  tbl <- data.frame(age = ages, n = n, nLx = as.numeric(data[[nLx]]),
                    asfr = as.numeric(data[[asfr]]), cx = cx)
  out <- list(r = r, b = b, d = d, NRR = NRR, mean_age = mean_age,
              iterations = it, converged = converged, radix = radix, table = tbl)
  class(out) <- "dem_stable"
  if (graph) {
    out$plot <- .plot_series(tbl$age, tbl$cx, ylab = "Proportion c(x)",
                             title = "Stable age distribution", geom = "col")
  }
  out
}

# Interval widths from a vector of interval lower bounds (last width repeats).
.interval_widths <- function(ages) {
  d <- diff(ages)
  if (length(d) == 0) return(rep(1, length(ages)))
  c(d, d[length(d)])
}

#' @export
print.dem_reproduction <- function(x, ...) {
  cat("Reproduction measures\n")
  cat(strrep("-", 40), "\n", sep = "")
  print(round(x$table, 6), row.names = FALSE)
  cat(strrep("-", 40), "\n", sep = "")
  cat(sprintf("TFR = %.3f   GRR = %.3f   NRR = %.3f\n", x$TFR, x$GRR, x$NRR))
  if (!is.na(x$NRR))
    cat(sprintf("Mean age of childbearing = %.2f   approx. intrinsic r = %.5f\n",
                x$mean_age, x$r_approx))
  invisible(x)
}

#' @export
print.dem_stable <- function(x, ...) {
  cat("Stable population model\n")
  cat(strrep("-", 46), "\n", sep = "")
  print(round(x$table, 6), row.names = FALSE)
  cat(strrep("-", 46), "\n", sep = "")
  cat(sprintf("intrinsic r = %.5f   birth rate b = %.5f   death rate d = %.5f\n",
              x$r, x$b, x$d))
  cat(sprintf("NRR = %.3f   mean age = %.2f   (%d iterations)\n",
              x$NRR, x$mean_age, x$iterations))
  invisible(x)
}

#' @export
plot.dem_reproduction <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call reproduction(..., graph = TRUE).")
  x$plot
}

#' @export
plot.dem_stable <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call stable_population(..., graph = TRUE).")
  x$plot
}
