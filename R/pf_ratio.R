# Brass P/F ratio method for fertility estimation.
# Source: United Nations (1983), Manual X, Chapter II (Brass; Coale-Trussell
# interpolation coefficients, Table 7).

# Interpolation coefficients (Manual X, Table 7) for estimating the average
# parity equivalents F(i) from the cumulated period fertility schedule.
.pf_coef <- function(timing) {
  if (timing == "end_of_period") {          # births in 12 months, age at end
    list(a = c(2.531, 3.321, 3.265, 3.442, 3.518, 3.862, 3.828),
         b = c(-0.188, -0.754, -0.627, -0.563, -0.763, -2.481, 0.016),
         cc = c(0.0024, 0.0161, 0.0145, 0.0029, 0.0006, -0.0001, -0.0002))
  } else {                                  # births tabulated by age at delivery
    list(a = c(2.147, 2.838, 2.760, 2.949, 3.029, 3.419, 3.535),
         b = c(-0.244, -0.758, -0.594, -0.566, -0.823, -2.966, -0.007),
         cc = c(0.0034, 0.0162, 0.0133, 0.0025, 0.0006, -0.0001, -0.0002))
  }
}

#' Brass P/F Ratio Method for Fertility Estimation
#'
#' Estimates the level of fertility by the Brass P/F ratio method (United
#' Nations, Manual X, Chapter II). The age pattern of current period fertility
#' is assumed reliable but its level suspect; the level is corrected using the
#' average parities of younger women, which are assumed accurate.
#'
#' From the average parities \eqn{P(i) = CEB(i)/W(i)} and the period fertility
#' rates \eqn{f(i) = B(i)/W(i)}, the cumulated schedule
#' \eqn{\phi(i) = 5 \sum_{j \le i} f(j)} is formed and interpolated to the
#' parity-equivalent \eqn{F(i) = \phi(i-1) + a_i f(i) + b_i f(i+1) + c_i \phi(7)}
#' using the Coale-Trussell coefficients of Manual X Table 7. The adjustment
#' factor \eqn{K} is the average of the \eqn{P(i)/F(i)} ratios over the age
#' groups given by `k_ages` (20-34 by default, the most reliable), and the
#' adjusted total fertility rate is \eqn{K \times TFR}.
#'
#' @param data A data frame with the seven reproductive age groups (lower
#'   bounds 15, 20, ..., 45).
#' @param age Column name for the (lower bound of the) age group.
#' @param women Column name for the number of women in each age group.
#' @param ceb Column name for children ever born in each age group.
#' @param births Column name for births in the year preceding the survey.
#' @param births_timing Either `"end_of_period"` (births in a 12-month period
#'   classified by mother's age at the survey; the default) or `"at_delivery"`
#'   (classified by age at the birth), selecting the Table 7 coefficient set.
#' @param k_ages Age-group lower bounds whose P/F ratios are averaged to form
#'   the adjustment factor `K` (default `c(20, 25, 30)`).
#' @param graph Logical; if `TRUE` (default) a plot of the P/F ratios by age is
#'   attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_pf`: a list with the adjustment factor `K`,
#'   the unadjusted and adjusted total fertility rates `TFR` and `TFR_adjusted`,
#'   and a `table` of `P`, `f`, `phi`, `F`, the `PF` ratios, and the adjusted
#'   rates `f_adjusted`.
#'
#' @examples
#' # Manual X, Chapter II worked example: Bangladesh, 1974 (K = 1.500)
#' bd <- data.frame(
#'   age    = seq(15, 45, 5),
#'   women  = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
#'   ceb    = c(1160919, 4901382, 9085852, 9910256, 10384001, 9164329, 6905673),
#'   births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125)
#' )
#' pf_ratio(bd, age = "age", women = "women", ceb = "ceb", births = "births",
#'          graph = FALSE)
#'
#' @references
#' United Nations (1983). \emph{Manual X: Indirect Techniques for Demographic
#' Estimation}. New York: United Nations, Chapter II. Coale, A. J., & Trussell,
#' T. J. (1974). Model fertility schedules. \emph{Population Index}, 40(2),
#' 185-258.
#'
#' @seealso [chm_brass()], [reproduction()]
#' @export
pf_ratio <- function(data, age, women, ceb, births,
                     births_timing = c("end_of_period", "at_delivery"),
                     k_ages = c(20, 25, 30), graph = TRUE) {
  births_timing <- match.arg(births_timing)
  ages <- as.numeric(data[[age]])
  if (length(ages) != 7 || !all(ages == seq(15, 45, 5)))
    stop("Need exactly the seven age groups 15, 20, ..., 45.")
  W <- as.numeric(data[[women]])
  P <- as.numeric(data[[ceb]]) / W          # reported average parities
  f <- as.numeric(data[[births]]) / W       # period fertility rates

  phi <- 5 * cumsum(f)                       # cumulated fertility to end of group
  phi_prev <- c(0, phi[-7])                  # ... to start of group
  phi7 <- phi[7]

  co <- .pf_coef(births_timing)
  F <- numeric(7)
  for (i in 1:6) F[i] <- phi_prev[i] + co$a[i] * f[i] + co$b[i] * f[i + 1] + co$cc[i] * phi7
  F[7] <- phi_prev[7] + co$a[7] * f[7] + co$b[7] * f[6] + co$cc[7] * phi7   # b(7) applies to f(6)

  PF <- P / F
  # K is an average over the age groups named in k_ages. An unmatched k_ages
  # otherwise averages an empty vector and every result downstream is NaN.
  sel <- ages %in% k_ages
  if (!any(sel)) {
    stop("'k_ages' matched no age group. Give age-group lower bounds drawn ",
         "from ", paste(ages, collapse = ", "), "; got ",
         paste(k_ages, collapse = ", "), ".")
  }
  unmatched <- setdiff(k_ages, ages)
  if (length(unmatched) > 0) {
    warning("Ignoring 'k_ages' value(s) that are not age-group lower bounds: ",
            paste(unmatched, collapse = ", "), ".", call. = FALSE)
  }
  K  <- mean(PF[sel])
  TFR      <- 5 * sum(f)
  TFR_adj  <- K * TFR

  tab <- data.frame(age = ages, P = P, f = f, phi = phi, F = F,
                    PF = PF, f_adjusted = K * f)
  out <- list(K = K, TFR = TFR, TFR_adjusted = TFR_adj,
              births_timing = births_timing, table = tab)
  class(out) <- "dem_pf"
  if (graph) {
    out$plot <- .plot_series(ages, PF, ylab = "P/F ratio",
                             title = "Brass P/F ratios by age", geom = "col")
  }
  out
}

#' @export
print.dem_pf <- function(x, ...) {
  cat("Brass P/F ratio fertility estimation\n")
  cat(strrep("-", 52), "\n", sep = "")
  print(round(x$table, 4), row.names = FALSE)
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("Adjustment factor K = %.3f\n", x$K))
  cat(sprintf("TFR (reported) = %.2f   TFR (adjusted) = %.2f\n",
              x$TFR, x$TFR_adjusted))
  invisible(x)
}

#' @export
plot.dem_pf <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call pf_ratio(..., graph = TRUE).")
  x$plot
}
