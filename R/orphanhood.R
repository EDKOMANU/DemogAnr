# Estimation of adult female mortality from maternal orphanhood.
# Source: United Nations (1983), Manual X, Chapter IV, section B.
#
# The coefficients below are transcribed from Manual X tables 86 and 88. They
# were checked in two ways: every row of table 86 increases monotonically with
# the mean age at maternity, and the whole chain reproduces the published
# Bolivia 1975 worked example (tables 91 and 92) exactly -- all seven
# weighting factors to four decimal places and all seven survivorship ratios
# to three. That example is kept as a test.

# Table 86: weighting factors W(n) for converting the proportion of
# respondents with mother alive into female survivorship. Rows are the age of
# the respondent, columns the mean age of mothers at maternity, M = 22..30.
.orph_W <- rbind(
  `10` = c(0.420, 0.470, 0.517, 0.557, 0.596, 0.634, 0.674, 0.717, 0.758),
  `15` = c(0.418, 0.489, 0.556, 0.618, 0.678, 0.738, 0.800, 0.863, 0.924),
  `20` = c(0.404, 0.500, 0.590, 0.673, 0.756, 0.838, 0.921, 1.004, 1.085),
  `25` = c(0.366, 0.485, 0.598, 0.704, 0.809, 0.913, 1.016, 1.118, 1.218),
  `30` = c(0.303, 0.445, 0.580, 0.708, 0.834, 0.957, 1.080, 1.203, 1.323),
  `35` = c(0.241, 0.401, 0.554, 0.701, 0.844, 0.986, 1.128, 1.270, 1.412),
  `40` = c(0.125, 0.299, 0.467, 0.630, 0.791, 0.950, 1.111, 1.274, 1.442),
  `45` = c(0.007, 0.186, 0.361, 0.535, 0.708, 0.884, 1.063, 1.250, 1.447),
  `50` = c(-0.190, -0.017, 0.158, 0.334, 0.514, 0.699, 0.890, 1.095, 1.318),
  `55` = c(-0.368, -0.220, -0.059, 0.101, 0.270, 0.456, 0.645, 0.856, 1.083),
  `60` = c(-0.466, -0.352, -0.217, -0.084, 0.053, 0.220, 0.378, 0.579, 0.800)
)
.orph_M <- 22:30          # columns of table 86
.orph_ages <- seq(10, 60, 5)

# Table 88: the standard function Z(x) used to locate an estimate in time,
# for exact ages 26 to 75.
.orph_Z <- c(
  0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.091,
  0.092, 0.093, 0.095, 0.099, 0.104, 0.109, 0.115, 0.122, 0.130, 0.139,
  0.149, 0.160, 0.171, 0.182, 0.193, 0.205, 0.218, 0.231, 0.245, 0.259,
  0.274, 0.289, 0.305, 0.321, 0.338, 0.356, 0.374, 0.392, 0.411, 0.431,
  0.452, 0.473, 0.495, 0.518, 0.542, 0.568, 0.595, 0.622, 0.650, 0.678
)
.orph_Zage <- 26:75

# Linear interpolation in table 86 at a given mean age at maternity.
.orph_weights <- function(M) {
  if (M < min(.orph_M) || M > max(.orph_M)) {
    stop(sprintf(
      "The mean age at maternity M = %.2f is outside the tabulated range %d to %d.",
      M, min(.orph_M), max(.orph_M)))
  }
  stats::setNames(
    apply(.orph_W, 1, function(row) stats::approx(.orph_M, row, xout = M)$y),
    rownames(.orph_W))
}

# Z(x) is tabulated for exact ages 26 to 75. Beyond that there is no value to
# interpolate, so the time reference is simply not available -- which is why
# Manual X's own worked example stops at the 45-49 age group.
.orph_Zfun <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x >= min(.orph_Zage) & x <= max(.orph_Zage)
  out[ok] <- stats::approx(.orph_Zage, .orph_Z, xout = x[ok])$y
  out
}


#' Adult Female Mortality from Maternal Orphanhood
#'
#' Estimates adult female survivorship from the proportion of respondents whose
#' mother is still alive, by the Brass method as set out in United Nations
#' \emph{Manual X}, Chapter IV. It is the adult counterpart of [chm_brass()]:
#' where that estimates child mortality from a mother's report on her children,
#' this estimates the mother's own mortality from her children's report on her.
#'
#' @details
#' A respondent aged \eqn{n} was, by definition, born to a mother who was alive
#' at the birth. The proportion of respondents in an age group whose mother is
#' still alive is therefore a measure of how many women have survived from the
#' mean age at maternity to that age plus \eqn{n}. Weighting factors \eqn{W(n)}
#' from \emph{Manual X} table 86 convert those proportions into conditional
#' survivorship from exact age 25,
#' \deqn{\frac{l(25+n)}{l(25)} = W(n)\,S(n-5) + (1 - W(n))\,S(n),}
#' where \eqn{S(n)} is the proportion with mother alive in the five-year age
#' group beginning at \eqn{n}, and \eqn{W(n)} depends on the mean age of
#' mothers at maternity, \eqn{M}.
#'
#' Because a respondent's mother has been exposed to the risk of dying ever
#' since the birth, each estimate refers to a period in the past rather than to
#' the survey date. That period is
#' \deqn{t(n) = \frac{n\,(1 - u(n))}{2}, \qquad
#'   u(n) = 0.3333 \log({}_{10}S_{n-5}) + Z(M+n) + 0.0037\,(27 - M),}
#' with \eqn{{}_{10}S_{n-5}} the proportion with mother alive over the ten-year
#' age group \eqn{n-5} to \eqn{n+4} and \eqn{Z} the standard function of table
#' 88. Estimates from maternal orphanhood typically refer to between 8 and 15
#' years before the survey. \eqn{Z} is tabulated only to exact age 75, so where
#' \eqn{M + n} exceeds that -- the oldest respondents -- no time reference can
#' be formed and `t` is `NA`. \emph{Manual X}'s own worked example stops at the
#' 45-49 age group for this reason.
#'
#' The estimate from the youngest respondents is affected by the adoption
#' effect -- orphaned children are more likely to be reported as the children
#' of whoever raised them -- which biases survivorship upwards, and \emph{Manual
#' X} advises treating the youngest age groups with caution.
#'
#' @param data A data frame with one row per five-year age group of respondent,
#'   in ascending order, starting at age 15.
#' @param age Column name for the lower bound of the respondent's age group.
#' @param alive Column name for the number of respondents whose mother is
#'   alive.
#' @param dead Column name for the number whose mother is dead. Alternatively
#'   give `prop` directly.
#' @param prop Optional column name for the proportion with mother alive,
#'   \eqn{S(n)}, if the counts are not available. Time references then require
#'   `alive` and `dead` as well, since they weight the ten-year groups by the
#'   number of respondents; without counts the groups are weighted equally.
#' @param prop10 Optional column name for the proportion with mother alive over
#'   the ten-year age group \eqn{n-5} to \eqn{n+4}, \eqn{{}_{10}S_{n-5}}, which
#'   the time reference needs. Survey tabulations often publish it directly;
#'   when it is not given it is formed from the two five-year groups, weighted
#'   by the number of respondents where the counts are known.
#' @param M Mean age of mothers at maternity, between 22 and 30. This is the
#'   mean age of the fertility schedule, and is commonly computed from births
#'   in the year before the survey by age of mother.
#' @param survey_year Optional decimal year of the survey. If given, each
#'   estimate is dated as `survey_year - t`.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   survivorship estimates against the date they refer to is attached as
#'   `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_orphanhood`: a list with a `table` of, for
#'   each age group, the proportion with mother alive `S`, the weighting factor
#'   `W`, the conditional survivorship `lx_ratio` = \eqn{l(25+n)/l(25)}, the
#'   time reference `t` in years before the survey and, where `survey_year` was
#'   given, the `ref_date`; plus the `M` used.
#'
#' @examples
#' # United Nations Manual X, Chapter IV: Bolivia, 1975.
#' # The published estimates are in tables 91 and 92.
#' bolivia <- data.frame(
#'   age   = seq(15, 50, 5),
#'   alive = c(5540, 3995, 2886, 1910, 1272, 855, 541, 234),
#'   dead  = c(448, 541, 769, 852, 945, 1059, 985, 788),
#'   # the ten-year proportions as published in Manual X table 92
#'   s10   = c(NA, 0.9060, 0.8401, 0.7474, 0.6313, 0.5175, 0.3996, NA)
#' )
#' orphanhood(bolivia, age = "age", alive = "alive", dead = "dead",
#'            prop10 = "s10", M = 28.8, survey_year = 1975.5, graph = FALSE)
#'
#' @references
#' United Nations (1983). \emph{Manual X: Indirect Techniques for Demographic Estimation}. Population Studies No. 81. New York: United Nations. (Chapter IV, section B; tables 86 and 88.)
#'
#' Brass, W., & Hill, K. (1973). Estimating adult mortality from orphanhood. In \emph{International Population Conference, Liege, 1973} (Vol. 3, pp. 111-123). Liege: IUSSP.
#'
#' Hill, K., & Trussell, T. J. (1977). Further developments in indirect mortality estimation. \emph{Population Studies}, 31(2), 313-334. \doi{10.1080/00324728.1977.10410432}
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: IUSSP. \url{https://demographicestimation.iussp.org/}
#'
#' @seealso [chm_brass()] for the child mortality counterpart, and
#'   [model_lifetable()] for turning a survivorship estimate into a life table.
#' @export
orphanhood <- function(data, age, alive = NULL, dead = NULL, prop = NULL,
                       prop10 = NULL, M, survey_year = NULL, graph = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!(age %in% names(data))) stop("Column not found in 'data': ", age)
  if (!is.numeric(M) || length(M) != 1 || !is.finite(M)) {
    stop("'M', the mean age at maternity, must be a single number.")
  }

  ages <- as.numeric(data[[age]])
  if (is.unsorted(ages, strictly = TRUE)) {
    stop("Ages must be the strictly increasing lower bounds of the age groups.")
  }

  n_resp <- NULL
  if (!is.null(prop)) {
    if (!(prop %in% names(data))) stop("Column not found in 'data': ", prop)
    S <- as.numeric(data[[prop]])
    if (!is.null(alive) && !is.null(dead)) {
      n_resp <- as.numeric(data[[alive]]) + as.numeric(data[[dead]])
    }
  } else {
    if (is.null(alive) || is.null(dead)) {
      stop("Give either 'prop', or both 'alive' and 'dead'.")
    }
    for (nm in c(alive, dead)) {
      if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
    }
    a <- as.numeric(data[[alive]]); d <- as.numeric(data[[dead]])
    if (any(a < 0) || any(d < 0)) stop("Counts must be non-negative.")
    n_resp <- a + d
    if (any(n_resp <= 0)) stop("Every age group must have at least one respondent.")
    S <- a / n_resp
  }
  if (any(S < 0 | S > 1, na.rm = TRUE)) {
    stop("The proportion with mother alive must be between 0 and 1.")
  }

  W_all <- .orph_weights(M)
  # Estimates run from the second age group on, since each needs S(n-5)
  idx <- which(ages >= 20 & as.character(ages) %in% names(W_all))
  idx <- idx[idx >= 2]
  if (length(idx) == 0) {
    stop("No age group can be estimated. Ages 20 to 60 are needed, each with ",
         "the group below it present.")
  }

  n     <- ages[idx]
  W     <- W_all[as.character(n)]
  ratio <- W * S[idx - 1] + (1 - W) * S[idx]

  # Ten-year proportion covering n-5 to n+4. Taken as given when supplied,
  # otherwise formed from the two five-year groups, weighted by respondents
  # where the counts are known and equally where they are not.
  S10 <- if (!is.null(prop10)) {
    if (!(prop10 %in% names(data))) stop("Column not found in 'data': ", prop10)
    # NA is allowed: it simply leaves that group without a time reference,
    # which is already the case wherever Z(M + n) runs off the end of table 88.
    v <- as.numeric(data[[prop10]])[idx]
    if (any(v <= 0 | v > 1, na.rm = TRUE)) {
      stop("'prop10' must be between 0 and 1.")
    }
    v
  } else if (!is.null(n_resp)) {
    (S[idx - 1] * n_resp[idx - 1] + S[idx] * n_resp[idx]) /
      (n_resp[idx - 1] + n_resp[idx])
  } else {
    (S[idx - 1] + S[idx]) / 2
  }

  u <- 0.3333 * log(S10) + .orph_Zfun(M + n) + 0.0037 * (27 - M)
  t <- n * (1 - u) / 2

  tab <- data.frame(
    age = n, S_prev = S[idx - 1], S = S[idx], S10 = S10,
    W = as.numeric(W), lx_ratio = as.numeric(ratio),
    u = u, t = t, row.names = NULL
  )
  if (!is.null(survey_year)) tab$ref_date <- survey_year - t

  out <- list(table = tab, M = M, survey_year = survey_year)
  class(out) <- "dem_orphanhood"
  if (graph) out$plot <- .plot_orphanhood(tab, survey_year)
  out
}

#' @export
print.dem_orphanhood <- function(x, ...) {
  cat("Adult female survivorship from maternal orphanhood (Brass; Manual X)\n")
  cat(sprintf("  mean age at maternity M = %.2f\n", x$M))
  cat(strrep("-", 64), "\n", sep = "")
  tb <- x$table
  disp <- data.frame(
    age = tb$age,
    S = sprintf("%.4f", tb$S),
    W = sprintf("%.4f", tb$W),
    `l(25+n)/l(25)` = sprintf("%.4f", tb$lx_ratio),
    `years before` = sprintf("%.1f", tb$t),
    check.names = FALSE)
  if (!is.null(tb$ref_date)) disp[["refers to"]] <- sprintf("%.1f", tb$ref_date)
  print(disp, row.names = FALSE)
  cat(strrep("-", 64), "\n", sep = "")
  cat("Estimates from the youngest respondents are biased upwards by the\n")
  cat("adoption effect and are conventionally discounted.\n")
  invisible(x)
}

#' @export
plot.dem_orphanhood <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call orphanhood(..., graph = TRUE).")
  x$plot
}
