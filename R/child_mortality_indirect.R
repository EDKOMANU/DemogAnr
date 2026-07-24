#' Indirect Child Mortality Estimation (Brass-Trussell method)
#'
#' Estimates child mortality indirectly from data on children ever born and
#' children surviving (or dead), classified by five-year age group of mother,
#' using the Trussell variant of the Brass method as set out in United Nations
#' \emph{Manual X} (1983). This is the standard technique for estimating child
#' mortality from census or survey data in populations that lack complete death
#' registration.
#'
#' @details
#' For each five-year age group of mother \eqn{i} (15-19, ..., 45-49) the
#' method computes the average parity \eqn{P(i) = CEB(i)/W(i)} and the
#' proportion of children dead \eqn{D(i) = CD(i)/CEB(i)}. Multipliers
#' \eqn{k(i)} that translate the reported proportions dead into life-table
#' probabilities of dying \eqn{q(x)} are obtained from
#' \deqn{k(i) = a(i) + b(i)\,P(1)/P(2) + c(i)\,P(2)/P(3),}
#' and \eqn{q(x) = k(i)\,D(i)}, where the age of the child \eqn{x} to which each
#' estimate refers is 1, 2, 3, 5, 10, 15 and 20 for \eqn{i = 1, \ldots, 7}. The
#' coefficients \eqn{a(i), b(i), c(i)} depend on the chosen Coale-Demeny model
#' life table family (`model`).
#'
#' Because a survey observes children born over a range of dates, each estimate
#' refers to a point in the past. The number of years before the survey to
#' which \eqn{q(x)} applies, the reference period \eqn{t(x)}, is estimated with
#' an analogous equation (Coale and Trussell 1977):
#' \deqn{t(x) = a'(i) + b'(i)\,P(1)/P(2) + c'(i)\,P(2)/P(3).}
#' If `survey_year` is supplied, the reference date of each estimate is returned
#' as `survey_year - t(x)`, so that a series of estimates from one survey traces
#' the trend in child mortality over the preceding two decades.
#'
#' The estimate derived from the youngest mothers (15-19) is well known to be
#' unreliable and is usually disregarded; see the references.
#'
#' @param data A data frame with seven rows, one per five-year age group of
#'   mother from 15-19 to 45-49, in ascending order.
#' @param women_col Column name for the number of women in each age group.
#' @param ceb_col Column name for the number of children ever born.
#' @param cd_col Column name for the number of children dead. Alternatively,
#'   supply `cs_col` (children surviving) instead.
#' @param cs_col Optional column name for the number of children surviving; if
#'   given, children dead are computed as `ceb - cs` and `cd_col` is ignored.
#' @param model Coale-Demeny model life table family: one of `"west"`
#'   (default), `"north"`, `"south"`, or `"east"`.
#' @param survey_year Optional numeric decimal year of the survey (e.g.
#'   `2010.4`). If given, a reference date is computed for each estimate.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @return A data frame (of class `dem_chm_indirect`) with one row per age
#'   group of mother, containing: `age_group`, index `i`, child age `x`, average
#'   parity `P`, proportion dead `D`, multiplier `k`, the probability of dying
#'   `qx` = \eqn{q(x)}, survivorship `lx` = \eqn{1 - q(x)}, the reference period
#'   `t` (years before the survey) and, if `survey_year` was given, the
#'   `ref_date`.
#'
#' @examples
#' # United Nations Manual X worked example: Panama, 1976 (both sexes)
#' panama <- data.frame(
#'   age_group = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
#'   women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
#'   ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
#'   cd    = c(40, 130, 312, 435, 636, 686, 689)
#' )
#' est <- chm_brass(panama, women_col = "women", ceb_col = "ceb",
#'                  cd_col = "cd", model = "west", survey_year = 1976.7)
#' est
#'
#' @references
#' United Nations (1983). \emph{Manual X: Indirect Techniques for Demographic Estimation}. Population Studies No. 81. New York: United Nations. (Chapter III; Tables 47 and 48.)
#'
#' Brass, W., & Coale, A. J. (1968). Methods of analysis and estimation. In W. Brass et al. (Eds.), \emph{The Demography of Tropical Africa} (pp. 88-139). Princeton: Princeton University Press.
#'
#' Coale, A. J., & Trussell, T. J. (1977). Estimating the time to which Brass estimates apply. \emph{Population Bulletin of the United Nations}, 10, 87-89.
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapter 11.)
#'
#' @export
chm_brass <- function(data, women_col, ceb_col, cd_col = NULL, cs_col = NULL,
                      model = c("west", "north", "south", "east"),
                      survey_year = NULL, verbose = FALSE) {
  model <- match.arg(model)

  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (nrow(data) != 7) {
    stop("'data' must have exactly 7 rows (age groups 15-19 to 45-49, ascending).")
  }
  if (!(women_col %in% names(data)) || !(ceb_col %in% names(data))) {
    stop("women_col and ceb_col must exist in 'data'.")
  }

  women <- as.numeric(data[[women_col]])
  ceb   <- as.numeric(data[[ceb_col]])

  if (!is.null(cs_col)) {
    if (!(cs_col %in% names(data))) stop("cs_col not found in 'data'.")
    cd <- ceb - as.numeric(data[[cs_col]])
  } else {
    if (is.null(cd_col) || !(cd_col %in% names(data))) {
      stop("Provide either cd_col (children dead) or cs_col (children surviving).")
    }
    cd <- as.numeric(data[[cd_col]])
  }

  if (any(women <= 0)) stop("Number of women must be positive in every age group.")
  if (any(ceb < 0) || any(cd < 0)) stop("Counts must be non-negative.")

  # Average parity and proportion dead
  P <- ceb / women
  D <- ifelse(ceb > 0, cd / ceb, 0)

  if (P[2] == 0 || P[3] == 0) stop("Average parities P(2) and P(3) must be positive.")
  P1P2 <- P[1] / P[2]
  P2P3 <- P[2] / P[3]

  # Coefficient tables (UN Manual X 1983, Tables 47 and 48), by CD model family.
  coef <- .chm_brass_coef(model)
  k <- coef$mult["a", ] + coef$mult["b", ] * P1P2 + coef$mult["c", ] * P2P3
  tref <- coef$time["a", ] + coef$time["b", ] * P1P2 + coef$time["c", ] * P2P3

  x <- c(1, 2, 3, 5, 10, 15, 20) # age of child to which q(x) refers
  qx <- as.numeric(k) * D

  if (verbose) {
    message(sprintf("chm_brass: model = %s, P(1)/P(2) = %.4f, P(2)/P(3) = %.4f",
                    model, P1P2, P2P3))
  }

  # Age-group labels: use an existing label column if present, else standard
  label_col <- intersect(c("age_group", "agegroup", "age", "Age"), names(data))
  age_labels <- if (length(label_col) > 0) {
    as.character(data[[label_col[1]]])
  } else {
    c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  }

  out <- data.frame(
    age_group = age_labels,
    i = 1:7,
    x = x,
    P = P,
    D = D,
    k = as.numeric(k),
    qx = qx,
    lx = 1 - qx,
    t = as.numeric(tref),
    stringsAsFactors = FALSE
  )
  if (!is.null(survey_year)) out$ref_date <- survey_year - out$t

  attr(out, "model") <- model
  attr(out, "survey_year") <- survey_year
  class(out) <- c("dem_chm_indirect", "data.frame")
  out
}

# Coale-Demeny Trussell coefficients (UN Manual X 1983).
# mult: multipliers k(i) = a + b*P1/P2 + c*P2/P3      (Table 47)
# time: reference period t(x) = a + b*P1/P2 + c*P2/P3 (Table 48)
.chm_brass_coef <- function(model) {
  m <- function(a, b, cc) rbind(a = a, b = b, c = cc)
  tbl <- list(
    north = list(
      mult = m(c(1.1119, 1.2390, 1.1884, 1.2046, 1.2586, 1.2240, 1.1772),
               c(-2.9287, -0.6865, 0.0421, 0.3037, 0.4236, 0.4222, 0.3486),
               c(0.8507, -0.2745, -0.5156, -0.5656, -0.5898, -0.5456, -0.4624)),
      time = m(c(1.0921, 1.3207, 1.5996, 2.0779, 2.7705, 4.1520, 6.9650),
               c(5.4732, 5.3751, 2.6268, -1.7908, -7.3403, -12.2448, -13.9160),
               c(-1.9672, 0.2133, 4.3701, 9.4126, 14.9352, 19.2349, 19.9542))
    ),
    south = list(
      mult = m(c(1.0819, 1.2846, 1.2223, 1.1905, 1.1911, 1.1564, 1.1307),
               c(-3.0005, -0.6181, 0.0851, 0.2631, 0.3152, 0.3017, 0.2596),
               c(0.8689, -0.3024, -0.4704, -0.4487, -0.4291, -0.3958, -0.3538)),
      time = m(c(1.0900, 1.3079, 1.5173, 1.9399, 2.6157, 4.0794, 7.1796),
               c(5.4443, 5.5568, 2.6755, -2.2739, -8.4819, -13.8308, -15.3880),
               c(-1.9721, 0.2021, 4.7471, 10.3876, 16.5153, 21.1866, 21.7892))
    ),
    east = list(
      mult = m(c(1.1461, 1.2231, 1.1593, 1.1404, 1.1540, 1.1336, 1.1201),
               c(-2.2536, -0.4301, 0.0581, 0.1991, 0.2511, 0.2556, 0.2362),
               c(0.6259, -0.2245, -0.3479, -0.3487, -0.3506, -0.3428, -0.3268)),
      time = m(c(1.0959, 1.2921, 1.5021, 1.9347, 2.6197, 4.1317, 7.3657),
               c(5.5864, 5.5897, 2.4692, -2.6419, -8.9693, -14.3550, -15.8083),
               c(-1.9949, 0.3631, 5.0927, 10.8533, 17.0981, 21.8247, 22.3005))
    ),
    west = list(
      mult = m(c(1.1415, 1.2563, 1.1851, 1.1720, 1.1865, 1.1746, 1.1639),
               c(-2.7070, -0.5381, 0.0633, 0.2341, 0.3080, 0.3314, 0.3190),
               c(0.7663, -0.2637, -0.4177, -0.4272, -0.4452, -0.4537, -0.4435)),
      time = m(c(1.0970, 1.3062, 1.5305, 1.9991, 2.7632, 4.3468, 7.5242),
               c(5.5628, 5.5677, 2.5528, -2.4261, -8.4065, -13.2436, -14.2013),
               c(-1.9956, 0.2962, 4.8962, 10.4282, 16.1787, 20.1990, 20.0162))
    )
  )
  tbl[[model]]
}

#' @export
print.dem_chm_indirect <- function(x, ...) {
  model <- attr(x, "model")
  cat(sprintf("Indirect child mortality estimates (Brass-Trussell, model %s)\n",
              toupper(model)))
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) round(v, 4))
  print(df, row.names = FALSE)
  # Headline: q(5) from mothers 30-34 (i=4) is a conventional under-5 indicator
  q5 <- x$qx[x$x == 5]
  if (length(q5) == 1) {
    cat(sprintf("\nImplied under-five mortality q(5): %.1f per 1,000",
                1000 * q5))
    if (!is.null(x$ref_date)) {
      cat(sprintf("  (reference date %.1f)", x$ref_date[x$x == 5]))
    }
    cat("\n")
  }
  invisible(x)
}
