#' Whipple's Index of Age Heaping
#'
#' Computes Whipple's index, a measure of preference for ages ending in the
#' digits 0 and 5, from single-year-of-age population counts.
#'
#' @details
#' Over the standard age range 23 to 62, Whipple's index is
#' \deqn{WI = 100 \times \frac{\sum_{a \in \{25,30,\dots,60\}} P_a}
#'   {\frac{1}{5}\sum_{a=23}^{62} P_a},}
#' i.e. the population at ages ending in 0 or 5 relative to one fifth of the
#' total over the range. A value of 100 indicates no preference; 500 means every
#' age was reported as a multiple of 5. The United Nations quality bands are:
#' under 105 highly accurate, 105-109.9 fairly accurate, 110-124.9 approximate,
#' 125-174.9 rough, and 175 or more very rough.
#'
#' @param data A data frame with single years of age.
#' @param age_col Column name for single year of age.
#' @param pop_col Column name for the population count.
#' @param lower,upper Age range over which to evaluate the index (defaults 23
#'   and 62, the UN standard).
#'
#' @return An object of class `dem_whipple`: a list with the index value, its
#'   quality band, and the numerator/denominator used.
#'
#' @examples
#' # Data with strong heaping on multiples of five
#' ages <- 20:65
#' pop <- ifelse(ages %% 5 == 0, 5000, 1000)
#' d <- data.frame(age = ages, pop = pop)
#' whipple(d, "age", "pop")
#'
#' @references
#' Shryock, H. S., & Siegel, J. S. (1976). \emph{The Methods and Materials of Demography}. New York: Academic Press.
#'
#' United Nations (1955). \emph{Manual II: Methods of Appraisal of Quality of Basic Data for Population Estimates}. New York: United Nations.
#'
#' @export
whipple <- function(data, age_col, pop_col, lower = 23, upper = 62) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!(age_col %in% names(data)) || !(pop_col %in% names(data))) {
    stop("age_col and pop_col must exist in 'data'.")
  }
  age <- as.numeric(data[[age_col]])
  pop <- as.numeric(data[[pop_col]])

  in_range <- age >= lower & age <= upper
  if (sum(in_range) == 0) stop("No ages fall within [lower, upper].")

  # Ages ending in 0 or 5, kept away from the range boundaries (25..60 for the
  # standard 23..62 range)
  mult5 <- age %% 5 == 0 & age >= (lower + 2) & age <= (upper - 2)
  numerator <- sum(pop[mult5])
  denominator <- sum(pop[in_range])

  wi <- 100 * numerator / (denominator / 5)

  band <- if (wi < 105) "highly accurate" else
    if (wi < 110) "fairly accurate" else
    if (wi < 125) "approximate" else
    if (wi < 175) "rough" else "very rough"

  out <- list(index = wi, quality = band, numerator = numerator,
              denominator = denominator, range = c(lower, upper))
  class(out) <- "dem_whipple"
  out
}

#' @export
print.dem_whipple <- function(x, ...) {
  cat(sprintf("Whipple's index (ages %d-%d): %.1f\n",
              x$range[1], x$range[2], x$index))
  cat(sprintf("  Data quality: %s\n", x$quality))
  invisible(x)
}


#' Myers' Blended Index of Age Heaping
#'
#' Computes Myers' blended index, which measures preference for or avoidance of
#' every terminal digit 0-9, from single-year-of-age population counts.
#'
#' @details
#' In the absence of digit preference, the population whose age ends in each
#' digit should be one tenth of the total. Because a population normally
#' declines with age, Myers "blends" the tabulation to remove that trend: each
#' age \eqn{a} in the range \eqn{[lower, upper]} receives weight 10, ages just
#' below the range (\eqn{lower} to \eqn{lower+8}) receive weights 1 to 9, and
#' ages just above (\eqn{upper+1} to \eqn{upper+9}) receive weights 9 to 1
#' (Rodriguez 2015; equivalent to Myers 1940). The blended counts by terminal
#' digit are expressed as percentages \eqn{p_d}, and the index is
#' \deqn{M = \tfrac{1}{2}\sum_{d=0}^{9} |p_d - 10|,}
#' ranging from 0 (no heaping) to 90 (all ages at one digit). Note that the data
#' should extend up to age \eqn{upper + 9} for the blending weights to apply.
#'
#' @param data A data frame with single years of age.
#' @param age_col Column name for single year of age.
#' @param pop_col Column name for the population count.
#' @param lower,upper The initial tabulation range (defaults 10 and 89). The
#'   blend extends nine years beyond `upper`.
#'
#' @return An object of class `dem_myers`: a list with the index value and a
#'   `table` of blended counts and percentages by terminal digit.
#'
#' @examples
#' # Uniform single-year population: no heaping
#' d <- data.frame(age = 0:99, pop = rep(1000, 100))
#' myers(d, "age", "pop", lower = 10, upper = 89)
#'
#' @references
#' Myers, R. J. (1940). Errors and bias in the reporting of ages in census data. \emph{Transactions of the Actuarial Society of America}, 41(2), 395-415.
#'
#' Rodriguez, G. (2015). \emph{Demographic Methods} (course notes). Princeton University. \url{https://grodri.github.io/demography/}
#'
#' @export
myers <- function(data, age_col, pop_col, lower = 10, upper = 89) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!(age_col %in% names(data)) || !(pop_col %in% names(data))) {
    stop("age_col and pop_col must exist in 'data'.")
  }
  age <- as.numeric(data[[age_col]])
  pop <- as.numeric(data[[pop_col]])
  ord <- order(age); age <- age[ord]; pop <- pop[ord]

  # Myers blending weights (Rodriguez 2015):
  # ramp 1..9 over [lower, lower+8], 10 over [lower+9, upper], 9..1 over
  # [upper+1, upper+9], and 0 elsewhere.
  w <- ifelse(age < lower, 0,
       ifelse(age < lower + 9, age - lower + 1,
       ifelse(age <= upper, 10,
       ifelse(age <= upper + 9, upper + 10 - age, 0))))

  digit <- age %% 10
  blended <- tapply(pop * w, digit, sum)
  blended <- blended[as.character(0:9)]
  blended[is.na(blended)] <- 0
  pct <- 100 * blended / sum(blended)
  index <- sum(abs(pct - 10)) / 2

  tbl <- data.frame(digit = 0:9, blended = as.numeric(blended),
                    percent = as.numeric(pct),
                    deviation = as.numeric(pct - 10))
  out <- list(index = index, table = tbl, range = c(lower, upper))
  class(out) <- "dem_myers"
  out
}

#' @export
print.dem_myers <- function(x, ...) {
  cat(sprintf("Myers' blended index (ages %d-%d): %.2f\n",
              x$range[1], x$range[2], x$index))
  cat("  (0 = no digit preference, 90 = all ages at one digit)\n")
  invisible(x)
}


#' Sex Ratios by Age Group
#'
#' Computes the sex ratio (males per 100 females) for each age group.
#'
#' @param data A data frame with one row per age group.
#' @param age_col Column name for the age group.
#' @param male_col,female_col Column names for the male and female counts.
#'
#' @return A data frame with the age group, male and female counts, and the sex
#'   ratio.
#'
#' @examples
#' d <- data.frame(age = c("0-4","5-9","10-14"),
#'                 m = c(5100, 4800, 4500), f = c(5000, 4700, 4600))
#' sex_ratio(d, "age", "m", "f")
#'
#' @references
#' Shryock, H. S., & Siegel, J. S. (1976). \emph{The Methods and Materials of Demography}. New York: Academic Press.
#'
#' @export
sex_ratio <- function(data, age_col, male_col, female_col) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  m <- as.numeric(data[[male_col]]); f <- as.numeric(data[[female_col]])
  data.frame(age = data[[age_col]], males = m, females = f,
             sex_ratio = 100 * m / f)
}


#' Age Ratios by Age Group
#'
#' Computes United Nations age ratios, \eqn{100 \times 2 P_x / (P_{x-n} +
#' P_{x+n})}, which measure how far each age group departs from the average of
#' its neighbours. A ratio near 100 indicates a smooth age distribution.
#'
#' @param data A data frame with one row per (equally spaced) age group.
#' @param age_col Column name for the age group lower bound.
#' @param pop_col Column name for the population count.
#' @param open_ended Logical; if `TRUE` (default), the last age group is treated
#'   as open-ended, so no age ratio is computed for it or for the group below it
#'   (whose upper neighbour would be the open interval).
#'
#' @return A data frame with the age group, population, age ratio, and its
#'   deviation from 100 (`NA` where the ratio is not defined).
#'
#' @examples
#' d <- data.frame(age = seq(0, 20, 5), pop = c(1000, 900, 820, 760, 700))
#' age_ratio(d, "age", "pop", open_ended = FALSE)
#'
#' @references
#' United Nations (1955). \emph{Manual II: Methods of Appraisal of Quality of Basic Data for Population Estimates}. New York: United Nations.
#'
#' @export
age_ratio <- function(data, age_col, pop_col, open_ended = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  age <- data[[age_col]]
  pop <- as.numeric(data[[pop_col]])
  k <- length(pop)
  ratio <- rep(NA_real_, k)
  last_valid <- if (open_ended) k - 2 else k - 1
  if (last_valid >= 2) {
    for (i in 2:last_valid) {
      ratio[i] <- 100 * 2 * pop[i] / (pop[i - 1] + pop[i + 1])
    }
  }
  data.frame(age = age, pop = pop, age_ratio = ratio, deviation = ratio - 100)
}


#' United Nations Age-Sex Accuracy Index
#'
#' Computes the United Nations joint age-sex accuracy index, a summary measure
#' of the internal consistency of an age-sex distribution built from the sex
#' ratio score and the male and female age ratio scores.
#'
#' @details
#' Let the sex ratio score (SRS) be the mean absolute difference between
#' successive sex ratios, and the age ratio scores (ARSM, ARSF) be the mean
#' absolute deviation of the age ratios from 100 for males and females
#' respectively. The joint index is
#' \deqn{UN = 3 \times SRS + ARSM + ARSF.}
#' Sex ratios use the closed age groups; age ratios exclude the first group, the
#' open-ended group, and (when `open_ended = TRUE`) the group adjacent to it.
#' United Nations guidance regards values under 20 as accurate, 20-40 as
#' inaccurate, and over 40 as highly inaccurate.
#'
#' @param data A data frame with one row per five-year age group.
#' @param age_col Column name for the age group.
#' @param male_col,female_col Column names for the male and female counts.
#' @param open_ended Logical; if `TRUE` (default), the last age group is treated
#'   as open-ended.
#'
#' @return An object of class `dem_unasa`: a list with the joint index and its
#'   components (`SRS`, `ARSM`, `ARSF`), the quality band, and the underlying
#'   `sex_ratios` and `age_ratios` tables.
#'
#' @examples
#' # United Nations / Kpedekpo worked example: Ghana, 1960 census
#' ghana1960 <- data.frame(
#'   age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
#'           "40-44","45-49","50-54","55-59","60-64","65+"),
#'   males = c(643041, 515520, 357831, 275542, 268336, 278601, 242515,
#'             198231, 168937, 122756, 96775, 59307, 63467, 113185),
#'   females = c(654258, 503070, 323460, 265534, 322576, 306329, 245883,
#'               179182, 145572, 95590, 81715, 48412, 54572, 100392)
#' )
#' un_age_sex_accuracy(ghana1960, "age", "males", "females")
#'
#' @references
#' United Nations (1955). \emph{Manual II: Methods of Appraisal of Quality of Basic Data for Population Estimates}. New York: United Nations.
#'
#' Kpedekpo, G. M. K. (1982). \emph{Essentials of Demographic Analysis for Africa}. London: Heinemann. (Chapter 3.)
#'
#' @export
un_age_sex_accuracy <- function(data, age_col, male_col, female_col,
                                open_ended = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  m <- as.numeric(data[[male_col]]); f <- as.numeric(data[[female_col]])
  k <- length(m)

  # Sex ratios over the closed age groups
  closed <- if (open_ended) 1:(k - 1) else 1:k
  sr <- 100 * m[closed] / f[closed]
  srs <- mean(abs(diff(sr)))

  # Age ratios (males and females), excluding first, open, and pre-open groups
  arm <- age_ratio(data.frame(a = data[[age_col]], p = m), "a", "p", open_ended)
  arf <- age_ratio(data.frame(a = data[[age_col]], p = f), "a", "p", open_ended)
  arsm <- mean(abs(arm$deviation), na.rm = TRUE)
  arsf <- mean(abs(arf$deviation), na.rm = TRUE)

  joint <- 3 * srs + arsm + arsf
  band <- if (joint < 20) "accurate" else
    if (joint <= 40) "inaccurate" else "highly inaccurate"

  out <- list(
    index = joint, SRS = srs, ARSM = arsm, ARSF = arsf, quality = band,
    sex_ratios = data.frame(age = data[[age_col]][closed], sex_ratio = sr),
    age_ratios = data.frame(age = data[[age_col]], male_ratio = arm$age_ratio,
                            female_ratio = arf$age_ratio)
  )
  class(out) <- "dem_unasa"
  out
}

#' @export
print.dem_unasa <- function(x, ...) {
  cat("United Nations age-sex accuracy index\n")
  cat(sprintf("  Sex ratio score (SRS):        %.2f\n", x$SRS))
  cat(sprintf("  Age ratio score, males:       %.2f\n", x$ARSM))
  cat(sprintf("  Age ratio score, females:     %.2f\n", x$ARSF))
  cat(sprintf("  Joint index (3*SRS+ARSM+ARSF): %.2f\n", x$index))
  cat(sprintf("  Assessment: %s\n", x$quality))
  invisible(x)
}
