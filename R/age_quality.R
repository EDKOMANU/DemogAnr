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
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} bar chart of the
#'   distribution by terminal digit is attached to the result and shown on print.
#'
#' @return An object of class `dem_whipple`: a list with the index value, its
#'   quality band, the numerator/denominator used, and a `table` giving the
#'   population and share by terminal digit over the range.
#'
#' @examples
#' # Ghana 2021 single-year ages, males
#' data(ghpop2021)
#' male <- subset(ghpop2021, Sex == "Male")
#' whipple(male, age_col = "Age", pop_col = "Population")
#'
#' @references
#' Shryock, H. S., & Siegel, J. S. (1976). \emph{The Methods and Materials of Demography}. New York: Academic Press.
#'
#' United Nations (1955). \emph{Manual II: Methods of Appraisal of Quality of Basic Data for Population Estimates}. New York: United Nations.
#'
#' @export
whipple <- function(data, age_col, pop_col, lower = 23, upper = 62, graph = TRUE) {
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
  expected <- denominator / 5
  wi <- 100 * numerator / expected

  band <- if (wi < 105) "highly accurate" else
    if (wi < 110) "fairly accurate" else
    if (wi < 125) "approximate" else
    if (wi < 175) "rough" else "very rough"

  # Distribution by terminal digit over the range (expected share 10% each)
  d <- age[in_range] %% 10
  by_digit <- tapply(pop[in_range], d, sum)
  by_digit <- by_digit[as.character(0:9)]
  by_digit[is.na(by_digit)] <- 0
  tbl <- data.frame(
    terminal_digit = 0:9,
    population = as.numeric(by_digit),
    percent = 100 * as.numeric(by_digit) / denominator
  )

  out <- list(index = wi, quality = band, numerator = numerator,
              denominator = denominator, expected = expected,
              range = c(lower, upper), table = tbl)
  if (graph) out$plot <- .plot_whipple(tbl)
  class(out) <- "dem_whipple"
  out
}

#' @export
plot.dem_whipple <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call whipple(..., graph = TRUE).")
  x$plot
}

#' @export
print.dem_whipple <- function(x, ...) {
  cat(sprintf("Whipple's index of age heaping (ages %d-%d)\n",
              x$range[1], x$range[2]))
  cat("\nPopulation by terminal digit over the range:\n")
  tb <- x$table
  tb$population <- format(round(tb$population), big.mark = ",", trim = TRUE, scientific = FALSE)
  tb$percent <- sprintf("%.2f", x$table$percent)
  print(tb, row.names = FALSE)
  cat(sprintf("\n  Population at ages ending 0 or 5 (%d-%d): %s\n",
              x$range[1] + 2, x$range[2] - 2,
              format(round(x$numerator), big.mark = ",", trim = TRUE, scientific = FALSE)))
  cat(sprintf("  One fifth of the total (%d-%d):          %s\n",
              x$range[1], x$range[2],
              format(round(x$expected), big.mark = ",", trim = TRUE, scientific = FALSE)))
  cat(sprintf("  Whipple's index:                         %.1f (%s)\n",
              x$index, x$quality))
  if (!is.null(x$plot)) print(x$plot)
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
#' ranging from 0 (no heaping) to 90 (all ages at one digit). The data should
#' extend up to age \eqn{upper + 9} for the blending weights to apply.
#'
#' @param data A data frame with single years of age.
#' @param age_col Column name for single year of age.
#' @param pop_col Column name for the population count.
#' @param lower,upper The initial tabulation range (defaults 10 and 89). The
#'   blend extends nine years beyond `upper`.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} bar chart of the
#'   deviation from 10% by terminal digit is attached to the result and shown
#'   on print.
#'
#' @return An object of class `dem_myers`: a list with the index value and a
#'   `table` giving, for each terminal digit, the reported and blended counts,
#'   the blended percentage, and its deviation from 10.
#'
#' @examples
#' # Ghana 2021 single-year ages, males (blend over 10-69)
#' data(ghpop2021)
#' male <- subset(ghpop2021, Sex == "Male")
#' myers(male, age_col = "Age", pop_col = "Population", lower = 10, upper = 69)
#'
#' @references
#' Myers, R. J. (1940). Errors and bias in the reporting of ages in census data. \emph{Transactions of the Actuarial Society of America}, 41(2), 395-415.
#'
#' Rodriguez, G. (2015). \emph{Demographic Methods} (course notes). Princeton University. \url{https://grodri.github.io/demography/}
#'
#' @export
myers <- function(data, age_col, pop_col, lower = 10, upper = 89, graph = TRUE) {
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
  reported <- tapply(pop * (w > 0), digit, sum)   # counts actually used
  blended <- tapply(pop * w, digit, sum)
  reported <- reported[as.character(0:9)]; reported[is.na(reported)] <- 0
  blended <- blended[as.character(0:9)];   blended[is.na(blended)] <- 0
  pct <- 100 * blended / sum(blended)
  index <- sum(abs(pct - 10)) / 2

  tbl <- data.frame(digit = 0:9,
                    reported = as.numeric(reported),
                    blended = as.numeric(blended),
                    percent = as.numeric(pct),
                    deviation = as.numeric(pct - 10))
  out <- list(index = index, table = tbl, range = c(lower, upper))
  if (graph) out$plot <- .plot_myers(tbl)
  class(out) <- "dem_myers"
  out
}

#' @export
plot.dem_myers <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call myers(..., graph = TRUE).")
  x$plot
}

#' @export
print.dem_myers <- function(x, ...) {
  cat(sprintf("Myers' blended index of age heaping (ages %d-%d)\n",
              x$range[1], x$range[2]))
  cat("\nBlended distribution by terminal digit:\n")
  tb <- x$table
  tb$reported <- format(round(tb$reported), big.mark = ",", trim = TRUE, scientific = FALSE)
  tb$blended <- format(round(tb$blended), big.mark = ",", trim = TRUE, scientific = FALSE)
  tb$percent <- sprintf("%.2f", x$table$percent)
  tb$deviation <- sprintf("%+.2f", x$table$deviation)
  print(tb, row.names = FALSE)
  cat(sprintf("\n  Myers' index (half the sum of |deviations|): %.2f\n", x$index))
  cat("  (0 = no digit preference, 90 = all ages at one digit)\n")
  if (!is.null(x$plot)) print(x$plot)
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
#' data(ghpop2021)
#' # collapse single years into five-year groups, by sex
#' g <- within(ghpop2021, grp <- ifelse(Age >= 80, 80, (Age %/% 5) * 5))
#' wide <- tapply(g$Population, list(g$grp, g$Sex), sum)
#' sr_data <- data.frame(age = rownames(wide),
#'                       m = wide[, "Male"], f = wide[, "Female"])
#' head(sex_ratio(sr_data, "age", "m", "f"))
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
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the male
#'   and female age ratios by age is attached to the result and shown on print.
#'
#' @return An object of class `dem_unasa`: a list with the joint index, its
#'   components (`SRS`, `ARSM`, `ARSF`), the quality band, and a full `table`
#'   with the male and female counts, their age ratios and deviations, the sex
#'   ratio and its successive differences.
#'
#' @examples
#' # Ghana 2021, five-year age groups by sex (from the single-year data)
#' data(ghpop2021)
#' g <- within(ghpop2021, grp <- ifelse(Age >= 80, 80, (Age %/% 5) * 5))
#' wide <- tapply(g$Population, list(g$grp, g$Sex), sum)
#' five <- data.frame(age = as.integer(rownames(wide)),
#'                    males = wide[, "Male"], females = wide[, "Female"])
#' un_age_sex_accuracy(five, "age", "males", "females")
#'
#' # Reproduces the published 1960 Ghana census worked example (Kpedekpo 1982)
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
                                open_ended = TRUE, graph = TRUE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  m <- as.numeric(data[[male_col]]); f <- as.numeric(data[[female_col]])
  k <- length(m)

  # Age ratios (males and females)
  arm <- age_ratio(data.frame(a = data[[age_col]], p = m), "a", "p", open_ended)
  arf <- age_ratio(data.frame(a = data[[age_col]], p = f), "a", "p", open_ended)
  arsm <- mean(abs(arm$deviation), na.rm = TRUE)
  arsf <- mean(abs(arf$deviation), na.rm = TRUE)

  # Sex ratios over the closed age groups, and their successive differences
  sr <- 100 * m / f
  sr_used <- if (open_ended) c(rep(TRUE, k - 1), FALSE) else rep(TRUE, k)
  sr[!sr_used] <- NA
  sr_diff <- rep(NA_real_, k)
  used_idx <- which(sr_used)
  if (length(used_idx) >= 2) {
    for (j in 2:length(used_idx)) {
      i <- used_idx[j]
      sr_diff[i] <- abs(sr[i] - sr[used_idx[j - 1]])
    }
  }
  srs <- mean(sr_diff, na.rm = TRUE)

  joint <- 3 * srs + arsm + arsf
  band <- if (joint < 20) "accurate" else
    if (joint <= 40) "inaccurate" else "highly inaccurate"

  tbl <- data.frame(
    age = data[[age_col]],
    males = m, male_ratio = arm$age_ratio, male_dev = arm$deviation,
    females = f, female_ratio = arf$age_ratio, female_dev = arf$deviation,
    sex_ratio = sr, sr_diff = sr_diff
  )

  out <- list(index = joint, SRS = srs, ARSM = arsm, ARSF = arsf,
              quality = band, table = tbl)
  if (graph) out$plot <- .plot_unasa(tbl)
  class(out) <- "dem_unasa"
  out
}

#' @export
plot.dem_unasa <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call un_age_sex_accuracy(..., graph = TRUE).")
  x$plot
}

#' @export
print.dem_unasa <- function(x, ...) {
  cat("United Nations age-sex accuracy index\n\n")
  tb <- x$table
  fmt_n <- function(v) format(round(v), big.mark = ",", trim = TRUE, scientific = FALSE)
  fmt_r <- function(v) ifelse(is.na(v), "", sprintf("%.1f", v))
  disp <- data.frame(
    age = as.character(tb$age),
    males = fmt_n(tb$males),
    m_ratio = fmt_r(tb$male_ratio),
    m_dev = ifelse(is.na(tb$male_dev), "", sprintf("%+.1f", tb$male_dev)),
    females = fmt_n(tb$females),
    f_ratio = fmt_r(tb$female_ratio),
    f_dev = ifelse(is.na(tb$female_dev), "", sprintf("%+.1f", tb$female_dev)),
    sex_ratio = fmt_r(tb$sex_ratio),
    sr_diff = fmt_r(tb$sr_diff),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  cat(sprintf("\n  Sex ratio score (SRS):           %.2f\n", x$SRS))
  cat(sprintf("  Age ratio score, males (ARSM):   %.2f\n", x$ARSM))
  cat(sprintf("  Age ratio score, females (ARSF): %.2f\n", x$ARSF))
  cat(sprintf("  Joint index (3*SRS + ARSM + ARSF): %.2f (%s)\n",
              x$index, x$quality))
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}
