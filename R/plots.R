# Internal plotting helpers and the exported pyramid() function.
# All analysis plots share a light theme and a common two-colour palette so
# that a demographer gets a clean, publication-ready figure with no extra work.

.dem_male   <- "#4E79A7"
.dem_female <- "#E15759"

# A minimal, readable theme. base_family lets users pass a font (e.g.
# "Century Gothic") if it is installed; the default uses the device font.
.dem_theme <- function(base_family = "") {
  theme_minimal(base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# Anchor a non-negative value axis at 0 (no padding below), with a little
# headroom above.
.y_from_zero <- function() {
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

# Order an age vector for plotting: numeric ages by value, "0-4"-style labels
# by their leading number.
.order_age <- function(x) {
  if (is.numeric(x)) return(factor(x, levels = sort(unique(x))))
  xx <- as.character(x)
  key <- suppressWarnings(as.numeric(sub("^[^0-9-]*(-?[0-9]+).*", "\\1", xx)))
  lv <- unique(xx[order(key, xx)])
  factor(xx, levels = lv)
}

# A single value-by-age plot (line for rates, column optional).
.plot_series <- function(age, value, ylab, title, geom = "line") {
  d <- data.frame(age = .order_age(age), value = as.numeric(value))
  p <- ggplot(d, aes(x = age, y = value, group = 1))
  if (geom == "col") {
    p <- p + geom_col(fill = .dem_male)
  } else {
    p <- p + geom_line(colour = .dem_male) + geom_point(colour = .dem_male)
  }
  p + .y_from_zero() +
    labs(x = "Age group", y = ylab, title = title) + .dem_theme()
}

.plot_survival <- function(age, lx, title = "Survival curve l(x)") {
  d <- data.frame(age = as.numeric(age), lx = as.numeric(lx))
  ggplot(d, aes(x = age, y = lx)) +
    geom_line(colour = .dem_male) + geom_point(colour = .dem_male) +
    .y_from_zero() +
    labs(x = "Age", y = "Survivors l(x)", title = title) + .dem_theme()
}

.plot_projection <- function(df, probs = c(0.1, 0.9)) {
  d <- data.frame(year = df$year, subregion = as.character(df$subregion),
                  lo = df$lower, hi = df$upper, med = df$median)
  ggplot(d, aes(x = year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = .dem_male) +
    geom_line(aes(y = med), colour = .dem_male) +
    facet_wrap(~ subregion, scales = "free_y") +
    .y_from_zero() +
    labs(x = "Year", y = "Population",
         title = sprintf("Projected population (median and %g%% interval)",
                         100 * (probs[2] - probs[1]))) +
    .dem_theme()
}

# Calibrated regional series. When the calibration is by age group the ages are
# summed to a regional total for display; the band is then the sum of the
# age-specific bands, which is an approximation used for the picture only (the
# tabled quantiles are the authoritative ones).
.plot_calibration <- function(tbl, deterministic = FALSE) {
  d <- stats::aggregate(
    data.frame(med = tbl$calibrated, lo = tbl$lower, hi = tbl$upper),
    list(year = tbl$year, region = as.character(tbl$region)), sum)
  p <- ggplot(d, aes(x = year, y = med, group = region))
  if (!deterministic) {
    p <- p + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = .dem_male)
  }
  p + geom_line(colour = .dem_male) +
    facet_wrap(~ region, scales = "free_y") +
    .y_from_zero() +
    labs(x = "Year", y = "Population",
         title = "Regional population calibrated to the national total") +
    .dem_theme()
}

.plot_le_decomp <- function(tbl) {
  d <- data.frame(age = .order_age(tbl$age),
                  contribution = tbl$contribution,
                  sign = tbl$contribution >= 0)
  ggplot(d, aes(x = age, y = contribution, fill = sign)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = .dem_male, `FALSE` = .dem_female),
                      guide = "none") +
    labs(x = "Age group", y = "Contribution to the difference (years)",
         title = "Age decomposition of the life-expectancy gap") +
    .dem_theme()
}

.plot_myers <- function(tbl) {
  d <- data.frame(digit = factor(tbl$digit, levels = 0:9),
                  deviation = tbl$deviation)
  ggplot(d, aes(x = digit, y = deviation)) +
    geom_col(fill = .dem_male) + geom_hline(yintercept = 0) +
    labs(x = "Terminal digit", y = "Deviation from 10%",
         title = "Myers' blended index: digit preference") + .dem_theme()
}

.plot_whipple <- function(tbl) {
  d <- data.frame(digit = factor(tbl$terminal_digit, levels = 0:9),
                  percent = tbl$percent)
  ggplot(d, aes(x = digit, y = percent)) +
    geom_col(fill = .dem_male) +
    geom_hline(yintercept = 10, linetype = 2) +
    .y_from_zero() +
    labs(x = "Terminal digit", y = "Percent of population",
         title = "Distribution by terminal digit (10% expected)") + .dem_theme()
}

.plot_unasa <- function(tbl) {
  d <- rbind(
    data.frame(age = tbl$age, ratio = tbl$male_ratio, sex = "Male"),
    data.frame(age = tbl$age, ratio = tbl$female_ratio, sex = "Female")
  )
  d <- d[!is.na(d$ratio), ]
  d$age <- .order_age(d$age)
  ggplot(d, aes(x = age, y = ratio, colour = sex, group = sex)) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 100, linetype = 2) +
    scale_colour_manual(values = c(Male = .dem_male, Female = .dem_female)) +
    .y_from_zero() +
    labs(x = "Age group", y = "Age ratio (100 = smooth)", colour = NULL,
         title = "Age ratios by sex") + .dem_theme()
}


#' Population Pyramid
#'
#' Draws a population pyramid (back-to-back horizontal bar chart of population
#' by age group and sex) using \pkg{ggplot2}.
#'
#' @param data A data frame with one row per age group and sex.
#' @param age_col Column name for the age group.
#' @param sex_col Column name for sex.
#' @param count_col Column name for the population count.
#' @param male_label,female_label The values in `sex_col` denoting males and
#'   females (defaults "Male" and "Female"). Males are drawn to the left.
#' @param male_color,female_color Bar fill colours.
#' @param title Optional plot title.
#' @param base_family Font family passed to the theme; leave as `""` to use the
#'   graphics device default, or set it to an installed font (e.g. "Century
#'   Gothic").
#'
#' @return A \pkg{ggplot2} object, which can be printed, saved with
#'   [ggplot2::ggsave()], or further customised by adding layers.
#'
#' @examples
#' # Ghana 2021: single-year ages grouped into five-year age groups, by sex
#' data(ghpop2021)
#' g <- within(ghpop2021, grp <- ifelse(Age >= 80, "80+",
#'             paste0((Age %/% 5) * 5, "-", (Age %/% 5) * 5 + 4)))
#' pyr <- aggregate(Population ~ grp + Sex, data = g, FUN = sum)
#' pyramid(pyr, age_col = "grp", sex_col = "Sex", count_col = "Population",
#'         title = "Ghana, 2021")
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#'
#' @export
pyramid <- function(data, age_col, sex_col, count_col,
                    male_label = "Male", female_label = "Female",
                    male_color = .dem_male, female_color = .dem_female,
                    title = NULL, base_family = "") {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (nm in c(age_col, sex_col, count_col)) {
    if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
  }
  df <- data.frame(
    age = .order_age(data[[age_col]]),
    sex = as.character(data[[sex_col]]),
    count = as.numeric(data[[count_col]]),
    stringsAsFactors = FALSE
  )
  if (!all(c(male_label, female_label) %in% df$sex)) {
    stop("male_label and female_label must match values in sex_col.")
  }
  # Males to the left (negative), females to the right
  df$y <- ifelse(df$sex == male_label, -df$count, df$count)

  ggplot(df, aes(x = age, y = y, fill = sex)) +
    geom_col() +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = stats::setNames(c(male_color, female_color),
                                               c(male_label, female_label))) +
    coord_flip() +
    labs(x = "Age group", y = "Population", fill = NULL, title = title) +
    .dem_theme(base_family)
}


# Growth balance diagnostic: the fitted line should pass through points that
# lie straight. Curvature is the signal that the method's assumptions fail.
.plot_ggb <- function(tbl, intercept, slope, age_range) {
  d <- data.frame(x = tbl$death_rate, y = tbl$y,
                  used = ifelse(tbl$fitted, "fitted", "excluded"))
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  ggplot(d, aes(x = x, y = y)) +
    ggplot2::geom_abline(intercept = intercept, slope = slope,
                         colour = .dem_male, linewidth = 0.6) +
    geom_point(aes(colour = used), size = 2) +
    scale_colour_manual(values = c(fitted = .dem_male, excluded = "grey65")) +
    labs(x = "Observed death rate above x,  d(x+)",
         y = "Entry rate minus growth rate,  b(x+) - r(x+)",
         colour = NULL,
         title = "Generalized growth balance",
         subtitle = sprintf("fitted over exact ages %g-%g; completeness = %.1f%%",
                            age_range[1], age_range[2], 100 / slope)) +
    .dem_theme()
}

# Extinct generations diagnostic: the age-specific completeness ratios should
# be flat. A trend means the estimate is not trustworthy.
.plot_seg <- function(tbl, completeness, age_range) {
  d <- data.frame(age = tbl$age, ratio = tbl$ratio,
                  used = ifelse(tbl$fitted, "averaged", "excluded"))
  d <- d[is.finite(d$ratio), ]
  ggplot(d, aes(x = age, y = ratio)) +
    geom_hline(yintercept = completeness, linetype = 2, colour = .dem_male) +
    geom_line(colour = "grey70") +
    geom_point(aes(colour = used), size = 2) +
    scale_colour_manual(values = c(averaged = .dem_male, excluded = "grey65")) +
    labs(x = "Exact age x", y = "Estimated completeness at x",
         colour = NULL,
         title = "Synthetic extinct generations",
         subtitle = sprintf("averaged over ages %g-%g; completeness = %.1f%%",
                            age_range[1], age_range[2], 100 * completeness)) +
    .dem_theme()
}
