# Smoothing of age distributions distorted by age misreporting.
#
# The five-year methods are those of the U.S. Census Bureau's Population
# Analysis System workbook AGESMTH (Carrier-Farrag, Karup-King-Newton, Arriaga,
# Arriaga's strong smoothing, and the United Nations formula); the single-year
# weighted moving average is that of the NewPAS workbook POP1SM.

# --- Shared input handling ----------------------------------------------------

# Lower bounds of the age groups, from numbers or from labels such as "0-4",
# "15-19", "65+" or "80 and over".
.age_lower <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  lab <- trimws(as.character(x))
  lo <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", lab)))
  if (anyNA(lo)) {
    stop("Could not read an age from the label '", lab[is.na(lo)][1], "'.")
  }
  lo
}

# Labels for printing: "10-14" for a closed group of width n, "80+" for the
# open-ended one.
.age_labels <- function(lower, width, open_ended) {
  lab <- paste0(lower, "-", lower + width - 1)
  if (open_ended) lab[length(lab)] <- paste0(lower[length(lower)], "+")
  lab
}

# Validate an age distribution and split it into its closed groups and the
# optional open-ended group. `width` is 5 for grouped data and 1 for single
# years of age.
.age_dist <- function(data, age_col, pop_col, open_ended, width) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (nm in c(age_col, pop_col)) {
    if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
  }
  if (!is.logical(open_ended) || length(open_ended) != 1 || is.na(open_ended)) {
    stop("'open_ended' must be TRUE or FALSE.")
  }
  lower <- .age_lower(data[[age_col]])
  pop <- as.numeric(data[[pop_col]])
  if (anyNA(lower)) stop("Ages must not be missing.")
  if (anyNA(pop)) stop("Population counts must not be missing.")
  if (any(pop < 0)) stop("Population counts must be non-negative.")
  if (is.unsorted(lower, strictly = TRUE)) {
    stop("Age groups must be in strictly increasing order of age.")
  }
  step <- diff(lower)
  if (any(step != width)) {
    bad <- which(step != width)[1]
    unit <- if (width == 1) "single years of age" else
      paste0("contiguous ", width, "-year age groups")
    stop("The data must be ", unit, ", but the groups starting at ",
         lower[bad], " and ", lower[bad + 1], " are ", step[bad],
         " years apart.")
  }
  k <- length(pop)
  n_closed <- if (open_ended) k - 1 else k
  list(lower = lower, pop = pop, k = k, n_closed = n_closed,
       closed = seq_len(n_closed), open_ended = open_ended,
       label = .age_labels(lower, width, open_ended))
}

# Which closed groups lie wholly inside [lower, upper]?
.in_range <- function(age, width, lower, upper) {
  lo <- if (is.null(lower)) -Inf else lower
  hi <- if (is.null(upper)) Inf else upper
  if (!is.numeric(lo) || !is.numeric(hi) || length(lo) != 1 ||
      length(hi) != 1 || lo > hi) {
    stop("'lower' and 'upper' must be single ages with lower <= upper.")
  }
  age >= lo & age + width - 1 <= hi
}

# Mean absolute deviation of the age ratios from 100, over the groups where
# both series have a ratio. This is the age ratio score of the United Nations
# accuracy index, used here to show how much irregularity was removed.
.age_ratio_score <- function(r) mean(abs(r - 100), na.rm = TRUE)


#' Smooth a Five-Year Age Distribution
#'
#' Corrects a population distribution in five-year age groups for the net
#' effect of age misreporting, using the classical smoothing formulas of the
#' U.S. Census Bureau's Population Analysis System (workbook AGESMTH).
#'
#' @details
#' Age misreporting shifts people from their true age group into a
#' neighbouring one, which gives the age distribution a saw-tooth pattern.
#' Every method replaces each group by a weighted combination of the group and
#' its neighbours. They differ in how hard they smooth and in whether they keep
#' the reported total. Writing \eqn{{}_{10}P_x} for the population aged
#' \eqn{x} to \eqn{x+9}:
#'
#' \describe{
#'   \item{`"carrier_farrag"`}{Carrier and Farrag (1959). The second half of
#'     each ten-year group is
#'     \deqn{{}_5P_{x+5} = \frac{{}_{10}P_x}{1 + ({}_{10}P_{x-10} /
#'     {}_{10}P_{x+10})^{1/4}},}
#'     and the first half is the remainder. It leaves unchanged any age
#'     distribution that declines geometrically, at a constant ratio from one
#'     group to the next.}
#'   \item{`"kkn"`}{Karup-King-Newton. The first half of each ten-year group is
#'     \deqn{{}_5P_x = \tfrac{1}{2}\,{}_{10}P_x + \tfrac{1}{16}({}_{10}P_{x-10}
#'     - {}_{10}P_{x+10}),}
#'     which leaves any quadratic age distribution unchanged.}
#'   \item{`"arriaga"`}{Arriaga (1968), the default. A second-degree polynomial
#'     is passed through the midpoints of three consecutive ten-year groups and
#'     integrated over the upper half of the middle one:
#'     \deqn{{}_5P_{x+5} = (-{}_{10}P_{x-10} + 11\,{}_{10}P_x +
#'     2\,{}_{10}P_{x+10}) / 24.}
#'     Unlike the two formulas above it also has versions for the youngest and
#'     oldest ten-year groups, \eqn{(8\,{}_{10}P_x + 5\,{}_{10}P_{x+10} -
#'     {}_{10}P_{x+20})/24} and its mirror image, so it can split every
#'     ten-year group. It leaves a linear age distribution unchanged.}
#'   \item{`"strong"`}{Arriaga's strong smoothing. The ten-year groups are
#'     first smoothed themselves,
#'     \deqn{{}_{10}P'_x = ({}_{10}P_{x-10} + 2\,{}_{10}P_x + {}_{10}P_{x+10})
#'     / 4,}
#'     then rescaled to the reported total over the ages smoothed, and finally
#'     split into five-year groups with the Arriaga formula. The first and last
#'     ten-year groups are left as reported.}
#'   \item{`"united_nations"`}{United Nations (1956). Works on the five-year
#'     groups directly,
#'     \deqn{{}_5P'_x = (-{}_5P_{x-10} + 4\,{}_5P_{x-5} + 10\,{}_5P_x +
#'     4\,{}_5P_{x+5} - {}_5P_{x+10}) / 16,}
#'     so the two youngest and two oldest groups are left as reported. It
#'     leaves a cubic age distribution unchanged.}
#' }
#'
#' Carrier-Farrag, Karup-King-Newton and Arriaga only redistribute people
#' between the two halves of each ten-year group, so every ten-year total, and
#' with it the overall total, is kept exactly. The strong and United Nations
#' formulas move people between ten-year groups and change the total slightly;
#' with `preserve_total = TRUE` (the default) the smoothed groups are rescaled
#' to the reported total over the same ages, as AGESMTH does for strong
#' smoothing.
#'
#' Smoothing is not always the right correction. It assumes that people were
#' reported at the wrong age, not that they were missed: where an age group was
#' under-enumerated, smoothing spreads the deficit to its neighbours. Genuine
#' irregularities, such as the cohorts born during a war or a famine, are
#' smoothed away along with the errors. And smoothing each sex separately does
#' not make the sex ratios plausible. Inspect the data first (for example with
#' [un_age_sex_accuracy()]), and use `lower` and `upper` to confine smoothing
#' to the ages that need it. AGESMTH itself applies the light formulas from
#' age 10 upwards, leaving ages under 10 as reported; `lower = 10` follows that
#' convention.
#'
#' The light formulas only move people between the two halves of a ten-year
#' group, so they cannot undo misreporting that carries people from one
#' ten-year group into another, such as young women reported at 20-29 who are
#' really 10-19; strong smoothing can, at the price of flattening genuine
#' features along with the errors.
#'
#' @param data A data frame with one row per five-year age group, in increasing
#'   order of age.
#' @param age_col Column name for the age group: either its lower bound (0, 5,
#'   10, ...) or a label such as `"0-4"` or `"80+"`.
#' @param pop_col Column name for the population count.
#' @param method One of `"arriaga"` (default), `"carrier_farrag"`, `"kkn"`,
#'   `"strong"` or `"united_nations"`; see Details.
#' @param lower,upper Optional age range to smooth. Only the age groups lying
#'   wholly within `[lower, upper]` are replaced (for the ten-year methods, both
#'   halves of a ten-year group must lie within it); the groups outside are
#'   returned as reported, though they still serve as neighbours in the
#'   formulas.
#' @param offset Where the ten-year groups start: `0` (default) pairs the
#'   five-year groups as 0-9, 10-19, ..., as AGESMTH does; `5` pairs them as
#'   5-14, 15-24, ..., as the NewPAS workbook ARRSM5 does. With heaping on ages
#'   ending in 0, people who should be just below such an age are counted at
#'   it, in the first half of the next ten-year group; pairing from 5 keeps
#'   both sides of the heaped age in the same ten-year group, so the split can
#'   move them back. Not used by the United Nations formula.
#' @param open_ended Logical; if `TRUE` (default) the last age group is
#'   open-ended. It is returned unchanged and is not used in the formulas.
#' @param preserve_total Logical; if `TRUE` (default), the strong and United
#'   Nations results are rescaled so that the smoothed age groups add up to the
#'   reported total over the same ages. The other methods keep the total
#'   exactly whatever this is set to.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   reported and smoothed distributions is attached as `$plot`; retrieve it
#'   with `plot()`.
#'
#' @return An object of class `dem_age_smooth`: a list with a `table` giving,
#'   for each age group, the `reported` and `smoothed` populations, whether the
#'   group was `changed`, and the age ratios before and after
#'   (`ratio_reported`, `ratio_smoothed`); the age ratio scores before and
#'   after (`score_reported`, `score_smoothed`, the mean absolute deviation of
#'   the age ratios from 100); the reported and smoothed totals; and the
#'   method settings.
#'
#' @examples
#' # Ghana, 1960 census (Kpedekpo 1982), a distribution with marked heaping
#' ghana1960 <- data.frame(
#'   age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
#'           "40-44","45-49","50-54","55-59","60-64","65+"),
#'   females = c(654258, 503070, 323460, 265534, 322576, 306329, 245883,
#'               179182, 145572, 95590, 81715, 48412, 54572, 100392)
#' )
#' s <- smooth_ages(ghana1960, "age", "females", method = "arriaga")
#' s
#'
#' # The methods side by side: the age ratio score falls further the harder
#' # the method smooths
#' sapply(c("carrier_farrag", "kkn", "arriaga", "united_nations", "strong"),
#'        function(m) smooth_ages(ghana1960, "age", "females", method = m,
#'                                graph = FALSE)$score_smoothed)
#'
#' # Leave the young ages, where the shape may be genuine, as reported
#' smooth_ages(ghana1960, "age", "females", lower = 10, upper = 59,
#'             graph = FALSE)
#'
#' @references
#' Arriaga, E. E., Johnson, P. D., & Jamison, E. (1994). \emph{Population Analysis with Microcomputers, Volume I: Presentation of Techniques}. Washington, DC: U.S. Bureau of the Census. (Chapter II: age smoothing.)
#'
#' Arriaga, E. E. (1968). \emph{New Life Tables for Latin American Populations in the Nineteenth and Twentieth Centuries}. Population Monograph Series No. 3. Berkeley: University of California.
#'
#' Carrier, N. H., & Farrag, A. M. (1959). The reduction of errors in census populations for statistically underdeveloped countries. \emph{Population Studies}, 12(3), 240-285.
#'
#' United Nations (1956). \emph{Manual III: Methods for Population Projections by Sex and Age}. Population Studies No. 25. New York: United Nations.
#'
#' U.S. Census Bureau. \emph{Population Analysis System (PAS) Software}: workbooks AGESMTH and ARRSM5. Washington, DC. \url{https://www.census.gov/data/software/pas.html}
#'
#' @seealso [smooth_single_ages()] for single years of age, [graduate()] for
#'   splitting the smoothed groups into single years, and [age_ratio()] and
#'   [un_age_sex_accuracy()] for deciding whether smoothing is needed.
#' @export
smooth_ages <- function(data, age_col, pop_col,
                        method = c("arriaga", "carrier_farrag", "kkn",
                                   "strong", "united_nations"),
                        lower = NULL, upper = NULL, offset = 0,
                        open_ended = TRUE, preserve_total = TRUE,
                        graph = TRUE) {
  method <- match.arg(method)
  if (!is.numeric(offset) || length(offset) != 1 || !(offset %in% c(0, 5))) {
    stop("'offset' must be 0 or 5.")
  }
  d <- .age_dist(data, age_col, pop_col, open_ended, width = 5)
  if (any(d$lower %% 5 != 0)) {
    stop("Five-year age groups must start at ages that are multiples of 5.")
  }
  P <- d$pop[d$closed]
  age <- d$lower[d$closed]
  inside <- .in_range(age, 5, lower, upper)

  new <- if (method == "united_nations") {
    .smooth_un(P, inside, preserve_total)
  } else {
    .smooth_pairs(P, age, method, offset, inside, preserve_total)
  }
  changed <- !is.na(new) & inside
  if (!any(changed)) {
    stop("No age group could be smoothed. The method needs neighbouring ",
         "groups on both sides, so check that 'lower' and 'upper' leave ",
         "room for them and that there are enough age groups.")
  }
  smoothed <- P
  smoothed[changed] <- new[changed]
  if (any(smoothed < 0)) {
    warning("Smoothing produced negative counts at age(s) ",
            paste(age[smoothed < 0], collapse = ", "),
            ". The formulas extrapolate badly where a group is very small ",
            "next to its neighbours; consider leaving those ages out with ",
            "'upper'.", call. = FALSE)
  }

  full_rep <- d$pop
  full_sm  <- c(smoothed, if (open_ended) d$pop[d$k])
  r_rep <- age_ratio(data.frame(a = d$lower, p = full_rep), "a", "p",
                     open_ended)$age_ratio
  r_sm  <- age_ratio(data.frame(a = d$lower, p = full_sm), "a", "p",
                     open_ended)$age_ratio
  both <- !is.na(r_rep) & !is.na(r_sm)

  tab <- data.frame(age = d$lower, interval = d$label,
                    reported = full_rep, smoothed = full_sm,
                    changed = c(changed, if (open_ended) FALSE),
                    ratio_reported = r_rep, ratio_smoothed = r_sm,
                    stringsAsFactors = FALSE)

  out <- list(table = tab, method = method,
              offset = if (method == "united_nations") NA else offset,
              score_reported = .age_ratio_score(r_rep[both]),
              score_smoothed = .age_ratio_score(r_sm[both]),
              total_reported = sum(full_rep), total_smoothed = sum(full_sm),
              preserve_total = preserve_total, open_ended = open_ended)
  class(out) <- "dem_age_smooth"
  if (graph) {
    out$plot <- .plot_age_smooth(tab[d$closed, ], .smooth_name(method))
  }
  out
}

.smooth_name <- function(method) {
  switch(method,
         carrier_farrag = "Carrier-Farrag",
         kkn = "Karup-King-Newton",
         arriaga = "Arriaga",
         strong = "Arriaga strong smoothing",
         united_nations = "United Nations")
}

# The four methods that work on ten-year groups. Returns the smoothed value of
# every closed five-year group, NA where the method cannot reach it.
.smooth_pairs <- function(P, age, method, offset, inside, preserve_total) {
  k <- length(P)
  # Each ten-year group is a five-year group starting at an age ending in
  # `offset` together with the group after it. The groups are contiguous and
  # start at multiples of 5, so these pairs never overlap; a group left over
  # at either end is simply not paired.
  first <- which(age %% 10 == offset & seq_len(k) < k)
  J <- length(first)
  if (J < 3) {
    stop("The ", .smooth_name(method), " method needs at least three ",
         "complete ten-year age groups (six five-year groups).")
  }
  second <- first + 1
  T10 <- P[first] + P[second]
  pair_in <- inside[first] & inside[second]

  # first half of each ten-year group; NA where the method does not apply
  h1 <- rep(NA_real_, J)
  mid <- 2:(J - 1)

  if (method == "carrier_farrag") {
    if (any(T10 <= 0)) {
      stop("The Carrier-Farrag method needs every ten-year group to be ",
           "positive.")
    }
    h1[mid] <- T10[mid] - T10[mid] / (1 + (T10[mid - 1] / T10[mid + 1])^0.25)
  } else if (method == "kkn") {
    h1[mid] <- T10[mid] / 2 + (T10[mid - 1] - T10[mid + 1]) / 16
  } else if (method == "arriaga") {
    h1 <- .arriaga_first_half(T10)
  } else if (method == "strong") {
    S <- mid[pair_in[mid]]
    if (length(S) == 0) return(rep(NA_real_, k))
    U <- T10
    U[S] <- (T10[S - 1] + 2 * T10[S] + T10[S + 1]) / 4
    if (preserve_total) U[S] <- U[S] * sum(T10[S]) / sum(U[S])
    split <- .arriaga_first_half(U)
    h1[S] <- split[S]
    T10[S] <- U[S]
  }

  out <- rep(NA_real_, k)
  ok <- !is.na(h1) & pair_in
  out[first[ok]] <- h1[ok]
  out[second[ok]] <- T10[ok] - h1[ok]
  out
}

# Arriaga's split of a series of ten-year totals: the first five-year half of
# each, using the end formulas for the youngest and oldest groups.
.arriaga_first_half <- function(T10) {
  J <- length(T10)
  mid <- 2:(J - 1)
  h1 <- numeric(J)
  h1[mid] <- T10[mid] - (-T10[mid - 1] + 11 * T10[mid] + 2 * T10[mid + 1]) / 24
  h1[1] <- T10[1] - (8 * T10[1] + 5 * T10[2] - T10[3]) / 24
  h1[J] <- (-T10[J - 2] + 5 * T10[J - 1] + 8 * T10[J]) / 24
  h1
}

.smooth_un <- function(P, inside, preserve_total) {
  k <- length(P)
  if (k < 5) {
    stop("The United Nations formula needs at least five closed age groups.")
  }
  out <- rep(NA_real_, k)
  i <- 3:(k - 2)
  out[i] <- (-P[i - 2] + 4 * P[i - 1] + 10 * P[i] + 4 * P[i + 1] -
               P[i + 2]) / 16
  use <- !is.na(out) & inside
  if (preserve_total && any(use)) {
    out[use] <- out[use] * sum(P[use]) / sum(out[use])
  }
  out
}

#' @export
print.dem_age_smooth <- function(x, ...) {
  cat(sprintf("Smoothed age distribution: %s\n", .smooth_name(x$method)))
  if (!is.na(x$offset)) {
    cat(sprintf("  ten-year groups start at ages ending in %d\n", x$offset))
  }
  cat(strrep("-", 66), "\n", sep = "")
  tb <- x$table
  fmt_n <- function(v) format(round(v), big.mark = ",", trim = TRUE,
                              scientific = FALSE)
  fmt_r <- function(v) ifelse(is.na(v), "", sprintf("%.1f", v))
  disp <- data.frame(
    age = tb$interval,
    reported = fmt_n(tb$reported),
    smoothed = fmt_n(tb$smoothed),
    `age ratio, reported` = fmt_r(tb$ratio_reported),
    `age ratio, smoothed` = fmt_r(tb$ratio_smoothed),
    check.names = FALSE, stringsAsFactors = FALSE)
  print(disp, row.names = FALSE)
  cat(strrep("-", 66), "\n", sep = "")
  if (any(!tb$changed)) {
    cat("  Left as reported:", paste(tb$interval[!tb$changed], collapse = ", "),
        "\n")
  }
  change <- x$total_smoothed - x$total_reported
  cat(sprintf("  Total: reported %s, smoothed %s",
              fmt_n(x$total_reported), fmt_n(x$total_smoothed)))
  if (abs(change) > 1e-6 * x$total_reported) {
    cat(sprintf(" (%+.2f%%)", 100 * change / x$total_reported))
  }
  cat("\n")
  cat(sprintf("  Age ratio score (mean |ratio - 100|): %.2f reported, %.2f smoothed\n",
              x$score_reported, x$score_smoothed))
  invisible(x)
}

#' @export
plot.dem_age_smooth <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call smooth_ages(..., graph = TRUE).")
  x$plot
}


#' Smooth Single Years of Age with a Weighted Moving Average
#'
#' Removes digit preference from a population by single years of age by
#' replacing each age with a weighted average of the ages around it, as the
#' U.S. Census Bureau's NewPAS workbook POP1SM does.
#'
#' @details
#' The smoothed population at age \eqn{x} is
#' \deqn{S(x) = \frac{\sum_{i=-h}^{h} w(i)\,P(x+i)}{\sum_{i=-h}^{h} w(i)},}
#' with symmetric weights \eqn{w}. The default is the POP1SM triangular kernel,
#' \eqn{w(i) \propto h + 1 - |i|} with \eqn{h = 4}: weights of 4, 8, 12, 16,
#' 20, 16, 12, 8 and 4 per cent. A wider kernel, or flatter weights, smooths
#' more.
#'
#' Near the ends of the age range there are not \eqn{h} ages on both sides, so
#' the window narrows: an age \eqn{d} years from the end uses half-width
#' \eqn{\min(h, d)}, and the first and last ages are left as reported. A named
#' kernel keeps its shape as it narrows (the triangular kernel at half-width 3
#' is 6.25, 12.5, 18.75, 25, 18.75, 12.5 and 6.25 per cent), whereas a numeric
#' `weights` vector is cut to the shorter window and rescaled. These are the
#' rules behind the trimmed weights that POP1SM tabulates for its default and
#' alternative weights respectively, and both are reproduced exactly.
#'
#' Heaping on particular digits is a feature of single years of age, and much
#' of it disappears once the data are grouped. If only five-year groups are
#' needed, [smooth_ages()] is usually the better tool.
#'
#' @param data A data frame with one row per single year of age, in increasing
#'   order.
#' @param age_col Column name for the age.
#' @param pop_col Column name for the population count.
#' @param weights Either `"triangular"` (default), the POP1SM kernel with
#'   half-width `span`, or a numeric vector of symmetric, non-negative weights
#'   of odd length, centred on the age being smoothed. POP1SM's alternative
#'   weights, for instance, are `c(4.5, rep(9, 9), 4.5)`. Weights need not sum
#'   to one.
#' @param span Half-width of the triangular kernel, in years (default 4, which
#'   averages over nine ages). Ignored when `weights` is numeric.
#' @param lower,upper Optional age range to smooth. Ages outside it are
#'   returned as reported, though they still enter the averages for the ages
#'   inside it.
#' @param open_ended Logical; if `TRUE` (default) the last row is an
#'   open-ended age group. It is returned unchanged and is not used in the
#'   averages.
#' @param preserve_total Logical; if `TRUE` (default) the smoothed ages are
#'   rescaled to the reported total over the same ages. A moving average keeps
#'   the total only approximately, because the windows narrow at the ends of
#'   the range. Set to `FALSE` for the moving average alone.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   reported and smoothed populations by age is attached as `$plot`; retrieve
#'   it with `plot()`.
#'
#' @return An object of class `dem_age_smooth_single`: a list with a `table` of
#'   the `reported` and `smoothed` population and whether each age was
#'   `changed` (with `open_ended = TRUE` its last row is the open-ended group,
#'   returned as reported); the interior `weights` used, as proportions; the reported and
#'   smoothed totals; and, where the data cover ages 23 to 62, Whipple's index
#'   before and after (`whipple_reported`, `whipple_smoothed`).
#'
#' @examples
#' # Ghana 2021, males: heaping on ages ending in 0 and 5
#' data(ghpop2021)
#' male <- subset(ghpop2021, Sex == "Male")
#' s <- smooth_single_ages(male, "Age", "Population")
#' s
#'
#' # POP1SM's flatter alternative weights smooth harder
#' smooth_single_ages(male, "Age", "Population",
#'                    weights = c(4.5, rep(9, 9), 4.5), graph = FALSE)
#'
#' @references
#' U.S. Census Bureau. \emph{Population Analysis System (PAS) Software}: workbook POP1SM. Washington, DC. \url{https://www.census.gov/data/software/pas.html}
#'
#' Arriaga, E. E., Johnson, P. D., & Jamison, E. (1994). \emph{Population Analysis with Microcomputers, Volume I: Presentation of Techniques}. Washington, DC: U.S. Bureau of the Census.
#'
#' @seealso [whipple()] and [myers()] for measuring digit preference, and
#'   [smooth_ages()] for five-year age groups.
#' @export
smooth_single_ages <- function(data, age_col, pop_col, weights = "triangular",
                               span = 4, lower = NULL, upper = NULL,
                               open_ended = TRUE, preserve_total = TRUE,
                               graph = TRUE) {
  d <- .age_dist(data, age_col, pop_col, open_ended, width = 1)
  P <- d$pop[d$closed]
  age <- d$lower[d$closed]
  k <- length(P)

  kern <- .ma_kernel(weights, span)
  H <- kern$H
  if (k < 3) stop("At least three single years of age are needed.")

  sm <- numeric(k)
  for (i in seq_len(k)) {
    hw <- min(H, i - 1, k - i)
    w <- kern$at(hw)
    sm[i] <- sum(w * P[(i - hw):(i + hw)]) / sum(w)
  }
  inside <- .in_range(age, 1, lower, upper)
  if (!any(inside)) stop("No age lies within [lower, upper].")
  if (preserve_total && sum(sm[inside]) > 0) {
    sm[inside] <- sm[inside] * sum(P[inside]) / sum(sm[inside])
  }
  smoothed <- P
  smoothed[inside] <- sm[inside]

  full_sm <- c(smoothed, if (open_ended) d$pop[d$k])
  tab <- data.frame(age = d$lower, reported = d$pop, smoothed = full_sm,
                    changed = c(inside, if (open_ended) FALSE))

  out <- list(table = tab, weights = kern$at(H) / sum(kern$at(H)),
              kernel = if (is.character(weights)) weights else "user-supplied",
              total_reported = sum(d$pop), total_smoothed = sum(full_sm),
              preserve_total = preserve_total)
  if (min(age) <= 23 && max(age) >= 62) {
    wr <- data.frame(a = age, p = P)
    ws <- data.frame(a = age, p = smoothed)
    out$whipple_reported <- whipple(wr, "a", "p", graph = FALSE)$index
    out$whipple_smoothed <- whipple(ws, "a", "p", graph = FALSE)$index
  }
  class(out) <- "dem_age_smooth_single"
  if (graph) out$plot <- .plot_age_smooth_single(tab[seq_len(k), ])
  out
}

# The moving-average kernel and how it narrows. Returns the interior
# half-width H and a function giving the weights at any half-width up to H.
.ma_kernel <- function(weights, span) {
  if (is.character(weights)) {
    if (!identical(weights, "triangular")) {
      stop("'weights' must be \"triangular\" or a numeric vector of weights.")
    }
    if (!is.numeric(span) || length(span) != 1 || span < 1 || span != round(span)) {
      stop("'span' must be a whole number of years, at least 1.")
    }
    # the triangle keeps its shape as the window narrows
    return(list(H = span, at = function(hw) hw + 1 - abs(-hw:hw)))
  }
  w <- as.numeric(weights)
  if (length(w) < 3 || length(w) %% 2 != 1) {
    stop("A numeric 'weights' vector must have an odd length of at least 3, ",
         "centred on the age being smoothed.")
  }
  if (anyNA(w) || any(w < 0)) stop("Weights must be non-negative numbers.")
  if (!isTRUE(all.equal(w, rev(w)))) stop("Weights must be symmetric.")
  H <- (length(w) - 1) / 2
  if (w[H + 1] <= 0) stop("The weight on the age being smoothed must be positive.")
  # user weights are cut to the narrower window and rescaled
  list(H = H, at = function(hw) w[(H + 1 - hw):(H + 1 + hw)])
}

#' @export
print.dem_age_smooth_single <- function(x, ...) {
  cat("Smoothed single-year age distribution (weighted moving average)\n")
  w <- x$weights
  cat(sprintf("  %s weights over %d ages, in per cent:\n",
              if (x$kernel == "triangular") "Triangular" else "User-supplied",
              length(w)))
  cat("   ", paste(sprintf("%.2f", 100 * w), collapse = "  "), "\n")
  tb <- x$table
  ch <- tb$age[tb$changed]
  cat(sprintf("  Ages smoothed: %g-%g\n", min(ch), max(ch)))
  fmt_n <- function(v) format(round(v), big.mark = ",", trim = TRUE,
                              scientific = FALSE)
  change <- x$total_smoothed - x$total_reported
  cat(sprintf("  Total: reported %s, smoothed %s", fmt_n(x$total_reported),
              fmt_n(x$total_smoothed)))
  if (abs(change) > 1e-6 * x$total_reported) {
    cat(sprintf(" (%+.2f%%)", 100 * change / x$total_reported))
  }
  cat("\n")
  if (!is.null(x$whipple_reported)) {
    cat(sprintf("  Whipple's index (23-62): %.1f reported, %.1f smoothed\n",
                x$whipple_reported, x$whipple_smoothed))
  }
  cat(sprintf("  (%d ages; the full series is in $table)\n", nrow(tb)))
  invisible(x)
}

#' @export
plot.dem_age_smooth_single <- function(x, ...) {
  if (is.null(x$plot)) {
    stop("No plot available; call smooth_single_ages(..., graph = TRUE).")
  }
  x$plot
}
