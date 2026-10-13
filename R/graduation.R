# Graduation of five-year age groups into single years of age with
# osculatory interpolation multipliers.
#
# Each method is a set of coefficient panels. A head panel (10 x 5) gives the
# single ages of the first two groups from the first five groups; a middle
# panel (5 x 5) gives the single ages of a group from it and the two groups on
# either side; and the tail panel is the head panel turned round, rows and
# columns reversed. Karup-King works the same way with panels three groups
# wide, and one group at each end.
#
# The tables are the standard four-decimal multipliers (Siegel and Swanson
# 2004, Appendix C; the Beers ordinary end panels are those of the U.S. Census
# Bureau's PAS workbooks). The Beers tables were transcribed from the DemoTools
# package (Riffe et al.), the Sprague table independently, and all of them are
# checked in the test suite against the identities that define them: every
# row sums to 1/5, each group's population is handed out exactly once, and,
# for all but Beers' modified formula, the single ages of a group add back to
# exactly that group.

.grad_panels <- list(
  sprague = list(
    head = matrix(c(
       0.3616, -0.2768,  0.1488, -0.0336,  0.0000,
       0.2640, -0.0960,  0.0400, -0.0080,  0.0000,
       0.1840,  0.0400, -0.0320,  0.0080,  0.0000,
       0.1200,  0.1360, -0.0720,  0.0160,  0.0000,
       0.0704,  0.1968, -0.0848,  0.0176,  0.0000,
       0.0336,  0.2272, -0.0752,  0.0144,  0.0000,
       0.0080,  0.2320, -0.0480,  0.0080,  0.0000,
      -0.0080,  0.2160, -0.0080,  0.0000,  0.0000,
      -0.0160,  0.1840,  0.0400, -0.0080,  0.0000,
      -0.0176,  0.1408,  0.0912, -0.0144,  0.0000), 10, byrow = TRUE),
    mid = matrix(c(
      -0.0128,  0.0848,  0.1504, -0.0240,  0.0016,
      -0.0016,  0.0144,  0.2224, -0.0416,  0.0064,
       0.0064, -0.0336,  0.2544, -0.0336,  0.0064,
       0.0064, -0.0416,  0.2224,  0.0144, -0.0016,
       0.0016, -0.0240,  0.1504,  0.0848, -0.0128), 5, byrow = TRUE)),

  beers_ordinary = list(
    head = matrix(c(
       0.3333, -0.1636, -0.0210,  0.0796, -0.0283,
       0.2595, -0.0780,  0.0130,  0.0100, -0.0045,
       0.1924,  0.0064,  0.0184, -0.0256,  0.0084,
       0.1329,  0.0844,  0.0054, -0.0356,  0.0129,
       0.0819,  0.1508, -0.0158, -0.0284,  0.0115,
       0.0404,  0.2000, -0.0344, -0.0128,  0.0068,
       0.0093,  0.2268, -0.0402,  0.0028,  0.0013,
      -0.0108,  0.2272, -0.0248,  0.0112, -0.0028,
      -0.0198,  0.1992,  0.0172,  0.0072, -0.0038,
      -0.0191,  0.1468,  0.0822, -0.0084, -0.0015), 10, byrow = TRUE),
    mid = matrix(c(
      -0.0117,  0.0804,  0.1570, -0.0284,  0.0027,
      -0.0020,  0.0160,  0.2200, -0.0400,  0.0060,
       0.0050, -0.0280,  0.2460, -0.0280,  0.0050,
       0.0060, -0.0400,  0.2200,  0.0160, -0.0020,
       0.0027, -0.0284,  0.1570,  0.0804, -0.0117), 5, byrow = TRUE)),

  beers_modified = list(
    head = matrix(c(
       0.3332, -0.1938,  0.0702, -0.0118,  0.0022,
       0.2569, -0.0753,  0.0205, -0.0027,  0.0006,
       0.1903,  0.0216, -0.0146,  0.0032, -0.0005,
       0.1334,  0.0969, -0.0351,  0.0059, -0.0011,
       0.0862,  0.1506, -0.0410,  0.0054, -0.0012,
       0.0486,  0.1831, -0.0329,  0.0021, -0.0009,
       0.0203,  0.1955, -0.0123, -0.0031, -0.0004,
       0.0008,  0.1893,  0.0193, -0.0097,  0.0003,
      -0.0108,  0.1677,  0.0577, -0.0153,  0.0007,
      -0.0159,  0.1354,  0.0972, -0.0170,  0.0003), 10, byrow = TRUE),
    mid = matrix(c(
      -0.0160,  0.0973,  0.1321, -0.0121, -0.0013,
      -0.0129,  0.0590,  0.1564,  0.0018, -0.0043,
      -0.0085,  0.0260,  0.1650,  0.0260, -0.0085,
      -0.0043,  0.0018,  0.1564,  0.0590, -0.0129,
      -0.0013, -0.0121,  0.1321,  0.0973, -0.0160), 5, byrow = TRUE))
)

# The (5k x k) matrix taking k closed five-year groups to 5k single years.
.grad_matrix <- function(method, k) {
  if (method == "karup_king") {
    head <- as.matrix(DemogAnr::first_coef)
    mid  <- as.matrix(DemogAnr::middle_coef)
    tail <- as.matrix(DemogAnr::last_coef)
  } else {
    head <- .grad_panels[[method]]$head
    mid  <- .grad_panels[[method]]$mid
    tail <- head[nrow(head):1, ncol(head):1]
  }
  w  <- ncol(mid)            # groups each panel draws on
  nh <- nrow(head) / 5       # groups covered by the head (and tail) panel
  half <- (w - 1) / 2
  M <- matrix(0, 5 * k, k)
  M[seq_len(5 * nh), seq_len(w)] <- head
  for (g in (nh + 1):(k - nh)) {
    M[(5 * g - 4):(5 * g), (g - half):(g + half)] <- mid
  }
  M[(5 * (k - nh) + 1):(5 * k), (k - w + 1):k] <- tail
  M
}

.grad_name <- function(method) {
  switch(method,
         sprague = "Sprague",
         beers_ordinary = "Beers ordinary",
         beers_modified = "Beers modified",
         karup_king = "Karup-King")
}


#' Graduate Five-Year Age Groups into Single Years of Age
#'
#' Splits a population (or deaths, or any count) in five-year age groups into
#' single years of age with the classical osculatory interpolation
#' multipliers of Sprague, Beers or Karup-King.
#'
#' @details
#' Each single-year value is a fixed linear combination of the five-year group
#' it belongs to and its neighbours. The weights come from a polynomial that
#' passes through the cumulated distribution and joins smoothly onto its
#' neighbours at the group boundaries, which is why the result has no breaks
#' where one group ends and the next begins. Separate sets of weights are used
#' for the two youngest and the two oldest groups, where there are not two
#' neighbours on each side.
#'
#' \describe{
#'   \item{`"sprague"`}{Sprague's (1880) fifth-difference formula, the
#'     default. The single years of each group add back exactly to the group,
#'     and a distribution that is a cubic in single years of age is reproduced
#'     exactly.}
#'   \item{`"beers_ordinary"`}{Beers' ordinary six-term formula, which
#'     minimises the fifth differences of the interpolated series. It has the
#'     same two properties as Sprague's and gives a markedly smoother curve: on
#'     census distributions its squared fifth differences are typically a
#'     tenth to a quarter of Sprague's.}
#'   \item{`"beers_modified"`}{Beers' (1945) modified formula, which combines
#'     interpolation with some smoothing. The single years of a group no
#'     longer add up to the reported group, so it suits data with age heaping
#'     that should not be passed straight into single years. The overall total
#'     is still kept exactly: each group's population is redistributed in
#'     full, partly into the neighbouring groups. It reproduces a quadratic
#'     exactly, but not a cubic.}
#'   \item{`"karup_king"`}{The Karup-King third-difference formula, with
#'     panels three groups wide. The group totals are kept, but it reproduces
#'     only a linear distribution exactly. These are the coefficients that
#'     [karup_king()] uses; `graduate()` returns them without rounding.}
#' }
#'
#' Graduation assumes that the five-year groups are themselves correct. If
#' they show the saw-tooth pattern of age misreporting, smooth them first with
#' [smooth_ages()]; graduating heaped groups faithfully reproduces the heaping
#' in single years.
#'
#' Where the population falls steeply from one group to the next, as at the
#' oldest ages, the negative coefficients can produce a negative count for a
#' single year. `graduate()` warns when this happens rather than altering the
#' result, which would break the group totals.
#'
#' @param data A data frame with one row per five-year age group, in increasing
#'   order of age.
#' @param age_col Column name for the age group: either its lower bound (0, 5,
#'   10, ...) or a label such as `"0-4"` or `"80+"`.
#' @param pop_col Column name for the count to be graduated.
#' @param method One of `"sprague"` (default), `"beers_ordinary"`,
#'   `"beers_modified"` or `"karup_king"`; see Details.
#' @param open_ended Logical; if `TRUE` (default) the last age group is
#'   open-ended. It is carried over as a single row and is not used in the
#'   interpolation.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   single-year values against the reported group averages is attached as
#'   `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_graduation`: a list with a `table` of single
#'   years of `age` and the graduated `population` (with `open_ended = TRUE`
#'   the last row is the open-ended group, carried over unchanged); a `groups`
#'   table comparing each `reported` group with the sum of its single years
#'   (`graduated`) and the `difference`; and the `method`.
#'
#' @examples
#' # The bundled five-year distribution
#' data(grouped_pop)
#' g <- graduate(grouped_pop, age_col = "age_col", pop_col = "2021",
#'               open_ended = FALSE)
#' g
#' head(g$table, 10)
#'
#' # Ghana, 1960 census: smooth the heaped groups, then graduate
#' ghana1960 <- data.frame(
#'   age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
#'           "40-44","45-49","50-54","55-59","60-64","65+"),
#'   females = c(654258, 503070, 323460, 265534, 322576, 306329, 245883,
#'               179182, 145572, 95590, 81715, 48412, 54572, 100392)
#' )
#' s <- smooth_ages(ghana1960, "age", "females", graph = FALSE)
#' graduate(s$table, "age", "smoothed", graph = FALSE)
#'
#' @references
#' Sprague, T. B. (1880). Explanation of a new formula for interpolation. \emph{Journal of the Institute of Actuaries}, 22(4), 270-285.
#'
#' Beers, H. S. (1945). Modified-interpolation formulas that minimize fourth differences. \emph{Record of the American Institute of Actuaries}, 34(69).
#'
#' Judson, D. H., & Popoff, C. L. (2004). Selected general methods. In J. S. Siegel & D. A. Swanson (Eds.), \emph{The Methods and Materials of Demography} (2nd ed., Appendix C). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (The multiplier tables.)
#'
#' Shryock, H. S., Siegel, J. S., & Larmon, E. A. (1973). \emph{The Methods and Materials of Demography}. Washington, DC: US Bureau of the Census. (Chapter 22 and Appendix C: interpolation and graduation.)
#'
#' @seealso [smooth_ages()] to correct the groups first, [karup_king()] for
#'   the rounded Karup-King split of several columns at once, and [whipple()]
#'   or [myers()] to check the single-year result.
#' @export
graduate <- function(data, age_col, pop_col,
                     method = c("sprague", "beers_ordinary", "beers_modified",
                                "karup_king"),
                     open_ended = TRUE, graph = TRUE) {
  method <- match.arg(method)
  d <- .age_dist(data, age_col, pop_col, open_ended, width = 5)
  P <- d$pop[d$closed]
  age <- d$lower[d$closed]
  k <- length(P)
  need <- if (method == "karup_king") 3 else 5
  if (k < need) {
    stop("The ", .grad_name(method), " method needs at least ", need,
         " closed five-year age groups; ", k, " given.")
  }

  single <- as.vector(.grad_matrix(method, k) %*% P)
  single_age <- age[1] + seq_len(5 * k) - 1
  if (any(single < 0)) {
    warning("Graduation produced negative values at age(s) ",
            paste(single_age[single < 0], collapse = ", "),
            ". This happens where a group is very small next to its ",
            "neighbours; consider a larger open-ended group.", call. = FALSE)
  }

  grp_sum <- colSums(matrix(single, nrow = 5))
  groups <- data.frame(age = age, interval = d$label[d$closed],
                       reported = P, graduated = grp_sum,
                       difference = grp_sum - P, stringsAsFactors = FALSE)

  tab <- data.frame(age = c(single_age, if (open_ended) d$lower[d$k]),
                    population = c(single, if (open_ended) d$pop[d$k]))

  out <- list(table = tab, groups = groups, method = method,
              open_ended = open_ended,
              open_label = if (open_ended) d$label[d$k] else NULL)
  class(out) <- "dem_graduation"
  if (graph) {
    out$plot <- .plot_graduation(tab[seq_len(5 * k), ], groups,
                                 .grad_name(method))
  }
  out
}

#' @export
print.dem_graduation <- function(x, ...) {
  cat(sprintf("Graduation into single years of age: %s multipliers\n",
              .grad_name(x$method)))
  cat(strrep("-", 58), "\n", sep = "")
  gr <- x$groups
  fmt_n <- function(v) format(round(v), big.mark = ",", trim = TRUE,
                              scientific = FALSE)
  disp <- data.frame(age = gr$interval, reported = fmt_n(gr$reported),
                     `sum of single years` = fmt_n(gr$graduated),
                     difference = sprintf("%+.1f", gr$difference),
                     check.names = FALSE, stringsAsFactors = FALSE)
  print(disp, row.names = FALSE)
  cat(strrep("-", 58), "\n", sep = "")
  n_single <- nrow(gr) * 5
  cat(sprintf("  %d single years, ages %g-%g", n_single, gr$age[1],
              gr$age[1] + n_single - 1))
  if (isTRUE(x$open_ended)) cat(sprintf(", plus the open group %s", x$open_label))
  cat("; see $table\n")
  if (x$method == "beers_modified") {
    cat("  Beers' modified formula smooths as it interpolates, so the groups\n")
    cat("  are not reproduced exactly; the overall total is.\n")
  }
  invisible(x)
}

#' @export
plot.dem_graduation <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call graduate(..., graph = TRUE).")
  x$plot
}
