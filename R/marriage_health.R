# Singulate mean age at marriage (Hajnal) and Sullivan health expectancy.

#' Singulate Mean Age at Marriage (SMAM)
#'
#' Estimates the mean age at first marriage from the proportions of people who
#' are single (never married) by age group, using Hajnal's (1953) method. No
#' marriage-registration data are required: the measure is derived entirely from
#' a census cross-tabulation of marital status by age.
#'
#' The proportion single at exact age 50 is taken as the average of the
#' proportions in the 45-49 and 50-54 groups and treated as the fraction who
#' never marry. The singulate mean age at marriage is
#' \deqn{SMAM = \frac{15 + 5 \sum_{x=15,\dots,45} s_x - 50\, s_{50}}{1 - s_{50}},}
#' the mean number of years lived single by those who marry before age 50.
#'
#' @param data A data frame with one row per five-year age group, spanning at
#'   least ages 15-19 through 50-54.
#' @param age Column name for the (lower bound of the) age group.
#' @param prop_single Column name for the proportion single. Supply this, or
#'   both `single` and `total`.
#' @param single,total Column names for the number single and the total number
#'   of people in each age group, used to form the proportion.
#'
#' @return An object of class `dem_smam`: a list with `SMAM`, the proportion
#'   never marrying `prop_never` (\eqn{s_{50}}), and the `table` of proportions
#'   single used.
#'
#' @examples
#' # Preston et al. (2001) Box 4.4: Turkish males, 1990 (SMAM = 25.0 years)
#' turkey <- data.frame(
#'   age    = seq(15, 50, 5),
#'   men    = c(3165061, 2581153, 2435765, 2096899, 1784121, 1418784,
#'              1111113, 980115),
#'   single = c(3030203, 1853222, 629077, 180767, 77134, 43412, 28627, 22527)
#' )
#' smam(turkey, age = "age", single = "single", total = "men")
#'
#' @references
#' Hajnal, J. (1953). Age at marriage and proportions marrying.
#' \emph{Population Studies}, 7(2), 111-136.
#'
#' @export
smam <- function(data, age, prop_single = NULL, single = NULL, total = NULL) {
  ages <- as.numeric(data[[age]])
  if (!is.null(prop_single)) {
    s <- as.numeric(data[[prop_single]])
  } else if (!is.null(single) && !is.null(total)) {
    s <- as.numeric(data[[single]]) / as.numeric(data[[total]])
  } else {
    stop("Supply 'prop_single', or both 'single' and 'total'.")
  }
  need <- c(seq(15, 45, 5), 50)
  if (!all(need %in% ages))
    stop("Ages 15-19 through 50-54 (lower bounds 15..50) must all be present.")

  s45 <- s[ages == 45]
  s50g <- s[ages == 50]
  s50 <- (s45 + s50g) / 2                       # proportion single at exact 50
  repro <- ages %in% seq(15, 45, 5)
  sum_s <- sum(s[repro])
  SMAM <- (15 + 5 * sum_s - 50 * s50) / (1 - s50)

  out <- list(SMAM = SMAM, prop_never = s50,
              table = data.frame(age = ages, prop_single = s))
  class(out) <- "dem_smam"
  out
}

#' Sullivan Health Expectancy
#'
#' Partitions life expectancy into years spent in good and in poor health using
#' the Sullivan (1971) method: age-specific prevalence of ill-health or
#' disability is applied to the person-years of a period life table. The health
#' expectancy at age \eqn{x} is
#' \deqn{HLE(x) = \frac{1}{\ell_x} \sum_{a \ge x} (1 - \pi_a)\, {}_nL_a,}
#' where \eqn{\pi_a} is the prevalence of ill-health in age group \eqn{a}. By
#' construction \eqn{HLE(x) + DLE(x) = e_x}.
#'
#' @param data A data frame with one row per age group, containing life-table
#'   `nLx` and `lx` columns and an age-specific prevalence column.
#' @param age Column name for the age group.
#' @param nLx,lx Column names for the life-table person-years and survivors.
#' @param prevalence Column name for the prevalence of ill-health/disability
#'   (a proportion between 0 and 1).
#' @param graph Logical; if `TRUE` (default) a plot of healthy vs. total life
#'   expectancy by age is attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_sullivan`: a list with `HLE0` (health
#'   expectancy at birth), `e0`, and a `table` of `ex`, `HLE`, `DLE`, and the
#'   proportion of life spent healthy by age.
#'
#' @examples
#' lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
#'                 graph = FALSE)$lifetable
#' lt$disab <- pmin(0.6, 0.02 + 0.006 * (seq_len(nrow(lt)) - 1))
#' sullivan_hle(lt, age = "Age", nLx = "Lx", lx = "lx",
#'              prevalence = "disab", graph = FALSE)
#'
#' @references
#' Sullivan, D. F. (1971). A single index of mortality and morbidity.
#' \emph{HSMHA Health Reports}, 86(4), 347-354.
#'
#' @seealso [lifetable()]
#' @export
sullivan_hle <- function(data, age, nLx, lx, prevalence, graph = TRUE) {
  ages <- as.numeric(data[[age]])
  L  <- as.numeric(data[[nLx]])
  l  <- as.numeric(data[[lx]])
  pi <- as.numeric(data[[prevalence]])
  if (any(pi < 0 | pi > 1, na.rm = TRUE)) stop("'prevalence' must be in [0, 1].")

  healthyL <- (1 - pi) * L
  Tx        <- rev(cumsum(rev(L)))
  Tx_health <- rev(cumsum(rev(healthyL)))
  ex  <- Tx / l
  HLE <- Tx_health / l
  DLE <- ex - HLE

  tab <- data.frame(age = ages, ex = ex, HLE = HLE, DLE = DLE,
                    prop_healthy = HLE / ex)
  out <- list(HLE0 = HLE[1], e0 = ex[1], table = tab)
  class(out) <- "dem_sullivan"
  if (graph) {
    out$plot <- .plot_series(ages, HLE, ylab = "Healthy life expectancy",
                             title = "Sullivan health expectancy by age",
                             geom = "col")
  }
  out
}

#' @export
print.dem_smam <- function(x, ...) {
  cat("Singulate mean age at marriage (Hajnal)\n")
  cat(strrep("-", 40), "\n", sep = "")
  print(round(x$table, 4), row.names = FALSE)
  cat(strrep("-", 40), "\n", sep = "")
  cat(sprintf("SMAM = %.2f years   proportion never marrying = %.3f\n",
              x$SMAM, x$prop_never))
  invisible(x)
}

#' @export
print.dem_sullivan <- function(x, ...) {
  cat("Sullivan health expectancy\n")
  cat(strrep("-", 46), "\n", sep = "")
  print(round(x$table, 3), row.names = FALSE)
  cat(strrep("-", 46), "\n", sep = "")
  cat(sprintf("At birth: e0 = %.2f, healthy = %.2f (%.1f%% of life)\n",
              x$e0, x$HLE0, 100 * x$HLE0 / x$e0))
  invisible(x)
}

#' @export
plot.dem_sullivan <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call sullivan_hle(..., graph = TRUE).")
  x$plot
}
