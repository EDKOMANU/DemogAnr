#' Compute Life Table Using Age and nqx
#'
#' This function computes a life table from a data frame containing age groups and
#' age-interval probabilities of dying (\code{nqx}). It calculates key life table
#' components such as \code{lx}, \code{dx}, \code{Lx}, \code{Tx}, and \code{ex}.
#'
#' Interval widths \eqn{n} are derived from the age column (e.g. ages 0, 1, 5,
#' 10, ... give widths 1, 4, 5, ...). Person-years lived in each closed
#' interval are approximated as \eqn{{}_nL_x = n(l_x - d_x/2)}, i.e. deaths are
#' assumed to occur at the midpoint of the interval. For the open-ended last
#' age group, person-years are approximated with the classical rule of thumb
#' \eqn{L_{x+} = l_x \log_{10}(l_x)} (equivalent to
#' \eqn{e_x \approx \log_{10} l_x}), which is used when no death rate for the
#' open interval is available. If death rates are available, [lifetable()]
#' provides a more precise treatment of the open interval.
#'
#' @param data A data frame containing age and nqx columns.
#' @param age A string specifying the column name for age groups (lower bounds,
#'   e.g. 0, 1, 5, 10, ...). Ages must be strictly increasing.
#' @param nqx A string specifying the column name for the probabilities of dying
#'   between age \code{x} and \code{x+n} (\code{nqx}). The last value should be
#'   1 (open-ended group).
#' @param l0 Radix of the life table; survivors at age 0 (default = 100,000).
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during computation.
#'
#' @return A data frame with the computed life table, including:
#'   \describe{
#'     \item{\code{Age}}{Age groups.}
#'     \item{\code{n}}{Width of each age interval (\code{Inf} for the last, open-ended group).}
#'     \item{\code{nqx}}{Probability of dying between age \code{x} and \code{x+n}.}
#'     \item{\code{lx}}{Number of survivors at the beginning of each age group.}
#'     \item{\code{dx}}{Number of deaths within each age group.}
#'     \item{\code{Lx}}{Person-years lived within each age group.}
#'     \item{\code{Tx}}{Cumulative person-years lived above a specific age group.}
#'     \item{\code{ex}}{Life expectancy at the beginning of each age group.}
#'   }
#'
#' @examples
#' # Example data
#' data <- data.frame(
#'   age = c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85),
#'   nqx = c(0.05, 0.01, 0.005, 0.002, 0.003, 0.004, 0.005, 0.007, 0.010, 0.015,
#'           0.020, 0.030, 0.050, 0.070, 0.100, 0.150, 0.200, 0.300, 1.000)
#' )
#' life_table <- lifetable_nqx(data, age = "age", nqx = "nqx")
#' print(life_table)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 3: The Life Table and Single Decrement Processes.)
#'
#' Chiang, C. L. (1984). \emph{The Life Table and Its Applications}. Malabar, FL: Robert E. Krieger Publishing.
#'
#' @export
lifetable_nqx <- function(data, age, nqx, l0 = 100000, verbose = FALSE) {
  # Validate input
  if (!is.data.frame(data)) stop("Input must be a data frame.")
  if (!age %in% colnames(data) || !nqx %in% colnames(data)) {
    stop("Specified columns for age and nqx must exist in the data frame.")
  }

  if (verbose) {
    message("Computing life table from nqx probabilities...")
  }

  # Extract age and nqx
  age <- data[[age]]
  nqx <- data[[nqx]]

  if (!is.numeric(age) || !is.numeric(nqx)) stop("Specified age and nqx columns must be numeric.")
  if (any(is.na(age)) || any(is.na(nqx))) stop("age and nqx must not contain missing values.")
  if (is.unsorted(age, strictly = TRUE)) stop("Ages must be strictly increasing.")
  if (any(nqx < 0 | nqx > 1)) stop("nqx values must be between 0 and 1.")

  n <- length(age)

  # Interval widths; the last age group is open-ended
  width <- c(diff(age), Inf)

  # Initialize life table components
  lx <- numeric(n)
  dx <- numeric(n)
  Lx <- numeric(n)
  Tx <- numeric(n)

  # Compute lx and dx
  lx[1] <- l0
  for (i in 1:(n - 1)) {
    dx[i] <- lx[i] * nqx[i]
    lx[i + 1] <- lx[i] - dx[i]
  }
  dx[n] <- lx[n] # Everyone dies in the open-ended last age group

  # Compute Lx for closed intervals: n * (lx - dx/2), deaths at mid-interval
  for (i in 1:(n - 1)) {
    Lx[i] <- width[i] * (lx[i] - dx[i] / 2)
  }
  # Open-ended group: classical approximation L = lx * log10(lx)
  Lx[n] <- if (lx[n] > 1) lx[n] * log10(lx[n]) else 0

  # Compute Tx
  Tx[n] <- Lx[n]
  for (i in (n - 1):1) {
    Tx[i] <- Tx[i + 1] + Lx[i]
  }

  # Compute ex
  ex <- ifelse(lx > 0, Tx / lx, NA_real_)

  # Return a data frame
  life_table <- data.frame(
    Age = age,
    n = width,
    nqx = nqx,
    lx = lx,
    dx = dx,
    Lx = Lx,
    Tx = Tx,
    ex = ex
  )

  return(life_table)
}
