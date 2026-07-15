#' @title lifetable
#' @name lifetable
#'
#' @description
#' This function computes a complete life table based on given mortality rates or raw mortality data.
#' It takes a dataframe with the necessary inputs and computes standard life table metrics.
#'
#' @details
#' Interval widths \eqn{n} are derived from the age column (e.g. ages 0, 1, 5,
#' 10, ... give widths 1, 4, 5, ...). Probabilities of dying are obtained from
#' the death rates with the standard conversion
#' \eqn{{}_nq_x = n \cdot {}_nM_x / (1 + (n - {}_na_x) \cdot {}_nM_x)}.
#' The average person-years lived by those dying in the interval
#' (\eqn{{}_na_x}) use the Coale-Demeny model values for ages 0 and 1-4
#' (which depend on \code{sex} and the level of infant mortality; see Preston
#' et al. 2001, Table 3.3), \eqn{n/2} for other closed intervals, and
#' \eqn{1/{}_nM_x} for the open-ended interval. Person-years lived are
#' \eqn{{}_nL_x = n \cdot l_{x+n} + {}_na_x \cdot {}_nd_x} for closed
#' intervals and \eqn{l_x / M_x} for the open-ended interval.
#'
#' @param data A data frame containing age-specific mortality data.
#' @param age String: name of the column with the lower bound of each age group, e.g., 0, 1, 5, 10, ...
#' @param nMx String (optional): name of the column with observed mortality rates. If not provided, `pop` and `Dx` must be given to compute it.
#' @param pop String (optional): name of the column with total population (exposure) in each age group. Used to compute `nMx` if not given.
#' @param Dx String (optional): name of the column with deaths in each age group. Used to compute `nMx` if not given.
#' @param sex Character: "male" (default) or "female". Determines the
#'   Coale-Demeny coefficients used for \eqn{a_0} and \eqn{{}_4a_1}.
#' @param radix Numeric: the life table radix, i.e. survivors at exact age 0 (default 100,000).
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during computation.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{`metrics`}{A list with the life expectancy at birth and the total life table deaths.}
#'   \item{`lifetable`}{A dataframe of the complete life table.}
#' }
#'
#' @examples
#' # Example computing rates from deaths and population
#' data(gphc2010)
#' res <- lifetable(data = gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
#' res$metrics
#' head(res$lifetable)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 3: The Life Table and Single Decrement Processes; Table 3.3 for the Coale-Demeny separation factors.)
#'
#' Chiang, C. L. (1984). \emph{The Life Table and Its Applications}. Malabar, FL: Robert E. Krieger Publishing.
#'
#' Coale, A. J., & Demeny, P. (1983). \emph{Regional Model Life Tables and Stable Populations} (2nd ed.). New York: Academic Press. (Source of the a0 and 4a1 separation factors.)
#'
#' @export
lifetable <- function(data,
                      age = "Age",
                      nMx = NULL,
                      pop = NULL,
                      Dx = NULL,
                      sex = c("male", "female"),
                      radix = 100000,
                      verbose = FALSE) {
  sex <- match.arg(sex)

  # Check if required columns are present
  if (is.null(nMx) && (is.null(pop) || is.null(Dx))) {
    stop("If nMx is not provided, both pop and Dx must be specified.")
  }

  if (verbose) {
    message("Computing complete life table...")
    if (is.null(nMx)) {
      message("  nMx not provided; computing from deaths and population.")
    } else {
      message("  Using observed nMx values.")
    }
  }

  # Extract columns
  age <- data[[age]]
  nMx <- if (!is.null(nMx)) {
    data[[nMx]]
  } else {
    # Compute nMx if not provided
    data[[Dx]] / data[[pop]]
  }

  if (any(is.na(age)) || any(is.na(nMx))) stop("age and nMx must not contain missing values.")
  if (is.unsorted(age, strictly = TRUE)) stop("Ages must be strictly increasing.")
  if (any(nMx < 0)) stop("Mortality rates must be non-negative.")

  Dx <- if (!is.null(Dx)) data[[Dx]] else NA

  k <- length(age)
  # Interval widths; the last age group is open-ended
  n <- c(diff(age), Inf)

  # Initialize data frame for life table metrics
  lifetable <- data.frame(
    Age = age,
    n = n,
    nDx = Dx,
    nMx = nMx,
    nax = NA_real_,
    nqx = NA_real_,
    lx = NA_real_,
    dx = NA_real_,
    Lx = NA_real_,
    Tx = NA_real_,
    ex = NA_real_
  )

  # Compute nax
  # Coale-Demeny model values for a0 and 4a1 (Preston et al. 2001, Table 3.3),
  # as functions of 1m0; n/2 for other closed intervals; 1/Mx for the open one.
  m0 <- nMx[1]
  if (sex == "male") {
    a0 <- if (m0 >= 0.107) 0.330 else 0.045 + 2.684 * m0
    a1 <- if (m0 >= 0.107) 1.352 else 1.651 - 2.816 * m0
  } else {
    a0 <- if (m0 >= 0.107) 0.350 else 0.053 + 2.800 * m0
    a1 <- if (m0 >= 0.107) 1.361 else 1.522 - 1.518 * m0
  }

  lifetable$nax <- n / 2
  if (age[1] == 0 && n[1] == 1) lifetable$nax[1] <- a0
  if (k >= 2 && age[2] == 1 && n[2] == 4) lifetable$nax[2] <- a1
  # Open-ended interval
  lifetable$nax[k] <- 1 / nMx[k]

  # Compute nqx: nqx = n*Mx / (1 + (n - nax)*Mx); the open interval has qx = 1
  lifetable$nqx <- (n * nMx) / (1 + (n - lifetable$nax) * nMx)
  lifetable$nqx[k] <- 1
  lifetable$nqx <- pmin(lifetable$nqx, 1) # Ensure nqx doesn't exceed 1

  # Compute lx and dx
  lifetable$lx[1] <- radix
  for (i in 2:k) {
    lifetable$lx[i] <- lifetable$lx[i - 1] * (1 - lifetable$nqx[i - 1])
  }
  lifetable$dx <- lifetable$lx * lifetable$nqx

  # Compute nLx = n*l(x+n) + nax*ndx for closed intervals
  for (i in 1:(k - 1)) {
    lifetable$Lx[i] <- n[i] * lifetable$lx[i + 1] + lifetable$nax[i] * lifetable$dx[i]
  }
  # Open-ended age group: Lx = lx / Mx
  lifetable$Lx[k] <- if (nMx[k] > 0) lifetable$lx[k] / nMx[k] else 0

  # Compute Tx and ex
  lifetable$Tx <- rev(cumsum(rev(lifetable$Lx)))
  lifetable$ex <- lifetable$Tx / lifetable$lx

  # Summary metrics
  metrics <- list(
    TotalDeaths = sum(lifetable$dx, na.rm = TRUE),
    LifeExpectancyAtBirth = lifetable$ex[1]
  )

  return(list(metrics = metrics, lifetable = lifetable))
}
