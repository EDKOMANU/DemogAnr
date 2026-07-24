#' @title lifetable
#' @name lifetable
#'
#' @description
#' This function computes a complete (or abridged) life table from given
#' mortality rates or from raw deaths and population counts. It returns the
#' standard life-table columns using textbook notation.
#'
#' @details
#' Interval widths \eqn{n} are derived from the age column (e.g. ages 0, 1, 5,
#' 10, ... give widths 1, 4, 5, ...). Probabilities of dying are obtained from
#' the death rates with the standard conversion
#' \eqn{{}_nq_x = n \cdot {}_nM_x / (1 + (n - {}_na_x) \cdot {}_nM_x)}, and
#' person-years lived are \eqn{{}_nL_x = n \cdot l_{x+n} + {}_na_x \cdot {}_nd_x}
#' for closed intervals and \eqn{l_x / M_x} for the open-ended interval.
#'
#' The average person-years lived by those dying in an interval,
#' \eqn{{}_na_x}, controls the accuracy of the \eqn{{}_nM_x \to {}_nq_x}
#' conversion. `lifetable()` offers three ways to set it, chosen with
#' `nax_method`:
#' \describe{
#'   \item{`"cd"` (default)}{Coale-Demeny model values for ages 0 and 1-4
#'     (sex-specific functions of \eqn{{}_1M_0}; Preston et al. 2001, Table
#'     3.3), \eqn{n/2} for other closed intervals, and \eqn{1/{}_nM_x} for the
#'     open interval. This reproduces the behaviour of earlier versions.}
#'   \item{`"keyfitz"`}{As `"cd"` for ages 0 and 1-4 and for the open interval,
#'     but for the interior five-year intervals \eqn{{}_na_x} is refined with
#'     the iterative Keyfitz (1966) formula
#'     \eqn{{}_na_x = n/2 + (n/24)({}_nd_{x+n} - {}_nd_{x-n})/{}_nd_x},
#'     which borrows information on the slope of the death distribution from
#'     the neighbouring age groups. The table is solved iteratively (a few
#'     iterations) because \eqn{{}_nd_x} itself depends on \eqn{{}_na_x}.}
#'   \item{user-supplied}{Pass a numeric vector to `nax` (one value per age
#'     group) to borrow \eqn{{}_na_x} values from any external model life table
#'     system, e.g. Coale-Demeny or the United Nations model life tables. This
#'     overrides `nax_method`.}
#' }
#'
#' @param data A data frame containing age-specific mortality data.
#' @param age String: name of the column with the lower bound of each age group, e.g., 0, 1, 5, 10, ...
#' @param nMx String (optional): name of the column with observed mortality rates. If not provided, `pop` and `Dx` must be given to compute it.
#' @param pop String (optional): name of the column with total population (exposure) in each age group. Used to compute `nMx` if not given.
#' @param Dx String (optional): name of the column with deaths in each age group. Used to compute `nMx` if not given.
#' @param sex Character: "male" (default) or "female". Determines the
#'   Coale-Demeny coefficients used for \eqn{a_0} and \eqn{{}_4a_1}.
#' @param nax_method Character: how to set \eqn{{}_na_x} for interior age
#'   groups. Either `"cd"` (Coale-Demeny young ages, \eqn{n/2} elsewhere;
#'   default) or `"keyfitz"` (iterative Keyfitz 1966 refinement for interior
#'   five-year intervals). Ignored if `nax` is supplied.
#' @param nax Optional numeric vector of \eqn{{}_na_x} values (one per age
#'   group) borrowed from an external model life table. Overrides `nax_method`.
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
#' # Refine nax for adult ages with the iterative Keyfitz method
#' res2 <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
#'                   nax_method = "keyfitz")
#' head(res2$lifetable[, c("Age", "nax", "nqx", "ex")])
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 3: The Life Table and Single Decrement Processes; Table 3.3 for the Coale-Demeny separation factors.)
#'
#' Keyfitz, N. (1966). A life table that agrees with the data. \emph{Journal of the American Statistical Association}, 61(314), 305-312. \doi{10.1080/01621459.1966.10480871}
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
                      nax_method = c("cd", "keyfitz"),
                      nax = NULL,
                      radix = 100000,
                      verbose = FALSE) {
  sex <- match.arg(sex)
  nax_method <- match.arg(nax_method)

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
    message("  nax method: ", if (!is.null(nax)) "user-supplied" else nax_method)
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

  if (!is.null(nax) && length(nax) != k) {
    stop("'nax' must have one value per age group (length ", k, ").")
  }

  # --- Baseline nax: Coale-Demeny young ages, n/2 elsewhere, 1/Mx open ---
  m0 <- nMx[1]
  if (sex == "male") {
    a0 <- if (m0 >= 0.107) 0.330 else 0.045 + 2.684 * m0
    a1 <- if (m0 >= 0.107) 1.352 else 1.651 - 2.816 * m0
  } else {
    a0 <- if (m0 >= 0.107) 0.350 else 0.053 + 2.800 * m0
    a1 <- if (m0 >= 0.107) 1.361 else 1.522 - 1.518 * m0
  }

  nax_vec <- n / 2
  if (age[1] == 0 && n[1] == 1) nax_vec[1] <- a0
  if (k >= 2 && age[2] == 1 && n[2] == 4) nax_vec[2] <- a1
  nax_vec[k] <- 1 / nMx[k] # open-ended interval

  # Helper: given an nax vector, return the life-table death counts ndx
  compute_dx <- function(nax_vec) {
    q <- (n * nMx) / (1 + (n - nax_vec) * nMx)
    q[k] <- 1
    q <- pmin(q, 1)
    lx <- numeric(k)
    lx[1] <- radix
    for (i in 2:k) lx[i] <- lx[i - 1] * (1 - q[i - 1])
    list(nqx = q, lx = lx, dx = lx * q)
  }

  if (!is.null(nax)) {
    nax_vec <- as.numeric(nax)
  } else if (nax_method == "keyfitz") {
    # Interior five-year intervals whose immediate neighbours are also width 5
    interior <- which(n == 5 &
                       c(FALSE, n[-k] == 5) &
                       c(n[-1] == 5, FALSE))
    interior <- interior[interior > 1 & interior < k]
    if (length(interior) > 0) {
      for (iter in 1:20) {
        dx <- compute_dx(nax_vec)$dx
        new_nax <- nax_vec
        for (i in interior) {
          if (dx[i] > 0) {
            val <- 2.5 + (5 / 24) * (dx[i + 1] - dx[i - 1]) / dx[i]
            new_nax[i] <- min(max(val, 0), n[i]) # keep within (0, n)
          }
        }
        if (max(abs(new_nax - nax_vec)) < 1e-7) {
          nax_vec <- new_nax
          break
        }
        nax_vec <- new_nax
      }
    }
  }

  lt <- compute_dx(nax_vec)

  lifetable <- data.frame(
    Age = age,
    n = n,
    nDx = Dx,
    nMx = nMx,
    nax = nax_vec,
    nqx = lt$nqx,
    lx = lt$lx,
    dx = lt$dx,
    Lx = NA_real_,
    Tx = NA_real_,
    ex = NA_real_
  )

  # Compute nLx = n*l(x+n) + nax*ndx for closed intervals
  for (i in 1:(k - 1)) {
    lifetable$Lx[i] <- n[i] * lifetable$lx[i + 1] + lifetable$nax[i] * lifetable$dx[i]
  }
  # Open-ended age group: Lx = lx / Mx
  lifetable$Lx[k] <- if (nMx[k] > 0) lifetable$lx[k] / nMx[k] else 0

  # Compute Tx and ex
  lifetable$Tx <- rev(cumsum(rev(lifetable$Lx)))
  lifetable$ex <- lifetable$Tx / lifetable$lx

  metrics <- list(
    TotalDeaths = sum(lifetable$dx, na.rm = TRUE),
    LifeExpectancyAtBirth = lifetable$ex[1]
  )

  return(list(metrics = metrics, lifetable = lifetable))
}
