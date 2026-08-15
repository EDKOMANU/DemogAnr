# Multiple-decrement (cause-specific) life tables and cause-deleted
# (associated single-decrement) life tables.
# Source: Preston, Heuveline & Guillot (2001), Chapter 4.

#' Multiple-Decrement (Cause-Specific) Life Table
#'
#' Splits an ordinary life table by cause of decrement. Given all-cause deaths
#' (or the sum of the cause-specific deaths) and the deaths attributable to each
#' cause, it distributes the life-table deaths \eqn{{}_nd_x} across causes in
#' proportion to the observed deaths, giving the cause-specific decrements
#' \eqn{{}_nd_x^i = {}_nd_x \cdot ({}_nD_x^i / {}_nD_x)} and probabilities
#' \eqn{{}_nq_x^i = {}_nd_x^i / \ell_x}. The proportion of a birth cohort that
#' will eventually die of each cause is \eqn{\sum_x {}_nd_x^i / \ell_0}.
#'
#' @param data A data frame with one row per age group.
#' @param age Column name for the (lower bound of the) age group.
#' @param causes Character vector of column names holding the cause-specific
#'   death counts.
#' @param pop Column name for the mid-year population (exposure), used to build
#'   the all-cause life table. Alternatively supply an existing life table
#'   through `lx` and `nqx`.
#' @param lx,nqx Optional column names of an existing all-cause life table
#'   (survivors and probabilities of dying). When given, these are used directly
#'   instead of building a life table from `pop`, which is then not required.
#' @param deaths Optional column name for all-cause deaths. Defaults to the row
#'   sum of the `causes` columns.
#' @param sex,nax_method,radix Passed to [lifetable()] for the all-cause table.
#' @param graph Logical; if `TRUE` (default) a plot of the cohort distribution
#'   of deaths by cause is attached and shown on printing.
#'
#' @return An object of class `dem_mdlt`: a list with the all-cause life
#'   expectancy `e0`, a named vector `cause_distribution` (proportion of the
#'   birth cohort dying of each cause, summing to 1), and a `table` of the
#'   life-table columns with the cause-specific decrements `d_<cause>`.
#'
#' @examples
#' # Preston et al. (2001) Box 4.1: US females, 1991, deaths from neoplasms
#' # against all other causes, using the published all-cause life table.
#' us91 <- data.frame(
#'   Age  = c(0, 1, seq(5, 85, 5)),
#'   Dall = c(15758, 3169, 1634, 1573, 3955, 4948, 6491, 9428, 12027, 15543,
#'            19264, 25384, 37211, 59431, 88087, 114693, 143554, 164986, 320578),
#'   Dneo = c(63, 275, 268, 217, 318, 467, 856, 1924, 3532, 5958,
#'            8434, 11673, 17078, 25263, 33534, 36695, 36571, 30220, 32739),
#'   lx   = c(100000, 99217, 99050, 98959, 98870, 98637, 98379, 98070, 97653,
#'            97083, 96289, 95008, 93018, 89882, 85249, 78711, 69618, 57486, 41756),
#'   nqx  = c(0.00783, 0.00168, 0.00092, 0.00090, 0.00236, 0.00262, 0.00314,
#'            0.00425, 0.00584, 0.00818, 0.01330, 0.02095, 0.03371, 0.05155,
#'            0.07669, 0.11552, 0.17427, 0.27363, 1.00000))
#' us91$Dother <- us91$Dall - us91$Dneo
#' # 21.2% of newborns would eventually die of neoplasms (Preston, Box 4.1)
#' multiple_decrement(us91, age = "Age", causes = c("Dneo", "Dother"),
#'                    deaths = "Dall", lx = "lx", nqx = "nqx", graph = FALSE)
#'
#' # Or build the life table from population and deaths (Ghana, 2010)
#' data(gphc2010)
#' g <- gphc2010
#' g$circ  <- round(g$Deaths * pmin(0.55, 0.04 + 0.012 * (seq_len(nrow(g)) - 1)))
#' g$other <- g$Deaths - g$circ
#' multiple_decrement(g, age = "Age", causes = c("circ", "other"),
#'                    pop = "Pop", graph = FALSE)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' Chapter 4.
#'
#' @seealso [cause_deleted_lt()], [lifetable()]
#' @export
multiple_decrement <- function(data, age, causes, pop = NULL, lx = NULL,
                               nqx = NULL, deaths = NULL,
                               sex = c("male", "female"),
                               nax_method = c("cd", "keyfitz"),
                               radix = 1e5, graph = TRUE) {
  sex <- match.arg(sex); nax_method <- match.arg(nax_method)
  if (!all(causes %in% names(data))) stop("Some 'causes' columns not found.")
  D <- as.matrix(data[causes]); storage.mode(D) <- "double"
  Dtot <- if (!is.null(deaths)) as.numeric(data[[deaths]]) else rowSums(D)

  if (!is.null(lx) && !is.null(nqx)) {
    # Use a supplied all-cause life table directly.
    l <- as.numeric(data[[lx]]); q <- as.numeric(data[[nqx]])
    lt <- data.frame(Age = as.numeric(data[[age]]), lx = l, dx = l * q,
                     nqx = q, ex = NA_real_)
    radix <- l[1]
  } else if (!is.null(pop)) {
    base_in <- data.frame(Age = as.numeric(data[[age]]),
                          Pop = as.numeric(data[[pop]]), nDx = Dtot)
    lt <- lifetable(base_in, age = "Age", pop = "Pop", Dx = "nDx", sex = sex,
                    nax_method = nax_method, radix = radix, graph = FALSE)$lifetable
  } else {
    stop("Supply either 'pop', or an existing life table via 'lx' and 'nqx'.")
  }

  R <- D / ifelse(Dtot == 0, NA, Dtot)      # proportion of deaths by cause
  R[is.na(R)] <- 0
  dxi <- lt$dx * R                          # cause-specific decrements
  qxi <- dxi / lt$lx
  cause_distribution <- colSums(dxi) / radix
  names(cause_distribution) <- causes

  tab <- data.frame(Age = lt$Age, lx = round(lt$lx, 1), dx = round(lt$dx, 1))
  for (cn in causes) tab[[paste0("d_", cn)]] <- round(dxi[, cn], 1)
  for (cn in causes) tab[[paste0("q_", cn)]] <- qxi[, cn]

  out <- list(e0 = lt$ex[1], cause_distribution = cause_distribution,
              radix = radix, table = tab)
  class(out) <- "dem_mdlt"
  if (graph) {
    out$plot <- .plot_series(names(cause_distribution), cause_distribution,
                             ylab = "Proportion of cohort",
                             title = "Distribution of deaths by cause",
                             geom = "col")
  }
  out
}

#' Cause-Deleted (Associated Single-Decrement) Life Table
#'
#' Builds the life table that would result if a single cause of death were
#' eliminated, and reports the resulting gain in life expectancy at birth.
#' Under the standard proportionality assumption, the cause-deleted survival
#' probability in an interval is
#' \eqn{{}_np_x^{-i} = ({}_np_x)^{\,1 - R_x^i}}, where
#' \eqn{R_x^i = {}_nD_x^i / {}_nD_x} is the share of deaths due to the cause;
#' the remaining life-table columns are rebuilt from these probabilities using
#' the all-cause separation factors \eqn{{}_na_x}.
#'
#' Alternative constructions exist, and the assumption is under the user's
#' control through `method`:
#' \describe{
#'   \item{`"proportional"`}{(default) the proportional-hazards form above,
#'     \eqn{{}_np_x^{-i} = ({}_np_x)^{1 - R_x^i}}, retaining the all-cause
#'     separation factors.}
#'   \item{`"rates"`}{subtracts the cause-specific death rates directly,
#'     \eqn{{}_nm_x^{-i} = {}_nm_x (1 - R_x^i)}, and rebuilds \eqn{{}_nq_x}
#'     from them by the usual conversion.}
#' }
#' A user-supplied `nax` vector overrides the separation factors in either
#' case. On the US 1991 female data of Preston et al. (2001, Box 4.2) both
#' methods give a life expectancy of 82.43-82.44 years in the absence of
#' neoplasms, against their published 82.46 (they additionally graduate the
#' cause-deleted \eqn{{}_na_x^{-i}}); the gain in \eqn{e_0} is 3.5 years either
#' way.
#'
#' @param data A data frame with one row per age group.
#' @param age Column name for the age group.
#' @param cause Column name of the cause to delete (its death counts).
#' @param pop Column name for the mid-year population (exposure). Alternatively
#'   supply an existing life table through `lx` and `nqx`.
#' @param lx,nqx Optional column names of an existing all-cause life table.
#'   When given, these are used directly instead of building one from `pop`.
#' @param nax_all Optional column name of the all-cause separation factors
#'   accompanying `lx`/`nqx` (defaults to half the interval width).
#' @param deaths Optional all-cause deaths column; defaults to `causes` sum if
#'   supplied, otherwise must be given.
#' @param causes Optional character vector of all cause columns, used to form
#'   all-cause deaths when `deaths` is not given.
#' @param method The cause-deletion assumption: `"proportional"` (default) or
#'   `"rates"`; see Details.
#' @param nax Optional vector of cause-deleted separation factors
#'   \eqn{{}_na_x^{-i}}, one per age group, overriding the all-cause values.
#' @param sex,nax_method,radix Passed to [lifetable()] for the all-cause table.
#' @param graph Logical; if `TRUE` (default) a plot of the years of life gained
#'   by age is attached and shown on printing.
#'
#' @return An object of class `dem_cause_deleted`: a list with `e0` (all-cause),
#'   `e0_deleted`, `gain` (`e0_deleted - e0`), the deleted `cause`, and a
#'   `table` comparing `nqx`/`ex` before and after deletion.
#'
#' @examples
#' # Preston et al. (2001) Box 4.2: US females 1991, deleting neoplasms.
#' # Published: e0 = 78.92 rises to 82.46 (a gain of 3.54 years).
#' us91 <- data.frame(
#'   Age  = c(0, 1, seq(5, 85, 5)),
#'   Dall = c(15758, 3169, 1634, 1573, 3955, 4948, 6491, 9428, 12027, 15543,
#'            19264, 25384, 37211, 59431, 88087, 114693, 143554, 164986, 320578),
#'   Dneo = c(63, 275, 268, 217, 318, 467, 856, 1924, 3532, 5958,
#'            8434, 11673, 17078, 25263, 33534, 36695, 36571, 30220, 32739),
#'   lx   = c(100000, 99217, 99050, 98959, 98870, 98637, 98379, 98070, 97653,
#'            97083, 96289, 95008, 93018, 89882, 85249, 78711, 69618, 57486, 41756),
#'   nqx  = c(0.00783, 0.00168, 0.00092, 0.00090, 0.00236, 0.00262, 0.00314,
#'            0.00425, 0.00584, 0.00818, 0.01330, 0.02095, 0.03371, 0.05155,
#'            0.07669, 0.11552, 0.17427, 0.27363, 1.00000),
#'   nax  = c(0.152, 1.605, 2.275, 2.843, 2.657, 2.547, 2.550, 2.616, 2.677,
#'            2.685, 2.681, 2.655, 2.647, 2.646, 2.631, 2.628, 2.618, 2.570, 6.539))
#' cause_deleted_lt(us91, age = "Age", cause = "Dneo", deaths = "Dall",
#'                  lx = "lx", nqx = "nqx", nax_all = "nax", graph = FALSE)
#'
#' # The competing assumption, on the same data
#' cause_deleted_lt(us91, age = "Age", cause = "Dneo", deaths = "Dall",
#'                  lx = "lx", nqx = "nqx", nax_all = "nax",
#'                  method = "rates", graph = FALSE)
#'
#' # Or from population and deaths (Ghana, 2010)
#' data(gphc2010)
#' g <- gphc2010
#' g$circ  <- round(g$Deaths * pmin(0.55, 0.04 + 0.012 * (seq_len(nrow(g)) - 1)))
#' g$other <- g$Deaths - g$circ
#' cause_deleted_lt(g, age = "Age", cause = "circ",
#'                  causes = c("circ", "other"), pop = "Pop", graph = FALSE)
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography:
#' Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers.
#' Chapter 4.
#'
#' @seealso [multiple_decrement()], [lifetable()]
#' @export
cause_deleted_lt <- function(data, age, cause, pop = NULL, lx = NULL,
                             nqx = NULL, nax_all = NULL, deaths = NULL,
                             causes = NULL,
                             method = c("proportional", "rates"), nax = NULL,
                             sex = c("male", "female"),
                             nax_method = c("cd", "keyfitz"),
                             radix = 1e5, graph = TRUE) {
  sex <- match.arg(sex); nax_method <- match.arg(nax_method)
  method <- match.arg(method)
  Dcause <- as.numeric(data[[cause]])
  if (!is.null(deaths)) {
    Dtot <- as.numeric(data[[deaths]])
  } else if (!is.null(causes)) {
    Dtot <- rowSums(as.matrix(data[causes]))
  } else {
    stop("Provide 'deaths' (all-cause) or 'causes' (to sum).")
  }

  if (!is.null(lx) && !is.null(nqx)) {
    # Use a supplied all-cause life table directly.
    l <- as.numeric(data[[lx]]); q <- as.numeric(data[[nqx]])
    ages <- as.numeric(data[[age]])
    w <- diff(ages); w <- c(w, w[length(w)])
    a_all <- if (!is.null(nax_all)) as.numeric(data[[nax_all]]) else w / 2
    Lx <- numeric(length(l))
    for (i in seq_len(length(l) - 1)) Lx[i] <- w[i] * l[i + 1] + a_all[i] * (l[i] * q[i])
    Lx[length(l)] <- l[length(l)] * a_all[length(l)]
    Tx <- rev(cumsum(rev(Lx)))
    base <- data.frame(Age = ages, n = w, nMx = (l * q) / ifelse(Lx == 0, NA, Lx),
                       nax = a_all, nqx = q, lx = l, dx = l * q,
                       Lx = Lx, Tx = Tx, ex = Tx / l)
    radix <- l[1]
  } else if (!is.null(pop)) {
    base_in <- data.frame(Age = as.numeric(data[[age]]),
                          Pop = as.numeric(data[[pop]]), nDx = Dtot)
    base <- lifetable(base_in, age = "Age", pop = "Pop", Dx = "nDx", sex = sex,
                      nax_method = nax_method, radix = radix, graph = FALSE)$lifetable
  } else {
    stop("Supply either 'pop', or an existing life table via 'lx' and 'nqx'.")
  }

  R <- Dcause / ifelse(Dtot == 0, NA, Dtot); R[is.na(R)] <- 0
  k <- nrow(base); n <- base$n
  nax_del <- if (!is.null(nax)) as.numeric(nax) else base$nax

  if (method == "proportional") {
    npx_del <- (1 - base$nqx)^(1 - R)         # proportional hazards
    qx_del  <- 1 - npx_del
  } else {                                     # "rates": subtract the rates
    m_del  <- base$nMx * (1 - R)
    qx_del <- (n * m_del) / (1 + (n - nax_del) * m_del)
    qx_del[k] <- 1
  }
  nax <- nax_del
  lx_del <- numeric(k); lx_del[1] <- radix
  for (i in 2:k) lx_del[i] <- lx_del[i - 1] * (1 - qx_del[i - 1])
  dx_del <- lx_del * qx_del
  Lx_del <- numeric(k)
  if (k > 1) for (i in 1:(k - 1)) Lx_del[i] <- n[i] * lx_del[i + 1] + nax[i] * dx_del[i]
  mMx_open <- base$nMx[k] * (1 - R[k])       # cause-deleted rate, open interval
  Lx_del[k] <- if (mMx_open > 0) lx_del[k] / mMx_open else lx_del[k] / base$nMx[k]
  Tx_del <- rev(cumsum(rev(Lx_del)))
  ex_del <- Tx_del / lx_del

  tab <- data.frame(Age = base$Age, nqx = base$nqx, nqx_deleted = qx_del,
                    ex = base$ex, ex_deleted = ex_del)
  out <- list(cause = cause, e0 = base$ex[1], e0_deleted = ex_del[1],
              gain = ex_del[1] - base$ex[1], radix = radix, table = tab)
  class(out) <- "dem_cause_deleted"
  if (graph) {
    out$plot <- .plot_series(base$Age, ex_del - base$ex,
                             ylab = "Years of life gained",
                             title = paste0("Life-expectancy gain by age from deleting '", cause, "'"),
                             geom = "col")
  }
  out
}

#' @export
print.dem_mdlt <- function(x, ...) {
  cat("Multiple-decrement (cause-specific) life table\n")
  cat(strrep("-", 52), "\n", sep = "")
  print(x$table, row.names = FALSE)
  cat(strrep("-", 52), "\n", sep = "")
  cat("Distribution of deaths by cause (per birth):\n")
  print(round(x$cause_distribution, 4))
  cat(sprintf("All-cause life expectancy at birth = %.2f\n", x$e0))
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}

#' @export
print.dem_cause_deleted <- function(x, ...) {
  cat(sprintf("Cause-deleted life table (deleting '%s')\n", x$cause))
  cat(strrep("-", 52), "\n", sep = "")
  print(round(x$table, 5), row.names = FALSE)
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("e0 (all causes) = %.2f   e0 (cause deleted) = %.2f   gain = %.2f years\n",
              x$e0, x$e0_deleted, x$gain))
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}

#' @export
plot.dem_mdlt <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call multiple_decrement(..., graph = TRUE).")
  x$plot
}

#' @export
plot.dem_cause_deleted <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call cause_deleted_lt(..., graph = TRUE).")
  x$plot
}
