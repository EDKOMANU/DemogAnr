# Model life tables: the Coale-Demeny and United Nations systems, indexed by
# family, sex and level, and matched to an observed child mortality.

#' Model Life Table
#'
#' Returns a life table from one of the two standard model life table systems:
#' the four Coale-Demeny regional families, and the five United Nations
#' patterns for developing countries. A model table supplies a complete age
#' pattern of mortality where the data give only a fragment of one, which is
#' the usual position when mortality is estimated indirectly.
#'
#' @details
#' A model table is picked in one of two ways. Give a level directly with `e0`,
#' or give an observed child mortality with `q1` or `q5` and let the level that
#' matches it be found. The second is what pairs with [chm_brass()], which
#' estimates \eqn{q(1)}, \eqn{q(2)}, \eqn{q(3)}, \eqn{q(5)} and so on for a
#' named Coale-Demeny family but stops short of a life table: passing its
#' \eqn{q(5)} here completes the estimate.
#'
#' The packaged rates are tabulated at half-year steps of \eqn{e_0} from 54 to
#' 76. A value between two steps is linearly interpolated on
#' \eqn{\log {}_nM_x}, which keeps the rates positive and the age pattern
#' smooth; a value outside the range is an error rather than an extrapolation.
#' The life table itself is then built by [lifetable()], so the treatment of
#' \eqn{{}_na_x} and of the open interval is exactly the same as for a table
#' computed from observed deaths.
#'
#' The realised \eqn{e_0} of the returned table runs a little above the `e0`
#' requested, because the tabulated rates stop at an open interval of 80 and
#' up, which [lifetable()] closes as \eqn{l_x/{}_nM_x}. The gap grows with the
#' level, since more of the cohort survives into that interval: about +0.1
#' years at \eqn{e_0 = 54}, +0.2 at 60, +0.5 at 70 and +1.2 at 76. Both
#' figures are returned, as `e0_requested` and `e0_realised`, and it is
#' `e0_requested` that identifies the published table. Where the level itself
#' must be exact, match on `q1` or `q5` instead, which is reproduced exactly.
#'
#' @param family Model life table family: one of `"CD_West"` (default),
#'   `"CD_North"`, `"CD_South"`, `"CD_East"`, `"UN_General"`,
#'   `"UN_Latin_American"`, `"UN_Chilean"`, `"UN_South_Asian"` or
#'   `"UN_Far_Eastern"`. The Coale-Demeny names may also be given in the short
#'   form used by [chm_brass()]: `"west"`, `"north"`, `"south"`, `"east"`.
#' @param sex `"male"` (default) or `"female"`.
#' @param e0 Life expectancy at birth identifying the level, between 54 and 76.
#' @param q1,q5 Alternatively, an observed probability of dying by exact age 1
#'   or by exact age 5, in which case the level whose model table matches it is
#'   found. Supply exactly one of `e0`, `q1` and `q5`.
#' @param radix Numeric: survivors at exact age 0 (default 100,000).
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} survival curve is
#'   attached as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_model_lt`, which inherits from
#'   `dem_lifetable`: a list with `metrics`, the `lifetable` data frame, and
#'   `model`, recording the `family`, `sex`, `e0_requested`, `e0_realised` and
#'   which quantity was matched on.
#'
#' @examples
#' # A Coale-Demeny West table for males at level e0 = 60
#' m <- model_lifetable("CD_West", sex = "male", e0 = 60, graph = FALSE)
#' m$metrics$LifeExpectancyAtBirth
#' head(m$lifetable[, c("Age", "nMx", "nqx", "lx", "ex")])
#'
#' # Completing an indirect estimate: Manual X, Panama 1976.
#' # chm_brass() gives q(x) but no life table; the model table that matches
#' # its q(5) supplies the whole age pattern.
#' panama <- data.frame(
#'   women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
#'   ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
#'   cd    = c(40, 130, 312, 435, 636, 686, 689))
#' est <- chm_brass(panama, "women", "ceb", "cd", model = "west")
#' q5 <- est$qx[est$x == 5]
#' q5
#'
#' fitted <- model_lifetable("west", sex = "male", q5 = q5, graph = FALSE)
#' fitted
#'
#' # The same child mortality implies different adult mortality in different
#' # families, which is the whole point of choosing one:
#' sapply(c("CD_West", "CD_North", "CD_South", "CD_East"), function(f)
#'   model_lifetable(f, sex = "male", q5 = q5,
#'                   graph = FALSE)$metrics$LifeExpectancyAtBirth)
#'
#' @references
#' Coale, A. J., & Demeny, P. (1983). \emph{Regional Model Life Tables and Stable Populations} (2nd ed.). New York: Academic Press.
#'
#' United Nations (1982). \emph{Model Life Tables for Developing Countries}. Population Studies No. 77. New York: United Nations.
#'
#' United Nations (1983). \emph{Manual X: Indirect Techniques for Demographic Estimation}. Population Studies No. 81. New York: United Nations.
#'
#' @seealso [model_lt] for the tabulated rates, [chm_brass()] for the child
#'   mortality these tables complete, and [brass_lifetable()] for the
#'   relational alternative to picking a family.
#' @export
model_lifetable <- function(family = "CD_West",
                            sex = c("male", "female"),
                            e0 = NULL, q1 = NULL, q5 = NULL,
                            radix = 1e5, graph = TRUE) {
  sex <- match.arg(sex)
  family <- .mlt_family(family)

  given <- sum(!is.null(e0), !is.null(q1), !is.null(q5))
  if (given != 1) {
    stop("Supply exactly one of 'e0', 'q1' and 'q5'; ", given, " were given.")
  }

  tab <- model_lt[model_lt$family == family & model_lt$sex == sex, ]
  if (nrow(tab) == 0) stop("No packaged table for family '", family, "'.")
  levels_e0 <- sort(unique(tab$e0))

  matched_on <- "e0"
  if (!is.null(q1) || !is.null(q5)) {
    matched_on <- if (!is.null(q1)) "q1" else "q5"
    target <- if (!is.null(q1)) q1 else q5
    if (!is.numeric(target) || length(target) != 1 || !is.finite(target) ||
        target <= 0 || target >= 1) {
      stop("'", matched_on, "' must be a single number strictly between 0 and 1.")
    }
    at_age <- if (!is.null(q1)) 1 else 5
    # q(x) falls as the level rises, so the relation is monotone
    qs <- vapply(levels_e0, function(L) {
      lt <- .mlt_build(tab, L, levels_e0, sex, radix)
      1 - lt$lx[lt$Age == at_age] / radix
    }, numeric(1))
    if (target > max(qs) || target < min(qs)) {
      stop(sprintf(
        "%s = %.4f is outside what the '%s' tables cover (%.4f to %.4f, for e0 %g to %g). %s",
        matched_on, target, family, min(qs), max(qs),
        min(levels_e0), max(levels_e0),
        "Use brass_lifetable(), which is not limited to a tabulated range."))
    }
    # A linear inversion of q(e0) lands within a few parts in ten thousand;
    # since the rates are interpolated continuously between levels, the level
    # that reproduces the target exactly can be solved for.
    start <- stats::approx(qs, levels_e0, xout = target)$y
    f <- function(L) {
      lt <- .mlt_build(tab, L, levels_e0, sex, radix)
      (1 - lt$lx[lt$Age == at_age] / radix) - target
    }
    lo <- max(min(levels_e0), start - 1)
    hi <- min(max(levels_e0), start + 1)
    e0 <- if (is.finite(f(lo)) && is.finite(f(hi)) && f(lo) * f(hi) <= 0) {
      stats::uniroot(f, lower = lo, upper = hi, tol = 1e-12)$root
    } else {
      start
    }
  }

  if (!is.numeric(e0) || length(e0) != 1 || !is.finite(e0)) {
    stop("'e0' must be a single finite number.")
  }
  if (e0 < min(levels_e0) || e0 > max(levels_e0)) {
    stop(sprintf("e0 = %g is outside the tabulated range %g to %g.",
                 e0, min(levels_e0), max(levels_e0)))
  }

  lt <- .mlt_build(tab, e0, levels_e0, sex, radix)

  out <- list(
    metrics = list(TotalDeaths = sum(lt$dx),
                   LifeExpectancyAtBirth = lt$ex[1]),
    lifetable = lt,
    model = list(family = family, sex = sex, e0_requested = e0,
                 e0_realised = lt$ex[1], matched_on = matched_on)
  )
  if (graph) out$plot <- .plot_survival(lt$Age, lt$lx)
  class(out) <- c("dem_model_lt", "dem_lifetable")
  out
}

# Accept both the packaged names and the short Coale-Demeny names used by
# chm_brass(), so that an estimate can be carried straight across.
.mlt_family <- function(family) {
  if (!is.character(family) || length(family) != 1) {
    stop("'family' must be a single family name.")
  }
  short <- c(west = "CD_West", north = "CD_North",
             south = "CD_South", east = "CD_East")
  if (tolower(family) %in% names(short)) return(unname(short[tolower(family)]))
  known <- levels(model_lt$family)
  hit <- known[tolower(known) == tolower(family)]
  if (length(hit) == 1) return(hit)
  stop("Unknown family '", family, "'. Choose one of: ",
       paste(known, collapse = ", "), " (or west/north/south/east).")
}

# Rates at an arbitrary level, interpolated on log(nMx) between the two
# tabulated levels either side, then turned into a life table by lifetable().
.mlt_build <- function(tab, e0, levels_e0, sex, radix) {
  ages <- sort(unique(tab$age))
  rate_at <- function(L) tab$nMx[tab$e0 == L][order(tab$age[tab$e0 == L])]

  if (e0 %in% levels_e0) {
    mx <- rate_at(e0)
  } else {
    lo <- max(levels_e0[levels_e0 <= e0])
    hi <- min(levels_e0[levels_e0 >= e0])
    w <- (e0 - lo) / (hi - lo)
    mx <- exp((1 - w) * log(rate_at(lo)) + w * log(rate_at(hi)))
  }
  lifetable(data.frame(Age = ages, nMx = mx), age = "Age", nMx = "nMx",
            sex = sex, radix = radix, graph = FALSE)$lifetable
}

#' @export
print.dem_model_lt <- function(x, ...) {
  m <- x$model
  cat(sprintf("Model life table: %s, %s\n", m$family, m$sex))
  if (m$matched_on == "e0") {
    cat(sprintf("  level e0 = %g requested; table gives e0 = %.2f\n",
                m$e0_requested, m$e0_realised))
  } else {
    cat(sprintf("  matched on %s; level e0 = %.2f, table gives e0 = %.2f\n",
                m$matched_on, m$e0_requested, m$e0_realised))
  }
  cat("\n")
  n_show <- min(10L, nrow(x$lifetable))
  print(x$lifetable[seq_len(n_show), ], row.names = FALSE)
  if (nrow(x$lifetable) > n_show) {
    cat(sprintf("... (%d age groups)\n", nrow(x$lifetable)))
  }
  invisible(x)
}
