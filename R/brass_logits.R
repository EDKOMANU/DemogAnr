#' Brass Relational Logit Model
#'
#' Fits the Brass relational logit model to an observed series of *cumulative*
#' probabilities of dying \eqn{q(x) = 1 - l(x)/l(0)} using a standard life
#' table. The observed logits are regressed on the logits of the standard, and
#' the fitted relationship (alpha, beta) is used to predict a complete,
#' smoothed series of \eqn{q(x)}, including for ages where the observed value
#' is missing (`NA`).
#'
#' The Brass logit of a survivorship value is
#' \eqn{Y(x) = 0.5 \log\{(1 - l(x))/l(x)\} = 0.5 \log\{q(x)/(1 - q(x))\}},
#' and predicted probabilities are recovered with the inverse transform
#' \eqn{q(x) = \exp(2Y)/(1 + \exp(2Y))}.
#'
#' Note that the input must be cumulative probabilities of dying by age
#' \eqn{x} (i.e. \eqn{1 - l(x)}), which increase with age, and not the
#' age-interval probabilities \eqn{{}_nq_x}. Values of exactly 0 or 1 (for
#' example \eqn{q(0) = 0}) are excluded from the fit.
#'
#' @param data A data frame containing the observed cumulative probabilities
#'   of dying \eqn{q(x)} and ages.
#' @param qx_col Column name for the cumulative probability of dying
#'   \eqn{q(x) = 1 - l(x)} in the data frame. May contain `NA` for ages to be
#'   estimated.
#' @param age_col Column name for age in the data frame.
#' @param standard Standard to use ("African" or "GeneralUN").
#' @param standards_data Standards dataset (default: included [standards]).
#' @param complete Logical. If `TRUE`, the fitted \eqn{\alpha} and
#'   \eqn{\beta} are carried straight through to the complete life table they
#'   imply, attached as `$lifetable`. The fit on its own predicts only at the
#'   ages supplied, so this is usually what is wanted; it is the same as
#'   passing the result to [brass_lifetable()]. Default `FALSE`.
#' @param sex,radix Passed to [brass_lifetable()] when `complete = TRUE`.
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console.
#' @return An object of class `dem_brass_fit`: a list with the fitted `model`,
#'   the `coefficients` (also given separately as `alpha` and `beta`), the
#'   `standard` used, and `data`: the input data with a `predicted_qx` column of
#'   smoothed/completed cumulative probabilities and a `predicted_lx` column of
#'   the implied survivorship \eqn{l(x) = 1 - q(x)}.
#'
#'   Note that the predictions cover only the ages supplied in `data`. To turn
#'   the fitted \eqn{\alpha} and \eqn{\beta} into a complete life table over
#'   the standard's whole age range, set `complete = TRUE`, which attaches it
#'   as `$lifetable`, or pass the result to [brass_lifetable()].
#' @importFrom stats lm predict coef
#'
#' @examples
#' # Observed cumulative probabilities of dying q(x) = 1 - l(x),
#' # with gaps (NA) to be filled by the model:
#' observed_data <- data.frame(
#'   age = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
#'   qx = c(0.10, 0.15, NA, 0.19, 0.21, NA, 0.26, 0.30, NA, 0.38, 0.44)
#' )
#'
#' model <- brass_logit(
#'   data = observed_data,
#'   qx_col = "qx",
#'   age_col = "age",
#'   standard = "African"
#' )
#'
#' model$coefficients
#' model$data
#'
#' # Or go straight to the life table the fit implies
#' full <- brass_logit(observed_data, qx_col = "qx", age_col = "age",
#'                     standard = "African", complete = TRUE)
#' full$lifetable$metrics$LifeExpectancyAtBirth
#'
#' @references
#' Brass, W. (1971). On the scale of mortality. In W. Brass (Ed.), \emph{Biological Aspects of Demography} (pp. 69-110). London: Taylor & Francis.
#'
#' Brass, W. (1975). \emph{Methods for Estimating Fertility and Mortality from Limited and Defective Data}. Chapel Hill: Carolina Population Center, University of North Carolina.
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: International Union for the Scientific Study of Population (IUSSP). \url{https://demographicestimation.iussp.org/}
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press. (Chapter on model life tables.)
#'
#' @export
brass_logit <- function(data, qx_col, age_col, standard, standards_data = NULL,
                        complete = FALSE, sex = c("male", "female"),
                        radix = 1e5, verbose = FALSE) {
  sex <- match.arg(sex)
  # Load the internal standards dataset if not provided
  if (is.null(standards_data)) {
    standards_data <- standards
  }

  # Validate inputs
  if (!is.data.frame(data)) stop("Input data must be a data frame.")
  if (!qx_col %in% colnames(data) || !age_col %in% colnames(data)) {
    stop("Specified columns for qx and age must exist in the data frame.")
  }
  if (!standard %in% c("African", "GeneralUN")) {
    stop("Invalid standard. Choose 'African' or 'GeneralUN'.")
  }

  # Extract observed cumulative probabilities and age
  qx <- data[[qx_col]]
  age <- data[[age_col]]

  if (!is.numeric(qx)) stop("The qx column must be numeric.")
  if (any(qx < 0 | qx > 1, na.rm = TRUE)) stop("qx values must be between 0 and 1.")

  # Warn if the observed series looks like interval nqx rather than cumulative q(x)
  non_na <- stats::na.omit(data.frame(age = age, qx = qx))
  if (nrow(non_na) >= 3 && all(diff(non_na$qx[order(non_na$age)]) <= 0)) {
    warning("Observed qx values decrease with age. brass_logit() expects *cumulative* probabilities q(x) = 1 - l(x), which increase with age.")
  }

  # Match age in standards data
  standard_logits <- standards_data[[paste0(standard, "_logit")]][match(age, standards_data$Age)]
  if (any(is.na(standard_logits))) {
    stop("Mismatch between age groups in the input data and standards data.")
  }

  # Brass (half) logit transformation of observed q(x); q of exactly 0 or 1
  # cannot be transformed and is excluded from the fit
  usable <- !is.na(qx) & qx > 0 & qx < 1
  logit_qx <- ifelse(usable, 0.5 * log(qx / (1 - qx)), NA_real_)

  if (sum(usable) < 2) {
    stop("At least two usable (non-NA, strictly between 0 and 1) qx values are required to fit the model.")
  }

  if (verbose) {
    message(sprintf("Fitting Brass Relational Logit Model using standard: '%s'...", standard))
    message(sprintf("  Number of observed age groups: %d (with %d usable values)", nrow(data), sum(usable)))
  }

  # Fit the model using usable values
  lm_fit <- lm(logit_qx ~ standard_logits, subset = usable)

  if (verbose) {
    message(sprintf("  Coefficients: Alpha (Intercept) = %.4f, Beta (Slope) = %.4f", coef(lm_fit)[1], coef(lm_fit)[2]))
  }

  # Predict logits for all ages
  predicted_logits <- predict(lm_fit, newdata = data.frame(standard_logits = standard_logits))

  # Inverse Brass logit: q(x) = exp(2Y) / (1 + exp(2Y))
  predicted_qx <- 1 / (1 + exp(-2 * predicted_logits))

  data$predicted_qx <- predicted_qx
  data$predicted_lx <- 1 - predicted_qx

  # Return results. 'standard' is recorded so that a fit can be handed
  # straight to brass_lifetable() without naming the standard again.
  out <- list(
    model = lm_fit,
    coefficients = coef(lm_fit),
    alpha = unname(coef(lm_fit)[1]),
    beta = unname(coef(lm_fit)[2]),
    standard = standard,
    standards_data = standards_data,
    data = data
  )
  class(out) <- "dem_brass_fit"

  # The fit predicts only at the ages supplied. complete = TRUE carries it
  # through to the life table the parameters imply over the standard's whole
  # age range, so the common path is a single call.
  if (isTRUE(complete)) {
    out$lifetable <- brass_lifetable(out, sex = sex, radix = radix,
                                     graph = FALSE)
    if (verbose) {
      message(sprintf("  Completed to a life table: e0 = %.2f",
                      out$lifetable$metrics$LifeExpectancyAtBirth))
    }
  }
  out
}


#' Life Table from Brass Relational Logit Parameters
#'
#' Generates a complete life table from a pair of Brass relational logit
#' parameters, \eqn{\alpha} and \eqn{\beta}, and a standard. Where
#' [brass_logit()] estimates \eqn{\alpha} and \eqn{\beta} from observed
#' survivorship and predicts only at the ages supplied, `brass_lifetable()`
#' runs the model in the other direction: it applies the fitted relation across
#' the standard's whole age range and returns the life-table columns, including
#' life expectancy.
#'
#' @details
#' The Brass relation transforms one survivorship curve into another through a
#' linear relation between their logits,
#' \deqn{Y(x) = \alpha + \beta\,Y^s(x), \qquad
#'   Y(x) = \tfrac{1}{2}\log\frac{1 - l(x)}{l(x)},}
#' so that survivorship is recovered as
#' \eqn{l(x) = 1 - \exp(2Y)/(1 + \exp(2Y))}. The level parameter \eqn{\alpha}
#' shifts mortality up or down at every age (\eqn{\alpha < 0} gives lighter
#' mortality than the standard), while \eqn{\beta} tilts the age pattern:
#' \eqn{\beta < 1} raises child mortality relative to adult mortality, and
#' \eqn{\beta > 1} does the reverse. Setting \eqn{\alpha = 0, \beta = 1}
#' reproduces the standard exactly.
#'
#' Survivorship at age 0 is 1 by definition and is not obtained from the
#' transform (the standard tabulates a placeholder logit there). From
#' \eqn{l(x)} the remaining columns follow exactly:
#' \eqn{{}_nd_x = l_x - l_{x+n}} and \eqn{{}_nq_x = {}_nd_x / l_x} need no
#' assumption, while \eqn{{}_nL_x = n\,l_{x+n} + {}_na_x\,{}_nd_x} uses
#' \eqn{{}_na_x = n/2} for closed intervals, with the Coale-Demeny separation
#' factor at age 0 (and for 1-4 where the standard is grouped that way). The
#' factor depends on \eqn{{}_1M_0}, which is itself solved for from
#' \eqn{{}_1q_0} by a short iteration.
#'
#' For the open-ended interval the death rate is extrapolated by continuing the
#' log-linear (Gompertz) trend of \eqn{{}_nM_x} through the last two closed
#' intervals, or set directly with `open_mx`. Because survivorship at the start
#' of that interval is minute -- about 0.0006 of the radix for the packaged
#' standards, which end at age 95 -- the choice moves \eqn{e_0} by well under a
#' hundredth of a year.
#'
#' @param alpha Either the level parameter \eqn{\alpha}, or a fitted object
#'   returned by [brass_logit()], in which case `beta` and `standard` are taken
#'   from it. Omit it when using `target_e0` or `target_q5`.
#' @param beta The slope parameter \eqn{\beta}, which must be positive
#'   (default 1). Ignored when `alpha` is a fitted object.
#' @param standard Standard to use, `"African"` or `"GeneralUN"`. Ignored when
#'   `alpha` is a fitted object.
#' @param sex Character: `"male"` (default) or `"female"`. Selects the
#'   Coale-Demeny separation factor at age 0.
#' @param radix Numeric: survivors at exact age 0 (default 100,000).
#' @param standards_data Standards dataset (default: the included [standards]).
#' @param open_mx Optional death rate for the open-ended interval. When `NULL`
#'   (default) it is extrapolated from the last two closed intervals.
#' @param target_e0,target_q5 Index the table by level instead of by
#'   \eqn{\alpha}: give a target life expectancy at birth, or a target
#'   under-five mortality \eqn{q(5)}, and the \eqn{\alpha} that attains it is
#'   solved for with `beta` held fixed. Supply exactly one of `alpha`,
#'   `target_e0` and `target_q5`.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} survival curve is
#'   attached to the result as `$plot`; retrieve it with `plot()`.
#'
#' @return An object of class `dem_brass_lt`, which inherits from
#'   `dem_lifetable`: a list with `metrics` (life expectancy at birth and total
#'   life-table deaths), the `lifetable` data frame, and `brass`, recording the
#'   `alpha`, `beta` and `standard` used.
#'
#' @examples
#' # Indirect child mortality, then the life table it implies.
#' # United Nations Manual X worked example: Panama, 1976.
#' panama <- data.frame(
#'   women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
#'   ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
#'   cd    = c(40, 130, 312, 435, 636, 686, 689)
#' )
#' est <- chm_brass(panama, women_col = "women", ceb_col = "ceb",
#'                  cd_col = "cd", model = "west")
#'
#' # Fit the relational model to the estimated q(x) ...
#' fit <- brass_logit(data.frame(age = est$x, qx = est$qx),
#'                    qx_col = "qx", age_col = "age", standard = "African")
#' fit$coefficients
#'
#' # ... and expand it into a complete life table
#' lt <- brass_lifetable(fit)
#' lt$metrics$LifeExpectancyAtBirth
#' head(lt$lifetable)
#'
#' # alpha = 0, beta = 1 returns the standard itself
#' std <- brass_lifetable(alpha = 0, beta = 1, standard = "African",
#'                        graph = FALSE)
#' head(std$lifetable[, c("Age", "lx", "nqx", "ex")])
#'
#' # A lighter-mortality variant of the same age pattern
#' brass_lifetable(alpha = -0.4, beta = 1,
#'                 standard = "African", graph = FALSE)$metrics
#'
#' # Index by level rather than by alpha: the table with e0 = 60
#' m <- brass_lifetable(target_e0 = 60, standard = "African", graph = FALSE)
#' m$brass$alpha
#' m$metrics$LifeExpectancyAtBirth
#'
#' @references
#' Brass, W. (1971). On the scale of mortality. In W. Brass (Ed.), \emph{Biological Aspects of Demography} (pp. 69-110). London: Taylor & Francis.
#'
#' Brass, W. (1975). \emph{Methods for Estimating Fertility and Mortality from Limited and Defective Data}. Chapel Hill: Carolina Population Center, University of North Carolina.
#'
#' Moultrie, T., Dorrington, R., Hill, A., Hill, K., Timaeus, I., & Zaba, B. (Eds.). (2013). \emph{Tools for Demographic Estimation}. Paris: IUSSP. \url{https://demographicestimation.iussp.org/}
#'
#' Coale, A. J., & Demeny, P. (1983). \emph{Regional Model Life Tables and Stable Populations} (2nd ed.). New York: Academic Press.
#'
#' @seealso [brass_logit()], which estimates `alpha` and `beta`; [lifetable()]
#'   for a life table built from observed death rates.
#' @export
brass_lifetable <- function(alpha = NULL, beta = 1,
                            standard = c("African", "GeneralUN"),
                            sex = c("male", "female"),
                            radix = 1e5, standards_data = NULL,
                            open_mx = NULL, target_e0 = NULL, target_q5 = NULL,
                            graph = TRUE) {
  sex <- match.arg(sex)

  n_given <- sum(!is.null(alpha), !is.null(target_e0), !is.null(target_q5))
  if (n_given != 1) {
    stop("Supply exactly one of 'alpha' (or a brass_logit() fit), 'target_e0', ",
         "or 'target_q5'; ", n_given, " were given.")
  }

  # A fitted object carries alpha, beta and the standard it was fitted against
  if (inherits(alpha, "dem_brass_fit")) {
    fit <- alpha
    beta <- fit$beta
    standard <- fit$standard
    if (is.null(standards_data)) standards_data <- fit$standards_data
    alpha <- fit$alpha
  } else {
    standard <- match.arg(standard)
  }
  if (is.null(standards_data)) standards_data <- standards

  if (!is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
    stop("'beta' must be a single positive number; got ", format(beta), ".")
  }

  # Indexing by level: solve for the alpha that hits a target, holding beta
  # (the age pattern) fixed. Both targets are monotone in alpha, so a single
  # root exists wherever the target is attainable from this standard.
  if (!is.null(target_e0) || !is.null(target_q5)) {
    alpha <- .brass_solve_alpha(target_e0, target_q5, beta, standard, sex,
                                radix, standards_data, open_mx)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha)) {
    stop("'alpha' must be a single finite number, or a brass_logit() fit.")
  }
  col <- paste0(standard, "_logit")
  if (!col %in% names(standards_data)) {
    stop("Standard '", standard, "' not found in the standards data.")
  }

  age <- as.numeric(standards_data$Age)
  ord <- order(age); age <- age[ord]
  Ys <- as.numeric(standards_data[[col]])[ord]
  k <- length(age)
  n <- c(diff(age), Inf)

  # l(0) = 1 by definition: the standard tabulates a placeholder logit at
  # age 0, since the logit of l = 1 is not finite.
  lprop <- numeric(k)
  Y <- alpha + beta * Ys
  lprop <- 1 - 1 / (1 + exp(-2 * Y))
  if (age[1] == 0) lprop[1] <- 1

  lx <- radix * lprop
  dx <- c(-diff(lx), lx[k])
  nqx <- dx / lx
  nqx[k] <- 1

  # nax: n/2, with the Coale-Demeny factor at age 0. 1M0 depends on a0 and a0
  # on 1M0, so solve the pair by a short iteration.
  nax <- n / 2
  if (age[1] == 0) {
    a0 <- 0.1
    for (i in 1:20) {
      L0 <- n[1] * lx[2] + a0 * dx[1]
      m0 <- dx[1] / L0
      new <- .cd_a0_a1(m0, sex)[["a0"]]
      if (abs(new - a0) < 1e-10) { a0 <- new; break }
      a0 <- new
    }
    nax[1] <- a0
    # a standard grouped as 0, 1-4, 5, ... also takes the CD factor at 1-4
    if (k >= 2 && age[2] == 1 && n[2] == 4) {
      nax[2] <- .cd_a0_a1(m0, sex)[["a1"]]
    }
  }

  Lx <- numeric(k)
  for (i in seq_len(k - 1)) Lx[i] <- n[i] * lx[i + 1] + nax[i] * dx[i]
  nMx <- c(dx[-k] / Lx[-k], NA_real_)

  # Open interval: continue the log-linear trend in nMx, or take open_mx
  mx_open <- if (!is.null(open_mx)) {
    as.numeric(open_mx)[1]
  } else if (k >= 3 && all(is.finite(nMx[c(k - 2, k - 1)])) &&
             all(nMx[c(k - 2, k - 1)] > 0)) {
    nMx[k - 1] * (nMx[k - 1] / nMx[k - 2])
  } else if (is.finite(nMx[k - 1]) && nMx[k - 1] > 0) {
    nMx[k - 1]
  } else {
    NA_real_
  }
  if (!is.finite(mx_open) || mx_open <= 0) mx_open <- 1
  nMx[k] <- mx_open
  nax[k] <- 1 / mx_open
  Lx[k] <- lx[k] / mx_open

  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx

  lt <- data.frame(Age = age, n = n, nMx = nMx, nax = nax, nqx = nqx,
                   lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex)

  out <- list(
    metrics = list(TotalDeaths = sum(dx), LifeExpectancyAtBirth = ex[1]),
    lifetable = lt,
    brass = list(alpha = alpha, beta = beta, standard = standard, sex = sex)
  )
  if (graph) out$plot <- .plot_survival(lt$Age, lt$lx)
  class(out) <- c("dem_brass_lt", "dem_lifetable")
  out
}

# Solve for the alpha that makes a Brass table hit a target, with beta held
# fixed. Life expectancy falls, and q(5) rises, monotonically with alpha, so
# uniroot has a single root wherever the target is reachable.
.brass_solve_alpha <- function(target_e0, target_q5, beta, standard, sex,
                               radix, standards_data, open_mx) {
  build <- function(a) {
    brass_lifetable(alpha = a, beta = beta, standard = standard, sex = sex,
                    radix = radix, standards_data = standards_data,
                    open_mx = open_mx, graph = FALSE)
  }
  if (!is.null(target_e0)) {
    if (!is.numeric(target_e0) || length(target_e0) != 1 ||
        !is.finite(target_e0) || target_e0 <= 0) {
      stop("'target_e0' must be a single positive number.")
    }
    f <- function(a) build(a)$metrics$LifeExpectancyAtBirth - target_e0
    what <- sprintf("a life expectancy at birth of %.2f", target_e0)
  } else {
    if (!is.numeric(target_q5) || length(target_q5) != 1 ||
        !is.finite(target_q5) || target_q5 <= 0 || target_q5 >= 1) {
      stop("'target_q5' must be a single number strictly between 0 and 1.")
    }
    f <- function(a) {
      lt <- build(a)$lifetable
      (1 - lt$lx[lt$Age == 5] / lt$lx[1]) - target_q5
    }
    what <- sprintf("an under-five mortality q(5) of %.4f", target_q5)
  }
  lo <- -5; hi <- 5
  flo <- f(lo); fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
    stop("Cannot reach ", what, " from the '", standard, "' standard with ",
         "beta = ", format(beta), ". The target lies outside the range this ",
         "standard can produce; try another standard or another beta.")
  }
  stats::uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

#' @export
print.dem_brass_lt <- function(x, ...) {
  b <- x$brass
  cat(sprintf("Brass relational life table (%s standard, %s)\n",
              b$standard, b$sex))
  cat(sprintf("  alpha = %.4f (level)   beta = %.4f (age pattern)\n",
              b$alpha, b$beta))
  cat(sprintf("  e0 = %.2f years\n\n", x$metrics$LifeExpectancyAtBirth))
  n_show <- min(10L, nrow(x$lifetable))
  print(x$lifetable[seq_len(n_show), ], row.names = FALSE)
  if (nrow(x$lifetable) > n_show) {
    cat(sprintf("... (%d age groups)\n", nrow(x$lifetable)))
  }
  invisible(x)
}
