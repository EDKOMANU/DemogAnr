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
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console.
#' @return A list containing the fitted model (`model`), the alpha/beta
#'   `coefficients`, and `data`: the input data with a `predicted_qx` column of
#'   smoothed/completed cumulative probabilities and a `predicted_lx` column of
#'   the implied survivorship \eqn{l(x) = 1 - q(x)}.
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
brass_logit <- function(data, qx_col, age_col, standard, standards_data = NULL, verbose = FALSE) {
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

  # Return results
  return(list(
    model = lm_fit,
    coefficients = coef(lm_fit),
    data = data
  ))
}
