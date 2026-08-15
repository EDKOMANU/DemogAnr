#' Smart Interpolation Function
#'
#' Performs interpolation on numeric data, supporting linear, quadratic,
#' Lagrange, and cubic spline methods. Can handle both vector and dataframe inputs.
#'
#' @param x Numeric vector of x-values.
#' @param y Numeric vector of y-values (ignored if using a dataframe).
#' @param x_new Numeric vector of new x-values to interpolate.
#' @param data Optional dataframe containing y-values to interpolate.
#' @param y_cols Character vector of column names in `data` containing y-values.
#' @param method Character. Interpolation method: `"auto"`, `"linear"`, `"quadratic"`, `"lagrange"`, or `"spline"`.
#' Default is `"auto"`: linear for 2 points, quadratic for 3, and a natural
#' cubic spline for more than 3 points. High-degree Lagrange polynomials are
#' numerically unstable for many points (Runge phenomenon), so `"lagrange"` is
#' only used when requested explicitly.
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during interpolation.
#'
#' @return If vectors `x` and `y` are provided, returns a named numeric vector of interpolated values.
#' If a dataframe is provided, returns the original dataframe with additional interpolated columns.
#'
#' @examples
#' # Intercensal estimation, the commonest use of interpolation in demography:
#' # estimating the population in years between two censuses. Ghana enumerated
#' # 18,912,079 in 2000 and 30,079,802 in 2010 (the bundled `region2000` and
#' # `gphc2010` totals).
#' census_years <- c(2000, 2010)
#' census_pop   <- c(18912079, 30079802)
#' interpolation(census_years, census_pop, x_new = c(2003, 2005, 2008),
#'               method = "linear")
#'
#' # The data-frame interface interpolates each row across the columns given in
#' # `y_cols`, which is what intercensal estimation for several areas looks
#' # like: one row per region, one column per census.
#' regions <- data.frame(
#'   Region  = c("Greater Accra", "Ashanti", "Northern"),
#'   pop2000 = c(2905726, 3612950, 1087606),
#'   pop2010 = c(4010054, 4780380, 1387744)
#' )
#' interpolation(x = c(2000, 2010), data = regions,
#'               y_cols = c("pop2000", "pop2010"),
#'               x_new = c(2005, 2008), method = "linear")
#'
#' @references
#' Shryock, H. S., Siegel, J. S., & Larmon, E. A. (1973). \emph{The Methods and Materials of Demography}. Washington, DC: US Bureau of the Census. (Appendix C: Interpolation and graduation.)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Appendix C: Selected General Methods.)
#'
#' Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). \emph{Numerical Recipes: The Art of Scientific Computing} (3rd ed.). Cambridge: Cambridge University Press. ISBN 978-0521880688. (Chapter 3: Interpolation and Extrapolation.)
#'
#' @export
interpolation <- function(x, y, x_new,
                                data = NULL, y_cols = NULL,
                                method = c("auto", "linear", "quadratic", "lagrange", "spline"),
                                verbose = FALSE) {
  method <- match.arg(method)  # Ensure method is valid

  # Validate input
  if (!is.null(data) && !is.null(y_cols)) {
    if (!all(y_cols %in% names(data))) stop("Some y_cols not found in data")
    mode <- "dataframe"
  } else if (is.numeric(x) && is.numeric(y)) {
    if (length(x) != length(y)) stop("x and y must have the same length")
    mode <- "vector"
  } else {
    stop("Invalid input: Provide either (x, y) vectors or (data, y_cols) with x")
  }

  if (!is.numeric(x) || !is.numeric(x_new)) stop("x and x_new must be numeric")
  if (any(duplicated(x))) stop("x contains duplicate values, which is not allowed")
  if (method == "quadratic" && length(x) != 3) stop("Quadratic interpolation requires exactly 3 points.")
  if (length(x) < 2) stop("At least 2 points are required for interpolation.")

  if (verbose) {
    message(sprintf("Smart interpolation initialized (mode: '%s', method: '%s').", mode, method))
    message(sprintf("  Number of input points: %d", length(x)))
    message(sprintf("  Number of new points to interpolate: %d", length(x_new)))
  }

  # Warn if extrapolating
  if (any(x_new < min(x)) || any(x_new > max(x))) {
    warning("Some x_new values are outside the range of x. Extrapolation may be unreliable.")
  }

  # Interpolation Methods
  linear_interpolate <- function(x0, x1, y0, y1, x) {
    y0 + (x - x0) * ((y1 - y0)/(x1 - x0))
  }

  quadratic_interpolate <- function(x0, x1, x2, y0, y1, y2, x) {
    L0 <- ((x - x1)*(x - x2))/((x0 - x1)*(x0 - x2))
    L1 <- ((x - x0)*(x - x2))/((x1 - x0)*(x1 - x2))
    L2 <- ((x - x0)*(x - x1))/((x2 - x0)*(x2 - x1))
    y0*L0 + y1*L1 + y2*L2
  }

  lagrange_interpolate <- function(x_points, y_points, x_val) {
    n <- length(x_points)
    y_new <- 0
    for(i in 1:n) {
      term <- y_points[i]
      for(j in 1:n) {
        if(i != j) {
          term <- term * (x_val - x_points[j])/(x_points[i] - x_points[j])
        }
      }
      y_new <- y_new + term
    }
    return(y_new)
  }

  # Auto-method selection
  interpolate_points <- function(x_points, y_points, x_val) {
    n_points <- length(x_points)

    if (method == "linear" || (method == "auto" && n_points == 2)) {
      return(approx(x_points, y_points, xout = x_val, rule = 2)$y) # rule=2 allows extrapolation
    } else if (method == "quadratic" || (method == "auto" && n_points == 3)) {
      return(quadratic_interpolate(x_points[1], x_points[2], x_points[3],
                                   y_points[1], y_points[2], y_points[3],
                                   x_val))
    } else if (method == "lagrange") {
      return(lagrange_interpolate(x_points, y_points, x_val))
    } else if (method == "spline" || (method == "auto" && n_points > 3)) {
      return(spline(x_points, y_points, xout = x_val, method = "natural")$y)
    }
  }

  # Handle different input modes
  if (mode == "vector") {
    results <- sapply(x_new, function(x_val) interpolate_points(x, y, x_val))
    names(results) <- paste0("y_", x_new)
    return(results)
  } else {
    # Dataframe interpolation (without dplyr)
    y_list <- lapply(x_new, function(x_val) {
      apply(data[y_cols], 1, function(y_vals) interpolate_points(x, as.numeric(y_vals), x_val))
    })

    y_df <- do.call(cbind, y_list)
    colnames(y_df) <- paste0("y_", x_new)

    return(cbind(data, y_df))
  }
}
