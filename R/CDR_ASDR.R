#' Calculate CDR and ASDR
#'
#' This function calculates the Crude Death Rate (CDR) and Age-Specific Death Rates (ASDR).
#'
#' @param data A data frame containing demographic data.
#' @param type A string or vector specifying the type(s) of death rate calculation.
#'   Options: "CDR", "ASDR". Use "all" to calculate both.
#' @param age_col The column name for age groups (used for "ASDR").
#' @param population_col The column name for total population.
#' @param deaths_col The column name for total deaths.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   age-specific death rates is attached to the result as `$plot`;
#'   retrieve it with `plot()`.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#'
#' @return A list of named outputs for the requested death rate calculations.
#' @examples
#' # Preston et al. (2001) Box 1.2: Sweden, 1988. With 96,756 deaths and
#' # 8,438,477 person-years lived, the crude death rate is 11.47 per 1000.
#' sweden88 <- data.frame(population = 8438477, deaths = 96756)
#' dem.cdr(sweden88, type = "CDR", population_col = "population",
#'         deaths_col = "deaths", graph = FALSE)
#'
#' # Age-specific death rates: Ghana, 2010 Population and Housing Census
#' data(gphc2010)
#' dem.cdr(gphc2010, type = "all", age_col = "Age",
#'         population_col = "Pop", deaths_col = "Deaths", graph = FALSE)
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 2: Age-Specific Rates and Probabilities.)
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press.
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Chapter on mortality measures.)
#'
#' @export
dem.cdr <- function(data, type, age_col = NULL, population_col, deaths_col, graph = TRUE, verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (identical(type, "all")) type <- c("CDR", "ASDR")
  if (!all(type %in% c("CDR", "ASDR"))) stop("Invalid 'type'. Choose from 'CDR' or 'ASDR'.")
  if (!all(c(population_col, deaths_col) %in% colnames(data))) stop("Specified columns not found in the dataframe.")
  if (any(type == "ASDR") && is.null(age_col)) stop("'age_col' must be specified for ASDR.")

  if (verbose) {
    message("dem.cdr: Starting death rate calculations for type(s): ", paste(type, collapse = ", "))
  }

  # Initialize result list
  results <- list()

  # Calculate CDR
  if ("CDR" %in% type) {
    if (verbose) message("dem.cdr: Calculating Crude Death Rate (CDR)...")
    cdr <- (sum(data[[deaths_col]]) / sum(data[[population_col]])) * 1000
    results$CDR <- cdr
    if (verbose) message("dem.cdr: CDR calculated successfully: ", round(cdr, 4))
  }

  # Calculate ASDR
  if ("ASDR" %in% type) {
    if (verbose) message("dem.cdr: Calculating Age-Specific Death Rate (ASDR)...")
    asdr <- (data[[deaths_col]] / data[[population_col]]) * 1000
    data$ASDR <- asdr  # Add ASDR to the dataframe
    results$ASDR <- asdr
    if (verbose) message("dem.cdr: ASDR calculated successfully for ", length(asdr), " age groups.")
  }

  out <- list(results = results, modified_data = data)
  if (graph && "ASDR" %in% names(results) && !is.null(age_col)) {
    out$plot <- .plot_series(data[[age_col]], results$ASDR,
                             "ASDR (per 1,000)", "Age-specific death rates")
  }
  class(out) <- "dem_cdr"
  return(out)
}

#' @export
print.dem_cdr <- function(x, ...) {
  if ("CDR" %in% names(x$results)) {
    cat("Crude Death Rate (CDR):\n")
    print(x$results$CDR)
  }
  if ("ASDR" %in% names(x$results)) {
    cat("Age-Specific Death Rate (ASDR):\n")
    print(x$results$ASDR)
  }
  invisible(x)
}

#' @export
plot.dem_cdr <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call dem.cdr(..., graph = TRUE).")
  x$plot
}



