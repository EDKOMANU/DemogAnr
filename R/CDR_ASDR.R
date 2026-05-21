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
#'
#' @return A list of named outputs for the requested death rate calculations.
#' @examples
#' data <- data.frame(
#'   age = c("0-4", "5-9", "10-14"),
#'   population = c(10000, 12000, 11000),
#'   deaths = c(100, 50, 30)
#' )
#' dem.cdr(data, type = "all", age_col = "age",
#'                     population_col = "population", deaths_col = "deaths")
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapters 2 and 5)
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press. (Chapters 4 and 6)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). Emerald Group Publishing. (Chapters 14 for Maternal Mortality)
#'
#' @export
dem.cdr <- function(data, type, age_col = NULL, population_col, deaths_col) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (is.character(type) && type == "all") type <- c("CDR", "ASDR")
  if (!all(type %in% c("CDR", "ASDR"))) stop("Invalid 'type'. Choose from 'CDR' or 'ASDR'.")
  if (!all(c(population_col, deaths_col) %in% colnames(data))) stop("Specified columns not found in the dataframe.")
  if (any(type == "ASDR") && is.null(age_col)) stop("'age_col' must be specified for ASDR.")

  # Initialize result list
  results <- list()

  # Calculate CDR
  if ("CDR" %in% type) {
    cdr <- (sum(data[[deaths_col]]) / sum(data[[population_col]])) * 1000
    results$CDR <- cdr
  }

  # Calculate ASDR
  if ("ASDR" %in% type) {
    asdr <- (data[[deaths_col]] / data[[population_col]]) * 1000
    data$ASDR <- asdr  # Add ASDR to the dataframe
    results$ASDR <- asdr
  }

  out <- list(results = results, modified_data = data)
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



