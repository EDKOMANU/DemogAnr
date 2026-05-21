#' Calculate Fertility Metrics
#'
#' This function calculates various fertility metrics (CBR, GFR, ASFR, TFR) based on user input
#' and a dataframe with demographic data.
#'
#' @param data A dataframe containing demographic data.
#' @param type A string or vector specifying the type(s) of fertility calculation.
#'   Options: "CBR", "GFR", "ASFR", "TFR". Use "all" to calculate all metrics.
#' @param age_col The column name for age groups (used for "ASFR" and "TFR").
#' @param population_col The column name for total population.
#' @param women_col The column name for number of women (used for "GFR", "ASFR", and "TFR").
#' @param births_col The column name for live births.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#' @return A list of named outputs for the requested fertility calculations and the modified dataframe with ASFR.
#' @examples
#' demo_data <- data.frame(
#'   age = c("15-19", "20-24", "25-29"),
#'   population = c(10000, 12000, 11000),
#'   women = c(5000, 6000, 5500),
#'   live_births = c(200, 400, 300)
#' )
#
#' # Calculate ASFR
#' dem.fert(demo_data, type = "ASFR", age_col = "age",
#'                     population_col = "population", women_col = "women", births_col = "live_births")
#' # Calculate all metrics
#' dem.fert(demo_data, type = "all", age_col = "age",
#'                     population_col = "population", women_col = "women", births_col = "live_births")
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapters 2 and 5)
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press. (Chapters 4 and 6)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). Emerald Group Publishing. (Chapters 14 for Maternal Mortality)
#'
#' @export
dem.fert <- function(data, type, age_col = NULL, population_col, women_col = NULL, births_col, verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (is.character(type) && type == "all") type <- c("CBR", "GFR", "ASFR", "TFR")
  if (!all(type %in% c("CBR", "GFR", "ASFR", "TFR"))) stop("Invalid 'type'. Choose from 'CBR', 'GFR', 'ASFR', or 'TFR'.")
  if (!all(c(population_col, births_col) %in% colnames(data))) stop("Specified columns not found in the dataframe.")
  if (any(type %in% c("GFR", "ASFR", "TFR")) && is.null(women_col)) stop("'women_col' must be specified for GFR, ASFR, or TFR.")
  if (any(type %in% c("ASFR", "TFR")) && is.null(age_col)) stop("'age_col' must be specified for ASFR or TFR.")

  if (verbose) {
    message("dem.fert: Starting fertility calculations for type(s): ", paste(type, collapse = ", "))
  }

  # Initialize result list
  results <- list()

  # Calculate metrics based on type
  if ("CBR" %in% type) {
    if (verbose) message("dem.fert: Calculating Crude Birth Rate (CBR)...")
    cbr <- (sum(data[[births_col]]) / sum(data[[population_col]])) * 1000
    results$CBR <- cbr
    if (verbose) message("dem.fert: CBR calculated successfully: ", round(cbr, 4))
  }
  if ("GFR" %in% type) {
    if (verbose) message("dem.fert: Calculating General Fertility Rate (GFR)...")
    gfr <- (sum(data[[births_col]]) / sum(data[[women_col]])) * 1000
    results$GFR <- gfr
    if (verbose) message("dem.fert: GFR calculated successfully: ", round(gfr, 4))
  }
  if ("ASFR" %in% type || "TFR" %in% type) {
    if (verbose) message("dem.fert: Calculating Age-Specific Fertility Rate (ASFR)...")
    asfr <- (data[[births_col]] / data[[women_col]]) * 1000
    data$ASFR <- asfr  # Add ASFR to the dataframe
    results$ASFR <- asfr
    if (verbose) message("dem.fert: ASFR calculated for ", length(asfr), " age groups.")
  }
  if ("TFR" %in% type) {
    if (verbose) message("dem.fert: Calculating Total Fertility Rate (TFR)...")
    # TFR is the sum of ASFRs, assuming age intervals of 5 years
    if (!"ASFR" %in% names(results)) {
      asfr <- (data[[births_col]] / data[[women_col]]) * 1000
      results$ASFR <- asfr
    }
    tfr <- sum(results$ASFR, na.rm = TRUE) * 5 / 1000
    results$TFR <- tfr
    if (verbose) message("dem.fert: TFR calculated successfully: ", round(tfr, 4))
  }

  out <- list(results = results, modified_data = data)
  class(out) <- "dem_fert"
  return(out)
}

#' @export
print.dem_fert <- function(x, ...) {
  if ("CBR" %in% names(x$results)) {
    cat("Crude Birth Rate (CBR):\n")
    print(x$results$CBR)
  }
  if ("GFR" %in% names(x$results)) {
    cat("General Fertility Rate (GFR):\n")
    print(x$results$GFR)
  }
  if ("ASFR" %in% names(x$results)) {
    cat("Age-Specific Fertility Rate (ASFR):\n")
    print(x$results$ASFR)
  }
  if ("TFR" %in% names(x$results)) {
    cat("Total Fertility Rate (TFR):\n")
    print(x$results$TFR)
  }
  invisible(x)
}
