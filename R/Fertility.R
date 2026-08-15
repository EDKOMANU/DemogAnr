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
#' @param age_interval Width (in years) of the age groups, used to compute TFR
#'   from the ASFRs (default 5, for standard 5-year age groups).
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} plot of the
#'   age-specific fertility rates is attached to the result as `$plot` and shown
#'   when the object is printed.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#' @return A list of named outputs for the requested fertility calculations and the modified dataframe with ASFR.
#' @examples
#' # United Nations (1983) Manual X, Ch. II: Bangladesh, 1974.
#' # Women and births in the year preceding the survey, by age of mother.
#' # The reported total fertility rate of this schedule is 4.83.
#' bangladesh <- data.frame(
#'   age         = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
#'   women       = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
#'   live_births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125),
#'   population  = 71315944
#' )
#'
#' # Age-specific fertility rates
#' dem.fert(bangladesh, type = "ASFR", age_col = "age",
#'          population_col = "population", women_col = "women",
#'          births_col = "live_births", graph = FALSE)
#'
#' # All fertility measures (TFR = 4.83 as published in Manual X)
#' dem.fert(bangladesh, type = "all", age_col = "age",
#'          population_col = "population", women_col = "women",
#'          births_col = "live_births", graph = FALSE)
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapter 5: Fertility and Reproduction.)
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press.
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Chapter on fertility measures.)
#'
#' @export
dem.fert <- function(data, type, age_col = NULL, population_col, women_col = NULL, births_col, age_interval = 5, graph = TRUE, verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (identical(type, "all")) type <- c("CBR", "GFR", "ASFR", "TFR")
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
    # TFR is the sum of ASFRs times the width of the age groups
    if (!"ASFR" %in% names(results)) {
      asfr <- (data[[births_col]] / data[[women_col]]) * 1000
      results$ASFR <- asfr
    }
    tfr <- sum(results$ASFR, na.rm = TRUE) * age_interval / 1000
    results$TFR <- tfr
    if (verbose) message("dem.fert: TFR calculated successfully: ", round(tfr, 4))
  }

  out <- list(results = results, modified_data = data)
  if (graph && "ASFR" %in% names(results) && !is.null(age_col)) {
    out$plot <- .plot_series(data[[age_col]], results$ASFR,
                             "ASFR (per 1,000 women)",
                             "Age-specific fertility rates")
  }
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
  if (!is.null(x$plot)) print(x$plot)
  invisible(x)
}

#' @export
plot.dem_fert <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call dem.fert(..., graph = TRUE).")
  x$plot
}
