#' Calculate Child Mortality Metrics by Age
#'
#' This function calculates child mortality metrics (neonatal, infant, child, under-5)
#' based on single-age data. It requires columns for age, denominator (live
#' births or population), and deaths.
#'
#' The result is expressed per 1,000 of the supplied denominator. If the
#' denominator is live births (as is conventional for neonatal and infant
#' mortality), the result is a rate per 1,000 live births. If a mid-year
#' population is supplied instead (often the only option for the child and
#' under-5 groups from census data), the result is a death *rate* per 1,000
#' population, which is not identical to the cohort probability of dying
#' (e.g. 5q0) reported by survey programmes such as the DHS.
#'
#' @param data A dataframe containing demographic data.
#' @param age_col The column name representing the age (in years or months for neonates).
#' @param live_births The column name for the denominator at each age
#'   (live births, or mid-year population if births are unavailable).
#' @param deaths_col The column name representing deaths at each age.
#' @param type The type of mortality calculation to perform:
#'             "neonatal", "infant", "child", or "under5".
#' @param age_in_months (Optional) The column name for age in months for neonatal mortality.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#' @return A numeric value representing the mortality rate per 1,000 live births.
#' @examples
#' # Ghana, 2021: deaths and population of children by age and sex.
#' # Where births are not registered, the population under age 1 is commonly
#' # used in place of live births as the denominator, as it is here.
#' data(ghmort2021)
#' males <- subset(ghmort2021, Sex == "Male")
#'
#' # Infant mortality (deaths under age 1)
#' dm.chm(males, age_col = "Age", live_births = "Population",
#'        deaths_col = "Deaths", type = "infant")
#'
#' # Under-five mortality
#' dm.chm(males, age_col = "Age", live_births = "Population",
#'        deaths_col = "Deaths", type = "under5")
#' @references
#' Pollard, A. H., Yusuf, F., & Pollard, G. N. (1990). \emph{Demographic Techniques} (3rd ed.). Sydney: Pergamon Press.
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161. (Chapters 1-2: rates versus probabilities of dying.)
#'
#' @export
dm.chm <- function(data, age_col, live_births, deaths_col,
                                             type = c("neonatal", "infant", "child", "under5"),
                                             age_in_months = NULL, verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")

  # Check required columns
  required_cols <- c(age_col, live_births, deaths_col)
  if (!all(required_cols %in% colnames(data))) stop("Required columns not found in the dataframe.")

  # Match type argument
  type <- match.arg(type)

  if (verbose) {
    message("dm.chm: Starting child mortality calculations for type: '", type, "'")
  }

  # Additional validation for neonatal mortality
  if (type == "neonatal" && is.null(age_in_months)) {
    stop("For neonatal mortality, the 'age_in_months' argument must be provided.")
  }
  if (!is.null(age_in_months) && !(age_in_months %in% colnames(data))) {
    stop("The specified 'age_in_months' column is not found in the dataframe.")
  }

  # Filter data based on mortality type
  if (type == "neonatal") {
    # Neonatal: deaths within the first month of life (age in months < 1)
    if (verbose) message("dm.chm: Filtering for neonatal age (age_in_months < 1)...")
    data_filtered <- data[data[[age_in_months]] < 1, ]
  } else if (type == "infant") {
    # Infant: Age in years should be < 1
    if (verbose) message("dm.chm: Filtering for infant age (age < 1)...")
    data_filtered <- data[data[[age_col]] < 1, ]
  } else if (type == "child") {
    # Child: Age in years should be between 1 and 4
    if (verbose) message("dm.chm: Filtering for child age (1 <= age < 5)...")
    data_filtered <- data[data[[age_col]] >= 1 & data[[age_col]] < 5, ]
  } else if (type == "under5") {
    # Under-5: Age in years should be < 5
    if (verbose) message("dm.chm: Filtering for under-5 age (age < 5)...")
    data_filtered <- data[data[[age_col]] < 5, ]
  }

  # Validate filtered data
  if (nrow(data_filtered) == 0) {
    stop("No data available for the specified mortality type.")
  }

  if (verbose) {
    message("dm.chm: Found ", nrow(data_filtered), " matching records in filtered data.")
  }

  # Calculate the mortality rate
  total_deaths <- sum(data_filtered[[deaths_col]])
  total_population <- sum(data_filtered[[live_births]])

  if (verbose) {
    message("dm.chm: Total deaths in filtered data: ", total_deaths)
    message("dm.chm: Total population/births in filtered data: ", total_population)
  }

  if (total_population == 0) {
    stop("Total population for the specified mortality type is zero. Calculation cannot proceed.")
  }

  mortality_rate <- (total_deaths / total_population) * 1000

  if (verbose) {
    message("dm.chm: Calculated mortality rate per 1,000: ", round(mortality_rate, 4))
  }

  class(mortality_rate) <- c("dem_chm", "numeric")
  attr(mortality_rate, "type") <- type
  return(mortality_rate)
}

#' @export
print.dem_chm <- function(x, ...) {
  type <- attr(x, "type")
  cat(sprintf("The %s mortality rate is: %.2f per 1,000\n", type, as.numeric(x)))
  invisible(x)
}
