#' Calculate Child Mortality Metrics by Age
#'
#' This function calculates child mortality metrics (neonatal, infant, child, under-5)
#' based on single-age data. It requires columns for age, population, and deaths.
#'
#' @param data A dataframe containing demographic data.
#' @param age_col The column name representing the age (in years or months for neonates).
#' @param live_births The column name representing the population at each age.
#' @param deaths_col The column name representing deaths at each age.
#' @param type The type of mortality calculation to perform:
#'             "neonatal", "infant", "child", or "under5".
#' @param age_in_months (Optional) The column name for age in months for neonatal mortality.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#' @return A numeric value representing the mortality rate per 1,000 live births.
#' @examples
#' demo_data <- data.frame(
#'   age = 0:10,
#'   population = c(1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500),
#'   deaths = c(30, 20, 15, 10, 5, 2, 1, 1, 0, 0, 0)
#' )
#' dm.chm(demo_data,age_col = "age",
#'                                  live_births = "population",
#'                                  deaths_col = "deaths",
#'                                  type = "under5")
#' @references
#' Pollard, A. H., Yusuf, F., & Pollard, G. N. (1990). \emph{Demographic Techniques} (3rd ed.). Pergamon Press.
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography}. (Chapter 1, distinguishing rates vs. probabilities).
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
    # Neonatal: Age in months should be <= 1 month
    if (verbose) message("dm.chm: Filtering for neonatal age (age_in_months <= 1)...")
    data_filtered <- data[data[[age_in_months]] <= 1, ]
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
  cat(sprintf("The %s mortality rate is: %.2f per 1,000 live births\n", type, as.numeric(x)))
  invisible(x)
}
#consustency in the naming of the functions
