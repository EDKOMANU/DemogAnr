#' Calculate Maternal Mortality Metrics
#'
#' This function calculates the Maternal Mortality Rate (MMR) and Maternal Mortality Ratio (MMR).
#'
#' @param data A dataframe containing demographic data.
#' @param deaths_col The column name for maternal deaths.
#' @param births_col The column name for live births.
#' @param pop_women The column containing the number of women in the reporductive age
#'
#' @return A list with Maternal Mortality Rate (MMR) and Maternal Mortality Ratio (MMR).
#' @examples
#' demo_data <- data.frame(
#'   age=c("15-24", "25-34", "35-44"),
#'   maternal_deaths = c(5, 10, 15),
#'   live_births = c(5000, 7000, 6000),
#'   women =c(20000, 22000, 32000)
#' )
#' dem.mmr(demo_data, deaths_col = "maternal_deaths", births_col = "live_births", pop_women = "women")
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapters 2 and 5)
#'
#' Newell, C. (1988). \emph{Methods and Models in Demography}. New York: Guilford Press. (Chapters 4 and 6)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). Emerald Group Publishing. (Chapters 14 for Maternal Mortality)
#'
#' @export
dem.mmr <- function(data, deaths_col, births_col, pop_women) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (!all(c(deaths_col, births_col, pop_women) %in% colnames(data))) stop("Specified columns not found in the dataframe.")

  # Calculate Maternal Mortality Rate (MMR)
  data$mm_rate <- (data[[deaths_col]] / data[[pop_women]]) * 100000
  mm_rate <- (sum(data[[deaths_col]]) / sum(data[[pop_women]])) * 100000
  # Calculate Maternal Mortality Ratio (MMR)
  data$mm_ratio <- (data[[deaths_col]] / data[[births_col]]) * 100000
  mm_ratio <- sum(data[[deaths_col]]) / sum(data[[births_col]]) * 100000
  results <- list(MMR = mm_rate, MMR_Ratio = mm_ratio)

  out <- list(results = results, modified_data = data)
  class(out) <- "dem_mmr"
  return(out)
}

#' @export
print.dem_mmr <- function(x, ...) {
  cat("Maternal Mortality Rate (MMR)  per 100,000 women:\n")
  print(x$results$MMR)
  cat("\nMaternal Mortality Ratio (MMR)  per 100,000 live births:\n")
  print(x$results$MMR_Ratio)
  cat("\nModified Data:\n")
  print(x$modified_data)
  invisible(x)
}

