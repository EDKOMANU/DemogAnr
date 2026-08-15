#' Calculate Maternal Mortality Metrics
#'
#' This function calculates the maternal mortality rate (maternal deaths per
#' 100,000 women of reproductive age) and the maternal mortality ratio
#' (maternal deaths per 100,000 live births).
#'
#' @param data A dataframe containing demographic data.
#' @param deaths_col The column name for maternal deaths.
#' @param births_col The column name for live births.
#' @param pop_women The column containing the number of women of reproductive age.
#' @param verbose Logical. If TRUE, prints progress messages during execution.
#'
#' @return An object of class `dem_mmr`: a list with `results` (the maternal
#'   mortality rate `MMR` per 100,000 women, and the maternal mortality ratio
#'   `MMR_Ratio` per 100,000 live births) and `modified_data` (the input data
#'   with per-row `mm_rate` and `mm_ratio` columns).
#' @examples
#' # Women of reproductive age in Ghana, 2021, by five-year age group, with an
#' # illustrative distribution of births and maternal deaths across those ages.
#' # The maternal mortality ratio is expressed per 100,000 live births.
#' data(ghpop2021)
#' women <- subset(ghpop2021, Sex == "Female" & Age >= 15 & Age <= 49)
#' women$grp <- paste0((women$Age %/% 5) * 5, "-", (women$Age %/% 5) * 5 + 4)
#' repro <- aggregate(Population ~ grp, data = women, FUN = sum)
#' repro$live_births     <- round(repro$Population * c(0.06, 0.16, 0.17, 0.14,
#'                                                     0.09, 0.04, 0.01))
#' repro$maternal_deaths <- round(repro$live_births * c(3, 2, 2, 3, 4, 6, 8) / 1e4)
#'
#' dem.mmr(repro, deaths_col = "maternal_deaths",
#'         births_col = "live_births", pop_women = "Population")
#' @references
#' World Health Organization (1992). \emph{International Statistical Classification of Diseases and Related Health Problems, 10th Revision (ICD-10)}. Geneva: World Health Organization. (Definitions of maternal death and maternal mortality measures.)
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. ISBN 978-0631226161.
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Chapter on health demography and maternal mortality.)
#'
#' @export
dem.mmr <- function(data, deaths_col, births_col, pop_women, verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) stop("Input 'data' must be a dataframe.")
  if (!all(c(deaths_col, births_col, pop_women) %in% colnames(data))) stop("Specified columns not found in the dataframe.")

  if (verbose) {
    message("dem.mmr: Starting maternal mortality calculations...")
    message("dem.mmr: Using columns - Deaths: '", deaths_col, "', Births: '", births_col, "', Women: '", pop_women, "'")
  }

  # Calculate Maternal Mortality Rate (MMR)
  if (verbose) message("dem.mmr: Calculating Maternal Mortality Rate (MMR)...")
  data$mm_rate <- (data[[deaths_col]] / data[[pop_women]]) * 100000
  mm_rate <- (sum(data[[deaths_col]]) / sum(data[[pop_women]])) * 100000
  if (verbose) message("dem.mmr: Overall MMR Rate per 100,000 women: ", round(mm_rate, 4))

  # Calculate Maternal Mortality Ratio (MMR)
  if (verbose) message("dem.mmr: Calculating Maternal Mortality Ratio (MMR)...")
  data$mm_ratio <- (data[[deaths_col]] / data[[births_col]]) * 100000
  mm_ratio <- sum(data[[deaths_col]]) / sum(data[[births_col]]) * 100000
  if (verbose) message("dem.mmr: Overall MMR Ratio per 100,000 live births: ", round(mm_ratio, 4))

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

