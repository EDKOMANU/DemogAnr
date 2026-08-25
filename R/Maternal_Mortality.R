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
#' @section Naming:
#' Almost everywhere in demography and public health, "MMR" means the maternal
#' mortality **ratio**: maternal deaths per 100,000 live births. In this
#' function's `results` list, however, the element `MMR` has always held the
#' maternal mortality **rate** (per 100,000 women of reproductive age), and
#' `MMR_Ratio` the ratio.
#'
#' Rather than change what `MMR` returns and silently alter existing results,
#' two unambiguous elements were added in version 0.3.0: `MMRate` for the rate
#' and `MMRatio` for the ratio. Prefer those. `MMR` and `MMR_Ratio` keep their
#' original values and will be removed in a future release.
#'
#' @return An object of class `dem_mmr`: a list with `results` -- `MMRate`
#'   (maternal deaths per 100,000 women of reproductive age), `MMRatio`
#'   (maternal deaths per 100,000 live births), and the deprecated aliases
#'   `MMR` (equal to `MMRate`) and `MMR_Ratio` (equal to `MMRatio`) -- and
#'   `modified_data`, the input with per-row `mm_rate` and `mm_ratio` columns.
#'
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

  # MMRate / MMRatio are the unambiguous names. MMR and MMR_Ratio are kept at
  # their original values for code written against earlier versions: note that
  # `MMR` here is the *rate*, which is the opposite of the usual convention,
  # which is why the clearer pair was added. See the Details section.
  results <- list(MMRate = mm_rate, MMRatio = mm_ratio,
                  MMR = mm_rate, MMR_Ratio = mm_ratio)

  out <- list(results = results, modified_data = data)
  class(out) <- "dem_mmr"
  return(out)
}

#' @export
print.dem_mmr <- function(x, ...) {
  cat(sprintf("Maternal mortality rate  (MMRate):  %.2f per 100,000 women\n",
              x$results$MMRate))
  cat(sprintf("Maternal mortality ratio (MMRatio): %.2f per 100,000 live births\n",
              x$results$MMRatio))
  n_show <- min(10L, nrow(x$modified_data))
  cat("\nBy row:\n")
  print(x$modified_data[seq_len(n_show), , drop = FALSE], row.names = FALSE)
  if (nrow(x$modified_data) > n_show) {
    cat(sprintf("... (%d rows)\n", nrow(x$modified_data)))
  }
  invisible(x)
}

