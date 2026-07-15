#' Split Age Groups into Single-Year Ages
#'
#' This function takes grouped population data (e.g., 5-year age groups) and
#' distributes the population into single-year age groups using the Karup-King
#' third-difference osculatory interpolation coefficients. The user can specify
#' the column containing age groups, one or more population columns, and
#' (optionally) custom coefficient matrices for the first, middle, and last
#' groups.
#'
#' Each coefficient matrix has 5 rows (one per single year within the 5-year
#' group) and 3 columns (weights applied to the neighbouring 5-year groups).
#' The defaults are the standard Karup-King multipliers (Siegel & Swanson 2004,
#' Appendix C): [first_coef] for the first group (columns: current, next,
#' next+1), [middle_coef] for interior groups (columns: previous, current,
#' next), and [last_coef] for the last group (columns: previous-1, previous,
#' current).
#'
#' @param df A data frame containing grouped population data.
#' @param age_col A string specifying the column in `df` that contains age group
#'   labels in the form "start-end" (e.g., "15-19"). Open-ended groups (e.g.
#'   "85+") are not supported and should be excluded or closed before calling.
#' @param pops A string or vector of strings specifying one or more columns in `df` that contain population data.
#' @param first_coef A 5x3 matrix of coefficients for splitting the first age
#'   group. Defaults to the packaged Karup-King multipliers ([first_coef]).
#' @param middle_coef A 5x3 matrix of coefficients for splitting middle age
#'   groups. Defaults to the packaged Karup-King multipliers ([middle_coef]).
#' @param last_coef A 5x3 matrix of coefficients for splitting the last age
#'   group. Defaults to the packaged Karup-King multipliers ([last_coef]).
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during interpolation.
#'
#' @return A data frame with single-year ages and corresponding population
#'   estimates for each specified population column. Ages beyond the range of
#'   the input age groups (overflow ages) are removed.
#'
#' @examples
#' data(data)
#' head(data)
#'
#' # Split the 5-year age groups into single years for the 2021 population,
#' # using the packaged Karup-King coefficients:
#' single <- karup_king(df = data, age_col = "age_col", pops = "2021")
#' head(single)
#'
#' # Population totals are (approximately) preserved:
#' sum(data[["2021"]])
#' sum(single[["2021"]])
#'
#' @references
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Appendix C: Karup-King osculatory interpolation multipliers.)
#'
#' Shryock, H. S., Siegel, J. S., & Larmon, E. A. (1973). \emph{The Methods and Materials of Demography}. Washington, DC: US Bureau of the Census. (Appendix C: Interpolation and graduation.)
#'
#' @export
karup_king <- function(df, age_col = "age_group", pops = "population",
                       first_coef = NULL, middle_coef = NULL, last_coef = NULL,
                       verbose = FALSE) {
  # Default to the packaged Karup-King coefficient matrices
  if (is.null(first_coef)) first_coef <- as.matrix(DemogAnr::first_coef)
  if (is.null(middle_coef)) middle_coef <- as.matrix(DemogAnr::middle_coef)
  if (is.null(last_coef)) last_coef <- as.matrix(DemogAnr::last_coef)
  first_coef <- as.matrix(first_coef)
  middle_coef <- as.matrix(middle_coef)
  last_coef <- as.matrix(last_coef)

  if (!is.data.frame(df)) stop("Input 'df' must be a dataframe.")
  if (!(age_col %in% names(df))) stop("age_col not found in df.")
  if (!all(pops %in% names(df))) stop("Some pops columns not found in df.")
  for (m in list(first_coef, middle_coef, last_coef)) {
    if (!all(dim(m) == c(5, 3))) stop("Coefficient matrices must be 5x3.")
  }
  if (any(!grepl("^\\s*\\d+\\s*-\\s*\\d+\\s*$", df[[age_col]]))) {
    stop("All values in age_col must be closed ranges of the form 'start-end' (e.g. '15-19'). Open-ended groups such as '85+' are not supported.")
  }
  if (nrow(df) < 3) stop("At least 3 age groups are required for Karup-King interpolation.")

  if (verbose) {
    message(sprintf("Running Karup-King interpolation on %d age groups...", nrow(df)))
  }

  # Split the specified age column into start and end ages
  df <- df |>
    tidyr::separate(col = {{age_col}}, into = c("start_age", "end_age"),
                    sep = "-", convert = TRUE) |>
    dplyr::arrange(start_age)

  single_ages <- list()

  for (i in 1:nrow(df)) {
    current_group <- df[i, ]

    # Determine coefficients and neighbors based on group position
    if (i == 1) {
      coef_matrix <- first_coef
      pop1 <- current_group[, pops, drop = FALSE]
      pop2 <- if (i + 1 <= nrow(df)) df[i + 1, pops, drop = FALSE] else 0 * pop1
      pop3 <- if (i + 2 <= nrow(df)) df[i + 2, pops, drop = FALSE] else 0 * pop1
    } else if (i == nrow(df)) {
      coef_matrix <- last_coef
      pop1 <- if (i - 2 >= 1) df[i - 2, pops, drop = FALSE] else 0 * current_group[, pops, drop = FALSE]
      pop2 <- if (i - 1 >= 1) df[i - 1, pops, drop = FALSE] else 0 * current_group[, pops, drop = FALSE]
      pop3 <- current_group[, pops, drop = FALSE]
    } else {
      coef_matrix <- middle_coef
      pop1 <- df[i - 1, pops, drop = FALSE]
      pop2 <- current_group[, pops, drop = FALSE]
      pop3 <- df[i + 1, pops, drop = FALSE]
    }

    # Calculate single-year populations for each population column
    for (year_idx in 1:5) {
      age <- current_group$start_age + (year_idx - 1)

      # Calculate population values using matrix operations
      pop_values <- round(coef_matrix[year_idx, 1] * pop1 +
        coef_matrix[year_idx, 2] * pop2 +
        coef_matrix[year_idx, 3] * pop3, 0)

      single_ages[[length(single_ages) + 1]] <- dplyr::bind_cols(
        data.frame(age = age),
        pop_values
      )
    }
  }

  # Combine results and remove overflow ages
  result <- dplyr::bind_rows(single_ages) |>
    dplyr::filter(age <= max(df$end_age)) |>
    dplyr::mutate(age = as.character(age))

  return(result)
}
