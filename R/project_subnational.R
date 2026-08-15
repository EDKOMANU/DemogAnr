#' Hierarchical Subnational Population Projection
#'
#' Projects lower-level geographies (e.g., Districts) while constraining them
#' to higher-level control totals (e.g., Regions or National) provided within the dataset.
#'
#' @param data Dataframe containing the subnational population data.
#' @param sub_group_col Character. Column name for the lower-level geography (e.g., "District").
#' @param parent_group_col Character. Optional. Column name for the higher-level geography (e.g., "Region").
#'        If NULL, assumes a single national parent group.
#' @param pop_cols Character vector of length 2. Column names for Base Year and Launch Year populations.
#' @param time_vals Numeric vector of length 2. The years corresponding to pop_cols.
#' @param target_times Numeric vector. The future years to project.
#' @param method Character. The projection method.
#' @param parent_totals_cols Character vector. Optional. Column names in `data` containing the pre-calculated
#'        projected totals for the parent group at `target_times`. Must match the length of `target_times`.
#' @param verbose Logical. If `TRUE`, prints detailed status messages to the console during projection.
#'
#' @return A dataframe appending the projected columns to the original data.
#' @examples
#' # Ghana's regions, projected from their 2000 census populations. A second
#' # observation is needed to establish a trend; here each region is grown at
#' # its own implied rate over the 2000-2010 decade.
#' data(region2000)
#' regions <- data.frame(
#'   Country  = region2000$Country,
#'   Region   = region2000$Region,
#'   pop2000  = region2000$base_pop,
#'   pop2010  = round(region2000$base_pop *
#'                    exp(10 * (region2000$cbr - region2000$cdr + region2000$nmr)))
#' )
#'
#' # Exponential extrapolation of each region to 2025 and 2030
#' math_project(regions,
#'   sub_group_col = "Region", parent_group_col = "Country",
#'   pop_cols = c("pop2000", "pop2010"), time_vals = c(2000, 2010),
#'   target_times = c(2025, 2030), method = "exponential"
#' )
#'
#' # Share-of-growth, constraining the regions to a national control total
#' regions$total2025 <- 40000000
#' math_project(regions,
#'   sub_group_col = "Region", parent_group_col = "Country",
#'   pop_cols = c("pop2000", "pop2010"), time_vals = c(2000, 2010),
#'   target_times = 2025, method = "share_of_growth",
#'   parent_totals_cols = "total2025"
#' )
#'
#' @references
#' Smith, S. K., Tayman, J., & Swanson, D. A. (2013). \emph{A Practitioner's Guide to State and Local Population Projections}. Dordrecht: Springer. \doi{10.1007/978-94-007-7551-0} (Trend extrapolation, constant-share, shift-share, and share-of-growth methods.)
#'
#' Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods and Materials of Demography} (2nd ed.). San Diego: Elsevier Academic Press. ISBN 978-0126419559. (Chapter on population projections.)
#'
#' Rowland, D. T. (2003). \emph{Demographic Methods and Concepts}. Oxford: Oxford University Press. ISBN 978-0198752639. (Chapter on population projections.)
#'
#' @export
math_project <- function(data, sub_group_col, parent_group_col = NULL,
                                 pop_cols, time_vals, target_times,
                                 method = c("linear", "exponential", "geometric",
                                            "constant_share", "shift_share", "share_of_growth"),
                                 parent_totals_cols = NULL, verbose = FALSE) {

  method <- match.arg(method)
  is_share_method <- method %in% c("constant_share", "shift_share", "share_of_growth")

  if (verbose) {
    message(sprintf("Starting hierarchical subnational population projection using method: '%s'", method))
  }

  # 1. Input Validation
  if (length(pop_cols) != 2 || length(time_vals) != 2) stop("pop_cols and time_vals must be length 2.")
  if (!all(pop_cols %in% names(data))) stop("pop_cols not found in data.")
  if (!(sub_group_col %in% names(data))) stop("sub_group_col not found in data.")

  if (is_share_method && !is.null(parent_totals_cols)) {
    if (length(parent_totals_cols) != length(target_times)) {
      stop("Length of parent_totals_cols must match length of target_times.")
    }
    if (!all(parent_totals_cols %in% names(data))) stop("parent_totals_cols not found in data.")
  }

  # Set a dummy parent if none is provided (defaults to treating everything as one national group)
  if (is.null(parent_group_col)) {
    data$Dummy_Parent <- "National"
    parent_group_col <- "Dummy_Parent"
    added_dummy <- TRUE
  } else {
    if (!(parent_group_col %in% names(data))) stop("parent_group_col not found in data.")
    added_dummy <- FALSE
  }

  t0 <- time_vals[1]
  tn <- time_vals[2]
  n_years <- tn - t0

  # 2. Mathematical Core Functions
  calc_linear <- function(P0, Pn, t_diff) Pn + ((Pn - P0) / n_years) * t_diff
  calc_exponential <- function(P0, Pn, t_diff) Pn * exp((log(Pn / P0) / n_years) * t_diff)
  calc_geometric <- function(P0, Pn, t_diff) Pn * (1 + ((Pn / P0)^(1 / n_years) - 1))^t_diff

  # 3. Split-Apply-Combine Processing
  # We split the data by the parent group (e.g., isolate all districts in Ashanti, then all in Greater Accra)
  data_split <- split(data, data[[parent_group_col]])

  processed_list <- lapply(data_split, function(grp_data) {
    if (verbose) {
      grp_name <- unique(grp_data[[parent_group_col]])[1]
      message(sprintf("  Processing parent group: %s (%d sub-groups)", grp_name, nrow(grp_data)))
    }

    P0 <- grp_data[[pop_cols[1]]]
    Pn <- grp_data[[pop_cols[2]]]

    Parent0 <- sum(P0, na.rm = TRUE)
    Parentn <- sum(Pn, na.rm = TRUE)

    # Extract parent control totals for the target years
    if (is_share_method && !is.null(parent_totals_cols)) {
      # Assumes the parent total is duplicated across the district rows, so we just take row 1
      grp_target_totals <- as.numeric(grp_data[1, parent_totals_cols])
    } else if (is_share_method) {
      # If no totals provided, default to extrapolating the parent sum exponentially
      grp_target_totals <- sapply(target_times, function(t) calc_exponential(Parent0, Parentn, t - tn))
    }

    # Process each target year
    proj_cols <- lapply(seq_along(target_times), function(i) {
      t_target <- target_times[i]
      t_diff <- t_target - tn

      if (method == "linear") res <- calc_linear(P0, Pn, t_diff)
      if (method == "exponential") res <- calc_exponential(P0, Pn, t_diff)
      if (method == "geometric") res <- calc_geometric(P0, Pn, t_diff)

      if (method == "constant_share") {
        S_n <- Pn / Parentn
        res <- S_n * grp_target_totals[i]
      }
      if (method == "shift_share") {
        S_0 <- P0 / Parent0
        S_n <- Pn / Parentn
        shift <- (S_n - S_0) / n_years
        S_t <- S_n + (shift * t_diff)

        S_t <- ifelse(S_t < 0, 0, S_t) # Floor negative shares
        S_t <- S_t / sum(S_t, na.rm = TRUE) # Re-normalize
        res <- S_t * grp_target_totals[i]
      }
      if (method == "share_of_growth") {
        S_growth <- (Pn - P0) / (Parentn - Parent0)
        Parent_proj_growth <- grp_target_totals[i] - Parentn
        res <- Pn + (S_growth * Parent_proj_growth)
      }

      out_col <- data.frame(res)
      colnames(out_col) <- paste0("proj_", t_target)
      return(out_col)
    })

    # Bind the projected columns back to this specific parent group's data
    cbind(grp_data, do.call(cbind, proj_cols))
  })

  # 4. Recombine and Clean Up
  final_data <- do.call(rbind, unname(processed_list))
  rownames(final_data) <- NULL

  if (added_dummy) final_data$Dummy_Parent <- NULL

  return(final_data)
}
