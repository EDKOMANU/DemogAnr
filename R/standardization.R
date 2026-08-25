#' Age Standardization: Comparing Two Populations
#'
#' Removes the effect of differing age structures when comparing the crude
#' rates of two populations, following Preston, Heuveline and Guillot (2001,
#' Chapter 2). Each population's age-specific rates are re-weighted by a common
#' standard age distribution -- by default the average of the two populations'
#' age compositions, their recommended choice for a two-population comparison.
#'
#' @details
#' For populations 1 and 2 with age-specific rates \eqn{M^1_i, M^2_i} and age
#' compositions \eqn{C^1_i = N^1_i / \sum N^1_i} and \eqn{C^2_i}, the crude rate
#' of each is \eqn{CDR = \sum_i C_i M_i}. The age-standardized crude rate applies
#' a common standard composition \eqn{C^s_i}:
#' \deqn{ASCDR = \sum_i M_i \, C^s_i .}
#' By default \eqn{C^s_i = (C^1_i + C^2_i)/2}; set `standard` to `"first"` or
#' `"second"` to use one population's composition, or supply an external
#' standard with `std_pop_col`. The function also reports the comparative
#' mortality ratio \eqn{CMR = \sum_i N^1_i M^1_i / \sum_i N^1_i M^2_i}
#' (population 1's deaths relative to those expected under population 2's rates).
#'
#' With the average standard, the difference between the two age-standardized
#' rates equals the rate component of the Kitagawa decomposition
#' ([decompose_rates()]) of the same difference.
#'
#' @param data A data frame with one row per age group.
#' @param age_col Column name for age groups.
#' @param rate1_col,rate2_col Columns of age-specific rates for populations 1 and 2.
#' @param pop1_col,pop2_col Columns of population counts (or proportions) for
#'   populations 1 and 2.
#' @param standard How to choose the standard age distribution: `"average"`
#'   (default), `"first"`, or `"second"`. Ignored if `std_pop_col` is supplied.
#' @param std_pop_col Optional column of an external standard population by age;
#'   overrides `standard`.
#' @param labels Length-2 character vector naming the two populations.
#' @param per Scaling of the returned rates (default 1000).
#'
#' @return An object of class `dem_std`: a list with the crude rates (`crude`),
#'   the age-standardized rates (`standardized`), their `ratio`, the comparative
#'   mortality ratio (`CMR`), and a per-age `table`.
#'
#' @examples
#' # Preston et al. (2001) Box 2.1: Sweden vs Kazakhstan, females, 1992
#' box21 <- data.frame(
#'   age = c("0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
#'           "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79",
#'           "80-84","85+"),
#'   pop_Sw = c(0.0136,0.0524,0.0559,0.0548,0.0604,0.0655,0.0709,0.0641,0.0654,
#'              0.0703,0.0730,0.0552,0.0481,0.0493,0.0512,0.0508,0.0420,0.0321,
#'              0.0251),
#'   pop_K  = c(0.0200,0.0868,0.1011,0.0929,0.0828,0.0716,0.0843,0.0842,0.0704,
#'              0.0561,0.0327,0.0579,0.0347,0.0430,0.0295,0.0178,0.0172,0.0102,
#'              0.0068),
#'   m_Sw = c(0.00467,0.00008,0.00013,0.00014,0.00023,0.00030,0.00032,0.00050,
#'            0.00069,0.00117,0.00201,0.00305,0.00461,0.00759,0.01226,0.02026,
#'            0.03664,0.06815,0.15729),
#'   m_K  = c(0.02137,0.00162,0.00045,0.00037,0.00078,0.00108,0.00103,0.00132,
#'            0.00182,0.00288,0.00430,0.00571,0.01082,0.01392,0.02679,0.03998,
#'            0.05469,0.10159,0.18030)
#' )
#' standardize(box21, age_col = "age",
#'             rate1_col = "m_Sw", rate2_col = "m_K",
#'             pop1_col = "pop_Sw", pop2_col = "pop_K",
#'             labels = c("Sweden", "Kazakhstan"))
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapter 2; Box 2.1.)
#'
#' @export
standardize <- function(data, age_col, rate1_col, rate2_col,
                        pop1_col, pop2_col,
                        standard = c("average", "first", "second"),
                        std_pop_col = NULL, labels = c("1", "2"), per = 1000) {
  standard <- match.arg(standard)
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (nm in c(age_col, rate1_col, rate2_col, pop1_col, pop2_col)) {
    if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
  }
  age <- data[[age_col]]
  M1 <- as.numeric(data[[rate1_col]]); M2 <- as.numeric(data[[rate2_col]])
  N1 <- as.numeric(data[[pop1_col]]);  N2 <- as.numeric(data[[pop2_col]])
  c1 <- N1 / sum(N1); c2 <- N2 / sum(N2)

  cs <- if (!is.null(std_pop_col)) {
    Ns <- as.numeric(data[[std_pop_col]]); Ns / sum(Ns)
  } else if (standard == "first") {
    c1
  } else if (standard == "second") {
    c2
  } else {
    (c1 + c2) / 2
  }

  cdr1 <- sum(c1 * M1) * per; cdr2 <- sum(c2 * M2) * per
  ascdr1 <- sum(M1 * cs) * per; ascdr2 <- sum(M2 * cs) * per
  cmr <- sum(N1 * M1) / sum(N1 * M2)

  out <- list(
    labels = labels, per = per,
    standard = if (!is.null(std_pop_col)) "external" else standard,
    crude = stats::setNames(c(cdr1, cdr2), labels),
    standardized = stats::setNames(c(ascdr1, ascdr2), labels),
    ratio = ascdr1 / ascdr2,
    CMR = cmr,
    table = data.frame(age = age, C1 = c1, C2 = c2, Cstd = cs,
                       M1 = M1, M2 = M2,
                       M1Cstd = M1 * cs * per, M2Cstd = M2 * cs * per)
  )
  class(out) <- "dem_std"
  out
}

#' @export
print.dem_std <- function(x, ...) {
  cat(sprintf("Age standardization: %s (1) vs %s (2)\n",
              x$labels[1], x$labels[2]))
  cat(sprintf("Standard age distribution: %s\n\n", x$standard))
  tb <- x$table
  disp <- data.frame(
    age = as.character(tb$age),
    C1 = round(tb$C1, 4), C2 = round(tb$C2, 4), Cstd = round(tb$Cstd, 4),
    M1 = round(tb$M1, 5), M2 = round(tb$M2, 5),
    M1Cstd = round(tb$M1Cstd, 4), M2Cstd = round(tb$M2Cstd, 4),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  cat(sprintf("\n  Crude rate:            (1) %.2f   (2) %.2f   per %d\n",
              x$crude[1], x$crude[2], x$per))
  cat(sprintf("  Age-standardized rate: (1) %.2f   (2) %.2f   per %d\n",
              x$standardized[1], x$standardized[2], x$per))
  cat(sprintf("  Ratio of standardized rates (1/2): %.3f\n", x$ratio))
  cat(sprintf("  Comparative mortality ratio (CMR, 1 vs 2): %.3f\n", x$CMR))
  invisible(x)
}


#' Kitagawa Decomposition of a Difference Between Two Rates
#'
#' Decomposes the difference between two populations' crude rates into a
#' component due to differences in age-specific rates and a component due to
#' differences in age composition (Kitagawa 1955; Preston et al. 2001, Section
#' 2.3).
#'
#' @details
#' With age composition \eqn{c_a = N_a / \sum N_a} and age-specific rates
#' \eqn{M_a}, the difference between the crude rates of populations 1 and 2 is
#' split as
#' \deqn{CDR_1 - CDR_2 = \sum_a (M^1_a - M^2_a)\frac{c^1_a + c^2_a}{2}
#'   + \sum_a (c^1_a - c^2_a)\frac{M^1_a + M^2_a}{2},}
#' the first term being the rate (or effect) component and the second the age
#' composition component.
#'
#' @param data A data frame with one row per age group.
#' @param age_col Column name for age groups.
#' @param rate1_col,rate2_col Columns of age-specific rates for populations 1 and 2.
#' @param pop1_col,pop2_col Columns of population counts for populations 1 and 2.
#' @param per Scaling of the returned components (default 1000).
#'
#' @return An object of class `dem_decomp`: a list with the `total` difference,
#'   the `rate_component`, the `composition_component`, and a per-age `table`.
#'
#' @examples
#' d <- data.frame(
#'   age   = c("0-14","15-44","45-64","65+"),
#'   m1    = c(0.002, 0.001, 0.008, 0.060),
#'   m2    = c(0.003, 0.0015, 0.010, 0.065),
#'   n1    = c(300000, 400000, 200000, 100000),
#'   n2    = c(200000, 350000, 250000, 200000)
#' )
#' decompose_rates(d, "age", "m1", "m2", "n1", "n2")
#'
#' @references
#' Kitagawa, E. M. (1955). Components of a difference between two rates. \emph{Journal of the American Statistical Association}, 50(272), 1168-1194. \doi{10.1080/01621459.1955.10501299}
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography}. Oxford: Blackwell Publishers. (Section 2.3.)
#'
#' @export
decompose_rates <- function(data, age_col, rate1_col, rate2_col,
                            pop1_col, pop2_col, per = 1000) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (nm in c(age_col, rate1_col, rate2_col, pop1_col, pop2_col)) {
    if (!(nm %in% names(data))) stop("Column not found in 'data': ", nm)
  }
  age <- data[[age_col]]
  M1 <- as.numeric(data[[rate1_col]]); M2 <- as.numeric(data[[rate2_col]])
  N1 <- as.numeric(data[[pop1_col]]);  N2 <- as.numeric(data[[pop2_col]])
  c1 <- N1 / sum(N1); c2 <- N2 / sum(N2)

  rate_contrib <- (M1 - M2) * (c1 + c2) / 2 * per
  comp_contrib <- (c1 - c2) * (M1 + M2) / 2 * per

  cdr1 <- sum(c1 * M1) * per
  cdr2 <- sum(c2 * M2) * per

  out <- list(
    total = cdr1 - cdr2,
    rate_component = sum(rate_contrib),
    composition_component = sum(comp_contrib),
    crude_rate1 = cdr1,
    crude_rate2 = cdr2,
    per = per,
    table = data.frame(age = age, rate1 = M1, rate2 = M2,
                       share1 = c1, share2 = c2,
                       rate_effect = rate_contrib,
                       composition_effect = comp_contrib)
  )
  class(out) <- "dem_decomp"
  out
}

#' @export
print.dem_decomp <- function(x, ...) {
  cat("Kitagawa decomposition of a difference between two rates\n\n")
  tb <- x$table
  disp <- data.frame(
    age = as.character(tb$age),
    rate1 = round(tb$rate1, 5), rate2 = round(tb$rate2, 5),
    share1 = round(tb$share1, 4), share2 = round(tb$share2, 4),
    rate_effect = round(tb$rate_effect, 4),
    comp_effect = round(tb$composition_effect, 4),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  cat(sprintf("\n  Crude rate 1: %.4f per %d\n", x$crude_rate1, x$per))
  cat(sprintf("  Crude rate 2: %.4f per %d\n", x$crude_rate2, x$per))
  cat(sprintf("  Total difference (1 - 2): %.4f\n", x$total))
  cat(sprintf("    Rate component (sum):        %.4f\n", x$rate_component))
  cat(sprintf("    Composition component (sum): %.4f\n", x$composition_component))
  invisible(x)
}


#' Arriaga Decomposition of a Difference in Life Expectancy
#'
#' Decomposes the difference in life expectancy at birth between two life tables
#' into the contributions of mortality differences in each age group, using
#' Arriaga's (1984) discrete method (Preston et al. 2001, equations 3.11-3.12).
#'
#' @details
#' For each closed age interval the contribution is
#' \deqn{{}_n\Delta_x = \frac{l^1_x}{l_0}\left(\frac{{}_nL^2_x}{l^2_x}
#'   - \frac{{}_nL^1_x}{l^1_x}\right)
#'   + \frac{T^2_{x+n}}{l_0}\left(\frac{l^1_x}{l^2_x}
#'   - \frac{l^1_{x+n}}{l^2_{x+n}}\right),}
#' the first term being the direct effect and the second the indirect plus
#' interaction effect. For the open-ended interval only the direct effect
#' applies. The contributions sum to \eqn{e_0^2 - e_0^1}. Superscripts 1 and 2
#' refer to `lt1` (baseline) and `lt2` (comparison).
#'
#' @param lt1,lt2 Data frames for the two life tables, each with columns for
#'   age and the life-table functions `lx`, `Lx` and `Tx` (as returned in the
#'   `lifetable` component of [lifetable()]). Both must share the same age
#'   groups and radix.
#' @param age_col,lx_col,Lx_col,Tx_col Column names in `lt1`/`lt2`.
#' @param graph Logical. If `TRUE` (default), a \pkg{ggplot2} bar chart of the
#'   age-specific contributions is attached to the result as `$plot`;
#'   retrieve it with `plot()`.
#'
#' @return An object of class `dem_le_decomp`: a list with the `total`
#'   difference in life expectancy at birth, its decomposition into `direct` and
#'   `indirect` effects, and a per-age `table` of contributions.
#'
#' @examples
#' # Male-female gap in life expectancy, Ghana 2021
#' data(ghmort2021)
#' male   <- subset(ghmort2021, Sex == "Male")
#' female <- subset(ghmort2021, Sex == "Female")
#' m <- lifetable(male,   age = "Age", pop = "Population", Dx = "Deaths",
#'                sex = "male",   nax_method = "keyfitz")
#' f <- lifetable(female, age = "Age", pop = "Population", Dx = "Deaths",
#'                sex = "female", nax_method = "keyfitz")
#' decompose_LE(m$lifetable, f$lifetable)
#'
#' @references
#' Arriaga, E. E. (1984). Measuring and explaining the change in life expectancies. \emph{Demography}, 21(1), 83-96. \doi{10.2307/2061029}
#'
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography}. Oxford: Blackwell Publishers. (Equations 3.11-3.12.)
#'
#' @export
decompose_LE <- function(lt1, lt2, age_col = "Age",
                         lx_col = "lx", Lx_col = "Lx", Tx_col = "Tx",
                         graph = TRUE) {
  for (df in list(lt1, lt2)) {
    for (nm in c(age_col, lx_col, Lx_col, Tx_col)) {
      if (!(nm %in% names(df))) stop("Column not found in a life table: ", nm)
    }
  }
  age <- lt1[[age_col]]
  if (!identical(age, lt2[[age_col]])) stop("The two life tables must share the same age groups.")
  k <- length(age)

  l1 <- as.numeric(lt1[[lx_col]]); l2 <- as.numeric(lt2[[lx_col]])
  L1 <- as.numeric(lt1[[Lx_col]]); L2 <- as.numeric(lt2[[Lx_col]])
  T1 <- as.numeric(lt1[[Tx_col]]); T2 <- as.numeric(lt2[[Tx_col]])
  # Arriaga's contributions are expressed per birth of a common radix. With
  # mismatched radices every term is scaled wrongly and the "difference" is
  # dominated by the radix ratio, so refuse rather than return a plausible
  # but meaningless number.
  if (!isTRUE(all.equal(l1[1], l2[1]))) {
    stop("The two life tables must share the same radix (l0): ",
         format(l1[1]), " in 'lt1' but ", format(l2[1]), " in 'lt2'. ",
         "Rebuild one of them with the other's 'radix'.")
  }
  l0 <- l1[1]

  direct <- numeric(k); indirect <- numeric(k)
  for (i in 1:(k - 1)) {
    direct[i] <- (l1[i] / l0) * (L2[i] / l2[i] - L1[i] / l1[i])
    indirect[i] <- (T2[i + 1] / l0) * (l1[i] / l2[i] - l1[i + 1] / l2[i + 1])
  }
  # Open-ended interval: direct effect only
  direct[k] <- (l1[k] / l0) * (T2[k] / l2[k] - T1[k] / l1[k])
  indirect[k] <- 0

  contribution <- direct + indirect
  e0_1 <- T1[1] / l0; e0_2 <- T2[1] / l0

  out <- list(
    total = e0_2 - e0_1,
    e0_1 = e0_1, e0_2 = e0_2,
    direct = sum(direct),
    indirect = sum(indirect),
    table = data.frame(age = age, direct = direct, indirect = indirect,
                       contribution = contribution)
  )
  if (graph) out$plot <- .plot_le_decomp(out$table)
  class(out) <- "dem_le_decomp"
  out
}

#' @export
print.dem_le_decomp <- function(x, ...) {
  cat("Arriaga decomposition of a difference in life expectancy at birth\n")
  cat(sprintf("  e0 (table 1): %.3f\n", x$e0_1))
  cat(sprintf("  e0 (table 2): %.3f\n", x$e0_2))
  cat(sprintf("  Total difference (2 - 1): %.3f years\n", x$total))
  cat(sprintf("    Direct effects:              %.3f\n", x$direct))
  cat(sprintf("    Indirect + interaction:      %.3f\n", x$indirect))
  cat("\nContribution by age group:\n")
  tb <- x$table
  tb[-1] <- lapply(tb[-1], function(v) round(v, 4))
  print(tb, row.names = FALSE)
  invisible(x)
}

#' @export
plot.dem_le_decomp <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot available; call decompose_LE(..., graph = TRUE).")
  x$plot
}
