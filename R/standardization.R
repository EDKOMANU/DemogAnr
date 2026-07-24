#' Age Standardization of Rates (direct and indirect)
#'
#' Computes age-standardized rates to make crude rates comparable across
#' populations with different age structures, using either the direct or the
#' indirect method (Preston, Heuveline and Guillot 2001, Chapter 2).
#'
#' @details
#' \strong{Direct standardization} applies the study population's age-specific
#' rates to a standard age distribution:
#' \deqn{ASR = \sum_a w_a M_a, \qquad w_a = N^{std}_a / \sum_a N^{std}_a,}
#' where \eqn{M_a} are the study rates and \eqn{N^{std}_a} the standard
#' population. It requires `rate_col` and `std_pop_col`.
#'
#' \strong{Indirect standardization} applies a standard schedule of rates to the
#' study population's age distribution. It yields the standardized mortality
#' ratio
#' \deqn{SMR = D / \sum_a N_a M^{std}_a}
#' (observed deaths over expected deaths) and, when a standard population is
#' also supplied, an indirectly standardized rate \eqn{SMR \times CDR^{std}}. It
#' requires `pop_col`, `deaths_col` and `std_rate_col`.
#'
#' @param data A data frame with one row per age group.
#' @param age_col Column name for age groups.
#' @param method `"direct"` (default) or `"indirect"`.
#' @param rate_col Column of study age-specific rates (direct method).
#' @param pop_col Column of study population by age (indirect method).
#' @param deaths_col Column of study deaths by age (indirect method).
#' @param std_pop_col Column of standard population by age (direct method, and
#'   optional for indirect to obtain the indirectly standardized rate).
#' @param std_rate_col Column of standard age-specific rates (indirect method).
#' @param per Scaling of the returned rate (default 1000).
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @return An object of class `dem_std`: a list with the standardized rate (and,
#'   for the indirect method, the `SMR`), the crude rate, and a per-age
#'   `table`.
#'
#' @examples
#' pop <- data.frame(
#'   age      = c("0-14","15-44","45-64","65+"),
#'   rate     = c(0.002, 0.001, 0.008, 0.060),   # study age-specific rates
#'   std_pop  = c(250000, 400000, 250000, 100000)
#' )
#' standardize(pop, age_col = "age", method = "direct",
#'             rate_col = "rate", std_pop_col = "std_pop")
#'
#' @references
#' Preston, S. H., Heuveline, P., & Guillot, M. (2001). \emph{Demography: Measuring and Modeling Population Processes}. Oxford: Blackwell Publishers. (Chapter 2: Age-standardization.)
#'
#' @export
standardize <- function(data, age_col, method = c("direct", "indirect"),
                        rate_col = NULL, pop_col = NULL, deaths_col = NULL,
                        std_pop_col = NULL, std_rate_col = NULL,
                        per = 1000, verbose = FALSE) {
  method <- match.arg(method)
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!(age_col %in% names(data))) stop("age_col not found in 'data'.")
  age <- data[[age_col]]

  if (method == "direct") {
    if (is.null(rate_col) || is.null(std_pop_col)) {
      stop("Direct standardization requires rate_col and std_pop_col.")
    }
    Mx <- as.numeric(data[[rate_col]])
    Ns <- as.numeric(data[[std_pop_col]])
    w <- Ns / sum(Ns)
    asr <- sum(w * Mx) * per
    crude_std <- sum(w * Mx) * per # crude rate of the standard under study rates
    if (verbose) message(sprintf("Direct age-standardized rate: %.4f per %d", asr, per))
    tbl <- data.frame(age = age, rate = Mx, std_pop = Ns, weight = w,
                      contribution = w * Mx * per)
    out <- list(method = "direct", standardized_rate = asr,
                crude_rate = asr, per = per, table = tbl)
  } else {
    if (is.null(pop_col) || is.null(deaths_col) || is.null(std_rate_col)) {
      stop("Indirect standardization requires pop_col, deaths_col and std_rate_col.")
    }
    Nx <- as.numeric(data[[pop_col]])
    Dx <- as.numeric(data[[deaths_col]])
    Ms <- as.numeric(data[[std_rate_col]])
    observed <- sum(Dx)
    expected <- sum(Nx * Ms)
    smr <- observed / expected
    crude_obs <- observed / sum(Nx) * per
    is_rate <- NA_real_
    if (!is.null(std_pop_col)) {
      Ns <- as.numeric(data[[std_pop_col]])
      cdr_std <- sum(Ns * Ms) / sum(Ns) * per
      is_rate <- smr * cdr_std
    }
    if (verbose) message(sprintf("SMR = %.4f (observed %.0f / expected %.1f)",
                                 smr, observed, expected))
    tbl <- data.frame(age = age, pop = Nx, deaths = Dx, std_rate = Ms,
                      expected_deaths = Nx * Ms)
    out <- list(method = "indirect", SMR = smr,
                standardized_rate = is_rate, crude_rate = crude_obs,
                observed_deaths = observed, expected_deaths = expected,
                per = per, table = tbl)
  }
  class(out) <- "dem_std"
  out
}

#' @export
print.dem_std <- function(x, ...) {
  if (x$method == "direct") {
    cat(sprintf("Direct age-standardized rate: %.4f per %d\n",
                x$standardized_rate, x$per))
  } else {
    cat(sprintf("Indirect standardization\n  SMR: %.4f (observed %.0f / expected %.1f)\n",
                x$SMR, x$observed_deaths, x$expected_deaths))
    if (!is.na(x$standardized_rate)) {
      cat(sprintf("  Indirectly standardized rate: %.4f per %d\n",
                  x$standardized_rate, x$per))
    }
    cat(sprintf("  Crude rate (observed): %.4f per %d\n", x$crude_rate, x$per))
  }
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
    table = data.frame(age = age, rate_effect = rate_contrib,
                       composition_effect = comp_contrib)
  )
  class(out) <- "dem_decomp"
  out
}

#' @export
print.dem_decomp <- function(x, ...) {
  cat("Kitagawa decomposition of a difference between two rates\n")
  cat(sprintf("  Crude rate 1: %.4f per %d\n", x$crude_rate1, x$per))
  cat(sprintf("  Crude rate 2: %.4f per %d\n", x$crude_rate2, x$per))
  cat(sprintf("  Total difference (1 - 2): %.4f\n", x$total))
  cat(sprintf("    Rate component:        %.4f\n", x$rate_component))
  cat(sprintf("    Composition component: %.4f\n", x$composition_component))
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
                         lx_col = "lx", Lx_col = "Lx", Tx_col = "Tx") {
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
