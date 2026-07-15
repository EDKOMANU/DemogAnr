#' Grouped population data by 5-year age groups
#'
#' Example population counts in 5-year age groups (0-4 to 75-79) for the years
#' 2021 to 2035, used to demonstrate [karup_king()].
#'
#' @format A data frame with 16 rows (age groups) and 16 variables:
#' \describe{
#'   \item{2021}{Numeric. Population in the age group in 2021.}
#'   \item{2022}{Numeric. Population in the age group in 2022.}
#'   \item{2023}{Numeric. Population in the age group in 2023.}
#'   \item{2024}{Numeric. Population in the age group in 2024.}
#'   \item{2025}{Numeric. Population in the age group in 2025.}
#'   \item{2026}{Numeric. Population in the age group in 2026.}
#'   \item{2027}{Numeric. Population in the age group in 2027.}
#'   \item{2028}{Numeric. Population in the age group in 2028.}
#'   \item{2029}{Numeric. Population in the age group in 2029.}
#'   \item{2030}{Numeric. Population in the age group in 2030.}
#'   \item{2031}{Numeric. Population in the age group in 2031.}
#'   \item{2032}{Numeric. Population in the age group in 2032.}
#'   \item{2033}{Numeric. Population in the age group in 2033.}
#'   \item{2034}{Numeric. Population in the age group in 2034.}
#'   \item{2035}{Numeric. Population in the age group in 2035.}
#'   \item{age_col}{Character. The 5-year age group, e.g. "0-4", "5-9", ..., "75-79".}
#' }
#' @source Simulated data based on population structures from the Ghana
#'   Statistical Service.
"data"

#' Karup-King coefficients for the first age group
#'
#' The standard Karup-King third-difference osculatory interpolation
#' multipliers used to split the *first* 5-year age group into single years.
#' Rows are the five single-year ages within the group; columns are the
#' weights applied to the first, second, and third 5-year groups.
#'
#' @format A 5 x 3 data frame of numeric coefficients.
#' @source Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods
#'   and Materials of Demography} (2nd ed.), Appendix C.
"first_coef"

#' Ghana 2010 Population and Housing Census mortality data
#'
#' Population and deaths by age group from the Ghana 2010 Population and
#' Housing Census, used to demonstrate [lifetable()].
#'
#' @format A data frame with 21 rows and 3 variables:
#' \describe{
#'   \item{Age}{Numeric. Lower bound of the age group (0, 1, 5, 10, ...).}
#'   \item{Pop}{Numeric. Population in the age group.}
#'   \item{Deaths}{Numeric. Deaths in the age group.}
#' }
#' @source Ghana Statistical Service, 2010 Population and Housing Census.
"gphc2010"

#' Karup-King coefficients for the last age group
#'
#' The standard Karup-King third-difference osculatory interpolation
#' multipliers used to split the *last* 5-year age group into single years.
#' Rows are the five single-year ages within the group; columns are the
#' weights applied to the antepenultimate, penultimate, and last 5-year groups.
#'
#' @format A 5 x 3 data frame of numeric coefficients.
#' @source Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods
#'   and Materials of Demography} (2nd ed.), Appendix C.
"last_coef"

#' Karup-King coefficients for middle age groups
#'
#' The standard Karup-King third-difference osculatory interpolation
#' multipliers used to split *interior* 5-year age groups into single years.
#' Rows are the five single-year ages within the group; columns are the
#' weights applied to the previous, current, and next 5-year groups.
#'
#' @format A 5 x 3 data frame of numeric coefficients.
#' @source Siegel, J. S., & Swanson, D. A. (Eds.). (2004). \emph{The Methods
#'   and Materials of Demography} (2nd ed.), Appendix C.
"middle_coef"

#' Regional projection parameters, Ghana 2000
#'
#' Base-year demographic parameters for the regions of Ghana in 2000, used to
#' demonstrate [project_population()].
#'
#' @format A data frame with 16 rows and 6 variables:
#' \describe{
#'   \item{Country}{Character. Country name.}
#'   \item{Region}{Character. Region name.}
#'   \item{base_pop}{Numeric. Base-year population.}
#'   \item{TFR}{Numeric. Total fertility rate.}
#'   \item{death_rate}{Numeric. Crude death rate.}
#'   \item{net_migration}{Numeric. Net migration (persons).}
#' }
#' @source Derived from Ghana Statistical Service census publications.
"region2000"

#' Brass standard life tables
#'
#' Survivorship values \eqn{l(x)} (radix 1) and their Brass logits for the
#' African and General United Nations standards, used by [brass_logit()].
#'
#' @format A data frame with 24 rows and 5 variables:
#' \describe{
#'   \item{Age}{Numeric. Exact age x.}
#'   \item{GeneralUN}{Numeric. l(x) of the General UN standard.}
#'   \item{GeneralUN_logit}{Numeric. Brass logit of the General UN standard.}
#'   \item{African}{Numeric. l(x) of the African standard.}
#'   \item{African_logit}{Numeric. Brass logit of the African standard.}
#' }
#' @source Brass, W. (1971). On the scale of mortality. In \emph{Biological
#'   Aspects of Demography}. London: Taylor & Francis. United Nations model
#'   life table publications.
"standards"

#' District-level demographic test data
#'
#' Simulated district-level demographic panel data (population, births,
#' deaths, migration) by age group and sex, for testing subnational
#' projection workflows such as [math_project()].
#'
#' @format A data frame with 1224 rows and 10 variables:
#' \describe{
#'   \item{region}{Character. Region name.}
#'   \item{district}{Character. District name.}
#'   \item{year}{Numeric. Calendar year.}
#'   \item{age_group}{Character. Age group label.}
#'   \item{sex}{Character. "Male" or "Female".}
#'   \item{population}{Numeric. Population count.}
#'   \item{births}{Numeric. Live births.}
#'   \item{deaths}{Numeric. Deaths.}
#'   \item{in_migration}{Numeric. In-migrants.}
#'   \item{out_migration}{Numeric. Out-migrants.}
#' }
#' @source Simulated data.
"testdata"
