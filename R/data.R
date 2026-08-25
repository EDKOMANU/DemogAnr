#' Grouped population data by 5-year age groups
#'
#' Example population counts in 5-year age groups (0-4 to 75-79) for the years
#' 2021 to 2035, used to demonstrate [karup_king()].
#'
#' This dataset was called `data` before version 0.3.0. That name shadowed
#' base R's [utils::data()] function whenever it was loaded, so it was renamed;
#' see [data] for the deprecated original.
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
"grouped_pop"

#' Grouped population data by 5-year age groups (deprecated name)
#'
#' `data` is the former name of [grouped_pop], kept so that code written
#' against earlier versions keeps working. The name shadowed base R's
#' [utils::data()] function once loaded, which made `data(...)` calls in the
#' same script fail in confusing ways. Use [grouped_pop] instead; `data` will
#' be removed in a future release.
#'
#' @format Identical to [grouped_pop]: a data frame with 16 rows (age groups)
#'   and 16 variables.
#' @seealso [grouped_pop], which this duplicates.
#' @keywords internal
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
#' demonstrate [project_population()]. The three per-capita rate columns
#' (`cbr`, `cdr`, `nmr`) are the inputs to the balancing-equation projection;
#' they are illustrative values derived from the raw fields so that the example
#' produces realistic growth (roughly 1--3\% per year): `cbr` rises with
#' fertility (\eqn{0.005 \times TFR + 0.010}), `cdr` rescales `death_rate`
#' (\eqn{death\_rate / 7}), and `nmr` spreads the migration stock into an
#' annual per-capita rate (\eqn{net\_migration / base\_pop / 40}).
#'
#' @format A data frame with 16 rows and 9 variables:
#' \describe{
#'   \item{Country}{Character. Country name.}
#'   \item{Region}{Character. Region name.}
#'   \item{base_pop}{Numeric. Base-year population.}
#'   \item{TFR}{Numeric. Total fertility rate.}
#'   \item{death_rate}{Numeric. Raw (illustrative) death index.}
#'   \item{net_migration}{Numeric. Net migration stock (persons).}
#'   \item{cbr}{Numeric. Crude birth rate (per-capita annual rate).}
#'   \item{cdr}{Numeric. Crude death rate (per-capita annual rate).}
#'   \item{nmr}{Numeric. Net-migration rate (per-capita annual rate).}
#' }
#' @source Population and fertility figures adapted from Ghana Statistical
#'   Service census publications; the `cbr`, `cdr`, and `nmr` rate columns are
#'   illustrative, constructed for the projection example.
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

#' Ghana age-specific mortality by sex, 2021
#'
#' Age-specific population and deaths for Ghana in 2021, separately for males
#' and females, taken as the base year of a national age-specific mortality
#' projection. Because it is split by sex, this data set is suitable for
#' demonstrating sex-specific life tables and the decomposition of the
#' male-female gap in life expectancy with [decompose_LE()].
#'
#' @format A data frame with 36 rows (18 age groups x 2 sexes) and 5 variables:
#' \describe{
#'   \item{Sex}{Character. "Male" or "Female".}
#'   \item{Age}{Integer. Lower bound of the age group (0, 1, 5, 10, ..., 80).}
#'   \item{AgeGroup}{Character. Age-group label ("<1", "1-4", ..., "80+").}
#'   \item{Population}{Integer. Mid-year population (exposure) in the age group.}
#'   \item{Deaths}{Integer. Deaths in the age group.}
#' }
#' @source Ghana Statistical Service, age-specific mortality projection
#'   (2021 base year).
"ghmort2021"

#' Ghana single-year age distribution by sex, 2021
#'
#' Population by single year of age (0 to 79, with 80 representing the open
#' age group 80+), separately for males and females, Ghana 2021. Single-year
#' age counts are the input required by the digit-preference indices
#' [whipple()] and [myers()], and can be grouped into five-year age groups for
#' [un_age_sex_accuracy()].
#'
#' @format A data frame with 162 rows (81 single ages x 2 sexes) and 3
#'   variables:
#' \describe{
#'   \item{Sex}{Character. "Male" or "Female".}
#'   \item{Age}{Integer. Single year of age; 0 is age <1 and 80 is the open
#'     interval 80+.}
#'   \item{Population}{Integer. Population count.}
#' }
#' @source Ghana Statistical Service, 2021 Population and Housing Census
#'   (national single-year age distribution).
"ghpop2021"

#' Coale-Demeny and United Nations model life table rates
#'
#' Age-specific death rates \eqn{{}_nM_x} for the four Coale-Demeny regional
#' model life table families and the five United Nations patterns for
#' developing countries, by sex and level. Used by [model_lifetable()], which
#' builds a life table from them and is the usual way to reach this data.
#'
#' @format A data frame with 14,580 rows and 5 variables:
#' \describe{
#'   \item{family}{Factor. One of `CD_West`, `CD_North`, `CD_South`, `CD_East`,
#'     `UN_General`, `UN_Latin_American`, `UN_Chilean`, `UN_South_Asian`,
#'     `UN_Far_Eastern`.}
#'   \item{sex}{Factor, `male` or `female`.}
#'   \item{e0}{Numeric. Level of the table, as life expectancy at birth, from
#'     54 to 76 in steps of half a year.}
#'   \item{age}{Integer. Lower bound of the age group: 0, 1, 5, 10, ..., 80.}
#'   \item{nMx}{Numeric. Central death rate in the age group.}
#' }
#'
#' @details
#' The tabulated levels run from \eqn{e_0} of 54 to 76, so the very high
#' mortality levels of the published Coale-Demeny system are not covered here.
#' For a level outside that range, [brass_lifetable()] fits a relational model
#' that is not restricted to a tabulated grid.
#'
#' The rates stop at age 80. The value given there behaves like a five-year
#' group rate rather than an aggregate rate for the open interval, so
#' [model_lifetable()] does not close the table on it by default; see that
#' function's Details.
#'
#' @source Generated with the \pkg{MortCast} package from the model life
#'   tables published by the United Nations Population Division
#'   (\url{https://www.un.org/development/desa/pd/data/model-life-tables}).
#'   The underlying systems are Coale, A. J., & Demeny, P. (1983),
#'   \emph{Regional Model Life Tables and Stable Populations} (2nd ed.), New
#'   York: Academic Press; and United Nations (1982), \emph{Model Life Tables
#'   for Developing Countries}, Population Studies No. 77, New York: United
#'   Nations.
#'
#' @seealso [model_lifetable()]
"model_lt"
