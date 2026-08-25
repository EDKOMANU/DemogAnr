# DemogAnr

<!-- badges: start -->
[![R-CMD-check](https://github.com/EDKOMANU/DemogAnr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EDKOMANU/DemogAnr/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Demographic analysis in R: the classical methods of applied demography, with
the arithmetic written down once and checked against the published worked
examples.

Applied demographic work is still done largely in spreadsheets — a life table
here, a total fertility rate there, an age distribution split into single
years. Spreadsheets are familiar, but their formulas are invisible until each
cell is inspected, and reproducing an analysis a year later often means
rebuilding it from memory. `DemogAnr` puts the same methods behind named
functions that take a data frame and return a documented object.

## Installation

```r
# from CRAN
install.packages("DemogAnr")

# development version
# install.packages("remotes")
remotes::install_github("EDKOMANU/DemogAnr")
```

## A first look

```r
library(DemogAnr)

# A life table from deaths and population, Ghana 2010 census
data(gphc2010)
lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
lt$metrics$LifeExpectancyAtBirth
#> [1] 75.42111

# Every analysis function can attach a ggplot2 figure
plot(lt)
```

Fertility measures from a standard schedule:

```r
bangladesh <- data.frame(
  age         = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
  women       = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
  live_births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125),
  population  = 71315944
)
dem.fert(bangladesh, type = "all", age_col = "age",
         population_col = "population", women_col = "women",
         births_col = "live_births", graph = FALSE)$results$TFR
#> [1] 4.830224
```

## What is covered

| Area | Functions |
| --- | --- |
| Rates and ratios | `dem.fert()`, `dem.cdr()`, `dem.mmr()`, `dem.chm()` |
| Life tables | `lifetable()`, `lifetable_nqx()`, `multiple_decrement()`, `cause_deleted_lt()` |
| Model life tables | `model_lifetable()` (Coale-Demeny and UN), `brass_logit()`, `brass_lifetable()` |
| Indirect estimation | `chm_brass()` (Brass–Trussell), `pf_ratio()` (Brass P/F) |
| Registration completeness | `ggb()`, `seg()`, `ggb_seg()` |
| Reproduction | `reproduction()`, `stable_population()` |
| Data quality | `whipple()`, `myers()`, `age_ratio()`, `sex_ratio()`, `un_age_sex_accuracy()` |
| Interpolation | `karup_king()`, `interpolation()` |
| Standardization | `standardize()`, `decompose_rates()` (Kitagawa), `decompose_LE()` (Arriaga) |
| Projection | `cohort_component()`, `project_population()`, `math_project()`, `calibrate_regions()` |
| Marriage and health | `smam()` (Hajnal), `sullivan_hle()` |
| Figures | `pyramid()`, and `plot()` on any result built with `graph = TRUE` |

## Design notes

**Results, then figures.** Analysis functions return a list (or data frame)
carrying the full working — not just the headline number — so the intermediate
columns can be checked or re-tabulated. With `graph = TRUE` a **ggplot2**
figure is attached as `$plot` and returned by `plot()`. Printing a result never
draws to a device.

**Checked against the literature.** Where a method has a published worked
example, the test suite reproduces it: Manual X (Bangladesh 1974, Panama 1976),
Preston, Heuveline and Guillot (Sweden 1988 and 1993, US 1991), and Kpedekpo
(Ghana 1960). The documentation for each function cites the source it follows.

**Ghana data included.** `gphc2010`, `ghpop2021`, `ghmort2021` and `region2000`
carry census and vital-statistics extracts used throughout the examples.

## Documentation

```r
vignette("DemogAnr")   # worked tour of the package
help(package = "DemogAnr")
```

## References

Brass, W. (1975). *Methods for Estimating Fertility and Mortality from Limited
and Defective Data*. Chapel Hill: Carolina Population Center.

Preston, S. H., Heuveline, P., & Guillot, M. (2001). *Demography: Measuring and
Modeling Population Processes*. Oxford: Blackwell.

Siegel, J. S., & Swanson, D. A. (Eds.) (2004). *The Methods and Materials of
Demography* (2nd ed.). San Diego: Elsevier.

United Nations (1983). *Manual X: Indirect Techniques for Demographic
Estimation*. New York: United Nations.

## License

MIT © Edward Owusu Manu
