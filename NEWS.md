# DemogAnr 0.3.0

This release extends the package from the everyday measures into the classical
models of formal demography: reproduction and the stable population, multiple
decrement and cause-deleted life tables, indirect fertility estimation, and
cohort-component projection. Wherever possible each new function reproduces the
published worked example from the standard texts, so that a reader can follow
the calculation in the book alongside the code.

## New features

* `reproduction()` computes the gross and net reproduction rates (GRR, NRR)
  with the total fertility rate, the mean age of childbearing, and an
  approximate intrinsic growth rate, from an age schedule of fertility and a
  female life table. Reproduces Preston et al. (2001) Box 5.5 (United States,
  1991: GRR 1.013, NRR 0.995).
* `stable_population()` solves Lotka's equation for the intrinsic growth rate
  by Coale's iterative procedure and derives the stable population: the
  intrinsic birth and death rates, the stable age distribution, and its mean
  age. Reproduces Preston et al. (2001) Box 7.1 (Egypt, 1997: intrinsic
  r = 0.01424).
* `multiple_decrement()` builds a cause-specific (multiple-decrement) life
  table, distributing the life-table deaths across causes and reporting the
  proportion of a birth cohort that will eventually die of each. Accepts either
  a population column, from which the all-cause life table is built, or a
  ready-made life table through `lx`, `nqx` and `nax_all`, so that the input can
  match whatever data is at hand. Reproduces Preston et al. (2001) Box 4.1 (US
  females 1991: 21,205 of 100,000 newborns eventually dying of neoplasms).
* `cause_deleted_lt()` builds the associated single-decrement (cause-deleted)
  life table, giving the life expectancy that would result if a cause were
  eliminated and the gain in e0. The cause-deletion assumption is the user's to
  choose through `method`, either `"proportional"` (proportional hazards, the
  default) or `"rates"` (direct subtraction of the cause-specific rates), with
  an optional user-supplied `nax`; it takes the same two input routes as
  `multiple_decrement()`. On Preston et al. (2001) Box 4.2 (US females 1991) it
  gives e0 = 78.92 rising to 82.43 without neoplasms, against their published
  82.46.
* `pf_ratio()` implements the Brass P/F ratio method for fertility estimation
  from children ever born and current births, with the Coale-Trussell
  interpolation coefficients of United Nations *Manual X*, Chapter II.
  Reproduces the *Manual X* Bangladesh 1974 worked example (K = 1.500, adjusted
  TFR = 7.24).
* `cohort_component()` performs cohort-component population projection through a
  Leslie matrix (Preston et al. 2001, Ch. 6), returning the matrix, the
  projected population by age and year, and the implied long-run growth rate.
  Reproduces Preston et al. (2001) Box 6.1 exactly (Sweden, females, 1993:
  4,449,570 in 1998 and 4,478,712 in 2003), and converges to the stable age
  structure and intrinsic growth rate given by `stable_population()`.
* `smam()` computes the singulate mean age at marriage from the proportions
  single by age (Hajnal 1953). Reproduces Preston et al. (2001) Box 4.4
  (Turkish males 1990: SMAM 25.0 years).
* `sullivan_hle()` computes Sullivan health expectancy, partitioning life
  expectancy into healthy and unhealthy years from an age-specific prevalence
  of ill-health.

## Documentation

* Examples throughout the package now use the published worked examples from
  the standard texts, so that a student can follow the calculation in the book
  alongside the code. `dem.fert()` uses the *Manual X* Bangladesh 1974 schedule
  (reproducing the published TFR of 4.83), `dem.cdr()` uses Preston et al.
  (2001) Box 1.2 (Sweden 1988, CDR 11.47 per 1000), and `lifetable_nqx()` uses
  the published US 1991 female probabilities of dying. The reproductions are
  covered by the test suite.
* The remaining examples have been moved off invented numbers onto the data
  bundled with the package: `dm.chm()` uses the Ghana 2021 child deaths,
  `dem.mmr()` the Ghana 2021 female population of reproductive age,
  `interpolation()` an intercensal estimation between the 2000 and 2010
  censuses, and `math_project()` the Ghana regional populations of
  `region2000`.
* A malformed roxygen line in `dem.fert()`'s documentation has been corrected,
  and `karup_king()`'s example no longer relies on the bundled data set that
  shares a name with base R's `data()` function, which could mask it.

## Changes to existing functions

* `project_population()` has been rebuilt on the demographic balancing
  equation. The annual growth rate is now `g = b - d + m` (crude birth rate
  minus crude death rate plus net-migration rate), and each simulated
  trajectory draws its three component rates **once** and holds them across the
  horizon, so that structural uncertainty persists and the projection fan
  widens with time. The distribution of every component is user-definable
  through the new `birth_dist`, `death_dist`, and `migration_dist` sampler
  arguments (defaulting to lognormal birth and death rates and a normal
  migration rate, with spread set by `cv`); an optional `annual_noise_sd` adds
  per-year innovations, and `probs` selects the summary quantiles. The previous
  fixed, ad hoc coefficient model has been removed.
  **Breaking change:** the inputs are now per-capita rates
  (`birth_rate_var`, `death_rate_var`, `net_migration_var`) rather than
  `TFR_var` and a migration count.
* The bundled `region2000` data set gains illustrative `cbr`, `cdr`, and `nmr`
  per-capita rate columns for the new `project_population()` interface.

# DemogAnr 0.2.0

* `lifetable()` gains a `nax_method` argument. In addition to the default
  Coale-Demeny young-age factors with `n/2` elsewhere (`"cd"`), it can now
  refine the average person-years lived in adult intervals with the iterative
  Keyfitz (1966) method (`"keyfitz"`). A user-supplied `nax` vector can also be
  passed to borrow separation factors from any external model life table
  system (e.g. Coale-Demeny or United Nations). Default behaviour is unchanged.
* New `chm_brass()`: indirect estimation of child mortality from children ever
  born and children surviving classified by age of mother, using the Trussell
  variant of the Brass method (United Nations Manual X, 1983). Supports all
  four Coale-Demeny model families (North, South, East, West), returns q(x) for
  x = 1, 2, 3, 5, 10, 15, 20 with reference periods and dates. Reproduces the
  Manual X Panama (1976) worked example to the published precision.
* New `standardize()`: compares two populations by age standardization
  following Preston et al. (2001, Chapter 2), using the average of the two age
  compositions as the default standard and reporting the crude and
  age-standardized rates plus the comparative mortality ratio (CMR). Reproduces
  their Box 2.1 (Sweden vs Kazakhstan) exactly.
* New `decompose_rates()`: Kitagawa decomposition of a difference between two
  crude rates into rate and age-composition components.
* New `decompose_LE()`: Arriaga decomposition of a difference in life
  expectancy at birth into age-specific direct and indirect/interaction
  effects.
* Visual output: `dem.fert`, `dem.cdr`, `lifetable`, `project_population`,
  `decompose_LE`, `myers`, `whipple` and `un_age_sex_accuracy` gain a
  `graph = TRUE` argument. When enabled, a 'ggplot2' plot of the analysis is
  attached to the result (`$plot`) and shown when the object is printed;
  `plot()` methods are provided for each. Scripts that assign the result draw
  nothing, so the behaviour is safe in batch code.
* New `pyramid()`: draws a population pyramid (back-to-back bar chart by age
  group and sex) with 'ggplot2'.
* New age-quality / age-heaping indices: `whipple()` (Whipple's index of
  digit preference for 0 and 5), `myers()` (Myers' blended index over all ten
  terminal digits), `sex_ratio()` and `age_ratio()` helpers, and
  `un_age_sex_accuracy()` (the United Nations joint age-sex accuracy index).
  `un_age_sex_accuracy()` reproduces the Ghana 1960 worked example of Kpedekpo
  (1982); Myers follows Rodriguez (2015) / Myers (1940).
* New bundled data set `ghpop2021`: Ghana 2021 single-year age distribution by
  sex, for the digit-preference indices.
* New bundled data set `ghmort2021`: age-specific population and deaths for
  Ghana in 2021 by sex, for demonstrating sex-specific life tables and the
  decomposition of the male-female gap in life expectancy.

# DemogAnr 0.1.0

* Initial CRAN release.
* Added functions for basic demographic analysis: fertility (`dem.fert`), mortality (`dem.cdr`, `dem.mmr`, `dm.chm`), and life tables (`lifetable`, `lifetable_nqx`).
* Implemented mathematical population projection (`math_project`), a stochastic simulation-based projection (`project_population`), and demographic interpolation models (`interpolation`, `karup_king`, `brass_logit`).
* Addressed CRAN reviewer feedback:
  - Formatted description and references in the `DESCRIPTION` file to match CRAN rules.
  - Suppressed all raw console printing in core demographic functions and implemented custom S3 print methods (`print.dem_cdr`, `print.dem_fert`, `print.dem_mmr`, and `print.dem_chm`) to align with standard R package behavior.
  - Added academic references to the documentation blocks of all core R files.
  - All invalid links have been removed from the documentation and Description file.
* Methodological corrections prior to release:
  - `brass_logit()` now uses the correct inverse of the Brass half-logit transform, exp(2Y)/(1 + exp(2Y)), and documents that its input is the cumulative probability of dying q(x) = 1 - l(x).
  - `lifetable_nqx()` and `lifetable()` now account for the width of each age interval (1, 4, 5, ... years) when computing nqx and person-years lived (nLx); `lifetable()` gains `sex` and `radix` arguments and uses Coale-Demeny separation factors for ages 0 and 1-4.
  - `karup_king()` defaults now work out of the box (packaged Karup-King coefficient matrices).
  - `project_population()` now propagates all simulated trajectories through the projection horizon and reports lower quartile, median, mean, and upper quartile per year; net out-migration now reduces growth.
  - `dm.chm()` neonatal filter now uses age < 1 month.
  - Removed unused dependencies (tidyverse, brms, logger, data.table) and leftover Stan scaffolding.
