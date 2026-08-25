# DemogAnr 0.3.0

## New features

### Model life tables and the completion of indirect estimates

* `model_lifetable()` returns a life table from the four Coale-Demeny regional
  families or the five United Nations patterns for developing countries,
  indexed either by level (`e0`) or, more usefully, by an observed child
  mortality (`q1` or `q5`). The tabulated rates ship as `model_lt`: nine
  families by two sexes by 45 levels of `e0` from 54 to 76, in half-year
  steps, generated with the MortCast package from the United Nations
  Population Division's published tables.

  The rates stop at age 80, and the value given there behaves like a
  five-year group rate rather than an aggregate rate for the open interval.
  Closing the table on it would let too much of the cohort survive above 80,
  overshooting the level's own label by up to 1.2 years. The label is itself
  published information, so it pins the open interval down exactly, and
  `calibrate_open = TRUE` (the default) uses it: the returned table
  reproduces its level exactly and implies a remaining life expectancy at 80
  of 5.4 to 6.9 years, which is the order real life tables show. Pass
  `calibrate_open = FALSE` for the uncalibrated figure.

  This closes a loop the package could not previously close. `chm_brass()`
  estimates q(x) for a named Coale-Demeny family but stops there; passing its
  q(5) to `model_lifetable()` now supplies the whole age pattern, and life
  expectancy with it.

* `brass_lifetable()` builds a complete life table from a pair of Brass
  relational logit parameters. Where `brass_logit()` estimates alpha and beta
  from observed survivorship and predicts only at the ages supplied, this runs
  the model the other way, across the standard's full age range. It takes a
  `brass_logit()` fit directly, or alpha and beta given by hand, or a
  `target_e0` or `target_q5` to solve the level for -- so a table can be
  generated with no data at all.

* `brass_logit()` gains `complete = TRUE`, which carries the fit straight
  through to that life table and attaches it as `$lifetable`, making the
  common path a single call. Its result is now an object of class
  `dem_brass_fit` recording `alpha`, `beta` and the `standard` used.

### Indirect adult mortality

* `orphanhood()` estimates adult female survivorship from the proportion of
  respondents whose mother is still alive, by the Brass method of *Manual X*
  Chapter IV. It is the adult counterpart of `chm_brass()`: that estimates
  child mortality from a mother's report on her children, this estimates the
  mother's own mortality from her children's report on her. It returns the
  conditional survivorship l(25+n)/l(25) for each age group of respondent,
  together with the number of years before the survey each estimate refers
  to, so that a single survey traces a mortality trend backwards.

  The weighting factors of *Manual X* table 86 and the standard function of
  table 88 are carried in the package. The whole chain reproduces the
  published Bolivia 1975 worked example exactly: all seven weighting factors
  to four decimal places (tables 91), all seven survivorship ratios to the
  three decimals published, and all six time references to the tenth of a
  year published (table 92). That example is a test.

### Completeness of death registration

* `ggb()`, `seg()` and `ggb_seg()` estimate how completely deaths are
  registered relative to an enumerated population, by the generalized growth
  balance, synthetic extinct generations, and the combined procedure. Where
  registration is incomplete a life table built from registered deaths
  understates mortality, and these give the factor to correct it by; `ggb()`
  also estimates how the coverage of two censuses differs. Each attaches a
  diagnostic plot, which should be inspected before the estimate is used.

  All three are validated in the test suite against a stable population built
  by numerical integration, so the true completeness is known exactly: they
  recover it to within 1 per cent over completeness from 0.3 to 1.0.

* A `README`, and a GitHub Actions `R-CMD-check` workflow covering macOS,
  Windows, and Linux on R-devel, release, oldrel and the declared minimum
  R 4.1.

### Documentation

* The vignette now covers the indirect estimation pipeline end to end: child
  mortality from `chm_brass()` completed into a life table with
  `model_lifetable()`, the Brass relational alternative through
  `brass_logit(complete = TRUE)` and `brass_lifetable()`, adult mortality from
  `orphanhood()`, and a section on assessing the completeness of death
  registration with `ggb()`, `seg()` and `ggb_seg()`.

## Renamed, with the old names kept working

* The bundled dataset `data` is now `grouped_pop`. Its former name shadowed
  base R's `data()` function as soon as it was loaded, which made `data(...)`
  calls later in the same script fail confusingly. `data` is still shipped and
  identical, and will be removed in a future release.

* `dm.chm()` is now `dem.chm()`, matching `dem.cdr()`, `dem.fert()` and
  `dem.mmr()`. `dm.chm()` still works, returns exactly the same value, and
  warns once per session.

* `dem.mmr()` gained `results$MMRate` (deaths per 100,000 women) and
  `results$MMRatio` (deaths per 100,000 live births). Note that the existing
  `results$MMR` holds the *rate*, which is the opposite of the near-universal
  convention that "MMR" means the ratio. Rather than change what `MMR` returns
  and silently alter existing results, the unambiguous pair was added
  alongside it; `MMR` and `MMR_Ratio` keep their original values and will be
  removed in a future release. `print()` now shows both, and no longer dumps
  the whole input data frame.

* `calibrate_regions()` is a new function for the disaggregation step that
  follows a national projection. Where `project_population()` projects each
  region independently from the balancing equation, `calibrate_regions()` takes
  an already-settled national projection as a control total and asks only how it
  is split between the regions: each region's base-year profile is allowed to
  deviate from the national one by a user-specified distribution, multiplied by
  the regional population, and then raked back onto the national figure. This is
  the pattern used in official projection work, where the national series is
  agreed first and the regional series must add up to it.

  Two things it does that independent projection cannot. It rakes: the regions
  sum exactly to the national control total in every year, and in every age
  group when `age_col` is supplied, with the raking factors reported in the
  output so the adjustment is auditable rather than silent. And its deviation is
  correlated: the regions move around a shared national path instead of
  wandering independently, since independent regional errors cancel on
  aggregation and understate uncertainty at the top.

  The deviation may be given as a single coefficient of variation, a per-region
  vector, or a sampler `function(mean, n)` in the same form as the samplers of
  `project_population()`; it is drawn once per region and trajectory and phased
  in over the horizon, so the base year is reproduced exactly and the fan widens
  with time. With no deviation the function reduces to the ordinary
  constant-share allocation and every raking factor is 1. The result is an
  object of class `dem_calibration` whose `print` method shows the full working
  table (base population, share, control total, growth factor, deviation,
  unraked allocation, raking factor, calibrated population and interval),
  the raking factors, and the additivity check.

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

* A package vignette, `vignette("DemogAnr")`, walks through a complete
  demographic workflow: rates, life tables, multiple-decrement and
  cause-deleted tables, reproduction and the stable population, indirect
  estimation, standardization and decomposition, projection, and data-quality
  assessment. Each section reproduces the published worked example for the
  method it demonstrates.
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
* `project_population()` also gains target rates. `birth_target`,
  `death_target` and `migration_target` set the value each component is assumed
  to reach at `future_year` — as a single number, or as a column giving a
  separate target for each region — and `path` chooses whether the rate moves
  there linearly or exponentially. Each trajectory draws both a start rate and
  an end rate, so the assumed trend and the uncertainty about it are carried
  together; `target_cv` sets the spread around the target, and is normally
  wider than `cv` because a rate two decades out is less well known than
  today's. Without a target the component is held constant and the previous
  behaviour is reproduced exactly. This turns the function from a
  constant-rate scenario into a projection in the usual demographic sense,
  where the analyst states where the rates are going.
* The bundled `region2000` data set gains illustrative `cbr`, `cdr`, and `nmr`
  per-capita rate columns for the new `project_population()` interface.

## Bug fixes

* `project_population()` now returns an object of class `dem_projection`
  whether or not `graph = TRUE`. Previously `graph = FALSE` returned a bare
  data frame, so `plot()` fell through to `plot.data.frame()` and drew a
  scatterplot matrix instead of raising an error.

* `project_population()` no longer collapses migration uncertainty to zero
  where the net migration rate is zero. A coefficient of variation is
  degenerate for a signed quantity (`cv * |0| = 0`), so a region with balanced
  migration carried no migration uncertainty at all. The spread is now an
  absolute standard deviation, settable with the new `migration_sd` argument
  and otherwise borrowed from the region's own birth and death rates, with a
  warning naming the regions concerned.

* The projection figure's title now names the interval it actually drew (for
  example "80% interval" for the default `probs = c(0.1, 0.9)`), instead of
  always claiming an interquartile band.

* `print()` methods no longer draw the attached figure. Printing a result in a
  non-interactive session opened a graphics device and wrote an `Rplots.pdf`
  into the working directory. Figures are unchanged and still reached with
  `plot()`, or as the `$plot` element.

* `project_population()` and `calibrate_regions()` now restore the caller's
  random stream when they return. Previously they called `set.seed()`
  unconditionally, so any simulation running around them silently changed
  its draws. Results for a given `random_seed` are unchanged; passing
  `random_seed = NULL` now skips seeding altogether.

* `decompose_LE()` now checks that the two life tables share a radix. Given
  tables built on different radices it returned a difference dominated by the
  radix ratio -- two tables with identical mortality reported a gap of 99
  years -- rather than the intended decomposition.

* `myers()` now warns when the data do not cover every single year of age from
  `lower` to `upper + 9`. Without the full span the blending weights are
  unbalanced and the index is biased; the warning names an `upper` the data
  can support.

* `karup_king()` now rejects age groups that are not contiguous five-year
  groups. Ten-year groups previously produced a distribution with half the age
  range missing while the column totals still matched the input.

* `stable_population()` now warns when the iterative solution of Lotka's
  equation stops at `max_iter` without converging, and reports `converged` in
  the returned object. The provisional `r` was previously returned in silence.

* `pf_ratio()` now raises an error when `k_ages` matches no age group (which
  produced a silent `NaN` for `K` and `TFR_adjusted`) and warns about
  individual `k_ages` values that are not age-group lower bounds.

This release extends the package from the everyday measures into the classical
models of formal demography: reproduction and the stable population, multiple
decrement and cause-deleted life tables, indirect fertility estimation,
cohort-component projection, and the regional disaggregation of a national
projection. Wherever possible each new function reproduces the published worked
example from the standard texts, so that a reader can follow the calculation in
the book alongside the code.

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
