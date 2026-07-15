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
