# DemogAnr 0.1.0

* Initial CRAN release.
* Added functions for basic demographic analysis: fertility (`dem.fert`), mortality (`dem.cdr`, `dem.mmr`, `dm.chm`), and life tables (`lifetable`, `lifetable_nqx`).
* Implemented mathematical population projection (`math_project`) and demographic interpolation models (`interpolation`, `karup_king`, `brass_logit`).
* Integrated Bayesian projection capabilities using `cmdstanr` and `brms`.
* Addressed CRAN reviewer feedback:
  - Formatted description and references in the `DESCRIPTION` file to match CRAN rules.
  - Suppressed all raw console printing in core demographic functions and implemented custom S3 print methods (`print.dem_cdr`, `print.dem_fert`, `print.dem_mmr`, and `print.dem_chm`) to align with standard R package behavior.
  - Added academic references to the documentation blocks of all core R files.
