# DemogAnr 0.1.0

* Initial CRAN release.
* Added functions for basic demographic analysis: fertility (`dem.fert`), mortality (`dem.cdr`, `dem.mmr`, `dm.chm`), and life tables (`lifetable`, `lifetable_nqx`).
* Implemented mathematical population projection (`math_project`) and demographic interpolation models (`interpolation`, `karup_king`, `brass_logit`).
* Integrated Bayesian projection capabilities using `cmdstanr` and `brms`.
