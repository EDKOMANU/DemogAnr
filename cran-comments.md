## Test environments

* local Windows 11 install, R 4.5.3
* win-builder (release and devel)
* Ubuntu 22.04, R 4.3.3

## R CMD check results

There were 0 errors, 0 warnings, and 1 note.

* checking for future file timestamps ... NOTE
  unable to verify current time (standard local check artifact)

## Submission notes

This is an update of an existing package, from version 0.2.0 to 0.3.0.

It extends the package from the everyday rate and life table measures into
the classical demographic models: multiple-decrement and cause-deleted life
tables, the Brass P/F ratio method for fertility, reproduction measures and
the stable population model, cohort-component projection, the singulate mean
age at marriage, Sullivan health expectancy, and the calibration of
subnational projections to a national control total.

It also adds the Coale-Demeny and United Nations model life table systems
(`model_lifetable()`, with the rates shipped as the `model_lt` dataset), a
life table generated from Brass relational logit parameters
(`brass_lifetable()`), and the three death distribution methods for estimating
the completeness of death registration (`ggb()`, `seg()`, `ggb_seg()`), and
indirect adult mortality from maternal orphanhood (`orphanhood()`).
Where a method has a published worked example the test suite reproduces it,
including from United Nations *Manual X*, Preston, Heuveline and Guillot
(2001), and Kpedekpo (1982).

### Deprecations

Two names are deprecated in this version. Both continue to work and return
exactly what they did before, so code written against version 0.2.0 is
unaffected:

* The bundled dataset `data` is now `grouped_pop`. The former name masked
  base R's `utils::data()` function once loaded, which made subsequent
  `data()` calls in the same script fail confusingly. The old dataset is
  still shipped and is identical to the new one; it is documented with
  `@keywords internal` and will be removed in a future release.

* `dm.chm()` is now `dem.chm()`, for consistency with `dem.cdr()`,
  `dem.fert()` and `dem.mmr()`. The old name delegates to the new one and
  warns once per session.

`dem.mmr()` gains two clearly named elements in its result, `MMRate` and
`MMRatio`. The existing `MMR` and `MMR_Ratio` are retained at their original
values rather than redefined, so that no existing result changes silently.

### Other notes

* `print()` methods no longer draw the attached figure. Figures are returned
  by `plot()`, so nothing is written to the user's working directory and no
  graphics device is opened as a side effect.

* `project_population()` and `calibrate_regions()` seed the random number
  generator for reproducibility, and restore the user's `.Random.seed` on
  exit so that the calling session's random stream is left undisturbed.
  Passing `random_seed = NULL` skips seeding entirely.

* All examples run in well under five seconds each and none are wrapped in
  `\donttest{}`.

* This package has no reverse dependencies.
