## Test environments
* local Windows 11 install, R 4.5.3
* win-builder (release and devel)

## R CMD check results
There were 0 errors, 0 warnings, and 1 note.

* checking for future file timestamps ... NOTE
  unable to verify current time (standard local check artifact)

## Submission Notes & Fixes
This is a resubmission addressing the reviewers' feedback:

1. **DESCRIPTION File:**
   - Rewrote the `Description` field to not start with the package name or "This package".
   - Added properly formatted methodological references in the form `authors (year) <doi:...>` or `authors (year, ISBN:...)` with no spaces after prefixes.
   - Removed unused dependencies (tidyverse, brms, logger, data.table).

2. **R Function References:**
   - Added comprehensive academic references (`@references`) to the roxygen2 documentation headers of R files to document the methods.

3. **Console Logging & Print/Cat Suppression:**
   - Removed direct console logging (`cat()` / `print()`) inside core computation functions.
   - Wrapped the outputs of `dem.cdr`, `dem.fert`, `dem.mmr`, and `dm.chm` in custom S3 classes with exported `print` methods.

4. **Methodological corrections:**
   - Corrected the inverse Brass logit transform in `brass_logit()`.
   - Corrected age-interval handling (widths 1, 4, 5, ...) in `lifetable()` and `lifetable_nqx()`.
   - `project_population()` now propagates full simulation trajectories and returns the documented quantile summaries.
   - Added a unit test suite covering all exported functions.
