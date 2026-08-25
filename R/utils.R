# Small internal utilities shared across the package.

# --- Random-number hygiene ---------------------------------------------------
# The stochastic functions seed the generator so that a projection is
# reproducible, but a package must not leave the user's random stream
# disturbed. The pattern is:
#
#   if (!is.null(random_seed)) {
#     old_rng <- .capture_seed()
#     on.exit(.restore_seed(old_rng), add = TRUE)
#     set.seed(random_seed)
#   }
#
# so that whatever the caller was drawing before the call carries on
# unaffected afterwards.

# Current RNG state, or NULL if the stream has not been initialised yet.
.capture_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

# Put back a state captured by .capture_seed(). A NULL state means the stream
# was not initialised before the call, so the variable is removed again.
.restore_seed <- function(state) {
  if (is.null(state)) {
    suppressWarnings(rm(".Random.seed", envir = globalenv()))
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}
