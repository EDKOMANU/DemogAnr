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

# --- Deprecation -------------------------------------------------------------
# Warn the first time a deprecated name is used in a session, and then keep
# quiet: a warning on every call of a function used in a loop is noise, but
# saying nothing at all leaves the rename undiscovered.
.dep_warned <- new.env(parent = emptyenv())

.deprecate_once <- function(id, message) {
  if (is.null(.dep_warned[[id]])) {
    assign(id, TRUE, envir = .dep_warned)
    warning(message, call. = FALSE)
  }
  invisible(NULL)
}

# --- Coale-Demeny separation factors ------------------------------------------
# a0 and 4a1: the average person-years lived in the interval by those dying in
# it, at ages 0 and 1-4, as functions of the infant death rate 1M0 (Coale &
# Demeny 1983; Preston et al. 2001, Table 3.3). Shared by lifetable(), which
# has 1M0 to hand, and brass_lifetable(), which solves for it from 1q0.
.cd_a0_a1 <- function(m0, sex = c("male", "female")) {
  sex <- match.arg(sex)
  if (sex == "male") {
    c(a0 = if (m0 >= 0.107) 0.330 else 0.045 + 2.684 * m0,
      a1 = if (m0 >= 0.107) 1.352 else 1.651 - 2.816 * m0)
  } else {
    c(a0 = if (m0 >= 0.107) 0.350 else 0.053 + 2.800 * m0,
      a1 = if (m0 >= 0.107) 1.361 else 1.522 - 1.518 * m0)
  }
}
