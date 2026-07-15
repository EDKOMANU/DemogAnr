test_that("brass_logit recovers the standard exactly (alpha=0, beta=1)", {
  data(standards, envir = environment())
  # Cumulative probabilities q(x) = 1 - l(x) of the African standard itself
  d <- data.frame(age = standards$Age[-1], qx = 1 - standards$African[-1])
  fit <- brass_logit(d, qx_col = "qx", age_col = "age", standard = "African")
  expect_equal(unname(fit$coefficients[1]), 0, tolerance = 1e-3)
  expect_equal(unname(fit$coefficients[2]), 1, tolerance = 1e-3)
  # Round trip: predicted qx must equal the input qx
  expect_equal(fit$data$predicted_qx, d$qx, tolerance = 1e-3)
  expect_equal(fit$data$predicted_lx, 1 - d$qx, tolerance = 1e-3)
})

test_that("brass_logit fills gaps and warns on decreasing (interval) input", {
  data(standards, envir = environment())
  d <- data.frame(age = standards$Age[-1], qx = 1 - standards$African[-1])
  d$qx[c(3, 7, 11)] <- NA
  fit <- brass_logit(d, qx_col = "qx", age_col = "age", standard = "African")
  expect_false(any(is.na(fit$data$predicted_qx)))
  # decreasing series looks like interval nqx, not cumulative q(x)
  d2 <- data.frame(age = c(1, 5, 10, 15), qx = c(0.10, 0.05, 0.02, 0.01))
  expect_warning(
    brass_logit(d2, qx_col = "qx", age_col = "age", standard = "African"),
    "cumulative"
  )
})
