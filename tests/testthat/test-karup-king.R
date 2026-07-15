test_that("karup_king works with packaged default coefficients", {
  data(data, envir = environment())
  out <- karup_king(df = data, age_col = "age_col", pops = "2021")
  # 16 groups of 5 years -> 80 single ages (0..79)
  expect_equal(nrow(out), 80)
  # totals approximately preserved (rounding to whole persons)
  expect_lt(abs(sum(out[["2021"]]) - sum(data[["2021"]])), nrow(out))
})

test_that("karup_king rejects open-ended age groups", {
  d <- data.frame(age_group = c("0-4", "5-9", "10+"), population = c(100, 90, 80))
  expect_error(karup_king(d, age_col = "age_group", pops = "population"), "Open-ended|form")
})

test_that("karup_king coefficient rows sum to 0.2 (fifths of the group)", {
  data(first_coef, envir = environment()); data(middle_coef, envir = environment())
  data(last_coef, envir = environment())
  expect_equal(unname(rowSums(as.matrix(first_coef))), rep(0.2, 5))
  expect_equal(unname(rowSums(as.matrix(middle_coef))), rep(0.2, 5))
  expect_equal(unname(rowSums(as.matrix(last_coef))), rep(0.2, 5))
})
