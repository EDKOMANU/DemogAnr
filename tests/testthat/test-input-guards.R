# Inputs that used to be accepted and quietly produce a wrong answer.

test_that("decompose_LE() refuses life tables with different radices", {
  mk <- function(radix) {
    lifetable(data.frame(Age = c(0, 1, seq(5, 85, 5)),
                         Pop = rep(1000, 19), Deaths = rep(10, 19)),
              age = "Age", pop = "Pop", Dx = "Deaths",
              radix = radix, graph = FALSE)$lifetable
  }
  big <- mk(100000); small <- mk(1000)

  # Identical mortality, so any honest decomposition is zero; before the guard
  # this returned a difference of -99 years.
  expect_error(decompose_LE(big, small, graph = FALSE), "same radix")

  same <- decompose_LE(big, mk(100000), graph = FALSE)
  expect_equal(same$total, 0)
})

test_that("myers() warns when the data do not cover the blending range", {
  heaped <- data.frame(Age = 0:99,
                       Pop = round(100000 * exp(-0.03 * (0:99))))
  heaped$Pop[heaped$Age %% 10 == 0] <- round(heaped$Pop[heaped$Age %% 10 == 0] * 1.3)
  short <- subset(heaped, Age <= 80)

  expect_silent(myers(heaped, "Age", "Pop", graph = FALSE))
  expect_warning(myers(short, "Age", "Pop", graph = FALSE),
                 "unbalanced")
  # the warning names a range the data can actually support
  expect_warning(myers(short, "Age", "Pop", graph = FALSE),
                 "upper = 71")
  # Taking the advice restores a balanced blend: it lands next to the
  # untruncated series, whereas the unbalanced default is far off.
  full      <- myers(heaped, "Age", "Pop", graph = FALSE)$index
  corrected <- myers(short, "Age", "Pop", upper = 71, graph = FALSE)$index
  unbalanced <- suppressWarnings(myers(short, "Age", "Pop", graph = FALSE)$index)
  expect_equal(corrected, full, tolerance = 0.01)
  expect_gt(abs(unbalanced - full), 10 * abs(corrected - full))
})

test_that("karup_king() rejects age groups that are not contiguous five-year groups", {
  ten <- data.frame(age_col = c("0-9", "10-19", "20-29", "30-39"),
                    pop = c(4000, 3000, 2000, 1000))
  expect_error(karup_king(ten, age_col = "age_col", pops = "pop"),
               "five-year age groups")

  gappy <- data.frame(age_col = c("0-4", "5-9", "15-19", "20-24"),
                      pop = c(4000, 3000, 2000, 1000))
  expect_error(karup_king(gappy, age_col = "age_col", pops = "pop"),
               "contiguous")

  ok <- data.frame(age_col = c("0-4", "5-9", "10-14", "15-19"),
                   pop = c(4000, 3000, 2000, 1000))
  expect_s3_class(karup_king(ok, age_col = "age_col", pops = "pop"),
                  "data.frame")
})

test_that("stable_population() warns when Lotka's equation has not converged", {
  d <- data.frame(age = seq(15, 45, 5), Lx = rep(90000, 7),
                  fx = c(0.02, 0.15, 0.16, 0.10, 0.05, 0.01, 0.001))

  expect_warning(s <- stable_population(d, "age", "fx", "Lx", max_iter = 2,
                                        tol = 1e-12, graph = FALSE),
                 "did not converge")
  expect_false(s$converged)

  s2 <- stable_population(d, "age", "fx", "Lx", tol = 1e-12, graph = FALSE)
  expect_true(s2$converged)
  expect_lt(s2$iterations, 100)
})

test_that("pf_ratio() rejects k_ages that match no age group", {
  bd <- data.frame(
    age    = seq(15, 45, 5),
    women  = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
    ceb    = c(1160919, 4901382, 9085852, 9910256, 10384001, 9164329, 6905673),
    births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125))

  expect_error(pf_ratio(bd, "age", "women", "ceb", "births",
                        k_ages = c(21, 26), graph = FALSE),
               "matched no age group")

  # a partly-valid k_ages warns but still uses the values that do match
  expect_warning(p <- pf_ratio(bd, "age", "women", "ceb", "births",
                               k_ages = c(20, 25, 31), graph = FALSE),
                 "not age-group lower bounds")
  expect_equal(p$K, mean((p$table$P / p$table$F)[p$table$age %in% c(20, 25)]))

  # the published Manual X answer is unaffected
  ok <- pf_ratio(bd, "age", "women", "ceb", "births", graph = FALSE)
  expect_equal(round(ok$K, 2), 1.50)
  expect_false(is.nan(ok$TFR_adjusted))
})
