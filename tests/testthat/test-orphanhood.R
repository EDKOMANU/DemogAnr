# Validated against the published worked example in United Nations Manual X,
# Chapter IV: Bolivia, 1975 (data in table 89, results in tables 91 and 92).

bolivia <- data.frame(
  age   = seq(15, 50, 5),
  alive = c(5540, 3995, 2886, 1910, 1272, 855, 541, 234),
  # the counts with mother dead are back-derived from the published
  # proportions and rounded to whole respondents, so S(n) is recovered only
  # to about a unit in the fourth decimal place
  dead  = c(448, 541, 769, 852, 945, 1059, 985, 788),
  # the ten-year proportions as published in table 92
  s10   = c(NA, 0.9060, 0.8401, 0.7474, 0.6313, 0.5175, 0.3996, NA)
)

test_that("the proportions with mother alive match table 89", {
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, graph = FALSE)
  # S(n) for the age groups being estimated, from table 89 column (3)
  expect_equal(round(o$table$S, 3),
               round(c(0.8807, 0.7896, 0.6915, 0.5737, 0.4467, 0.3546,
                       0.2289), 3))
})

test_that("the weighting factors reproduce Manual X table 91 exactly", {
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, graph = FALSE)
  expect_equal(o$table$age, seq(20, 50, 5))
  expect_equal(round(o$table$W, 4),
               c(0.9874, 1.0976, 1.1784, 1.2416, 1.2414, 1.2126, 1.0540))
})

test_that("the survivorship ratios reproduce Manual X table 91", {
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, graph = FALSE)
  expect_equal(round(o$table$lx_ratio, 3),
               c(0.925, 0.890, 0.807, 0.720, 0.604, 0.466, 0.361))
  # survivorship must fall as the age of the respondent rises
  expect_true(all(diff(o$table$lx_ratio) < 0))
})

test_that("the time references reproduce Manual X table 92", {
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, survey_year = 1975.5, graph = FALSE)
  # table 92 column (6); the oldest group has no entry because Z(x) is
  # tabulated only to age 75 and M + n = 78.8 there
  expect_equal(round(o$table$t[1:6], 1),
               c(8.6, 10.3, 11.8, 13.2, 14.3, 15.0))
  expect_true(is.na(o$table$t[7]))
  # and the intermediate u(n), table 92 column (5). Manual X carries Z(x) to
  # three decimals in its own arithmetic, so the published 0.1405 and the
  # 0.1402 obtained from the interpolated Z differ in the fourth place; both
  # give t = 8.6, which is the figure the method actually reports.
  expect_equal(round(o$table$u[1], 3), 0.140)

  expect_equal(round(o$table$ref_date[1:2], 1), c(1966.9, 1965.2))
  expect_true(is.na(o$table$ref_date[7]))
})

test_that("estimates refer to 8 to 15 years before the survey", {
  # Manual X notes this as the characteristic range of the method
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, graph = FALSE)
  tt <- o$table$t[!is.na(o$table$t)]
  expect_true(all(tt > 8 & tt < 16))
  expect_true(all(diff(tt) > 0))    # older respondents look further back
})

test_that("proportions can be supplied instead of counts", {
  b2 <- data.frame(age = seq(15, 50, 5),
                   S = c(0.9252, 0.8807, 0.7896, 0.6915, 0.5737, 0.4467,
                         0.3546, 0.2289),
                   s10 = bolivia$s10)
  o <- orphanhood(b2, "age", prop = "S", prop10 = "s10", M = 28.8,
                  graph = FALSE)
  expect_equal(round(o$table$W, 4),
               c(0.9874, 1.0976, 1.1784, 1.2416, 1.2414, 1.2126, 1.0540))
  expect_equal(round(o$table$lx_ratio, 3),
               c(0.925, 0.890, 0.807, 0.720, 0.604, 0.466, 0.361))
})

test_that("without prop10 the ten-year proportion is built from the counts", {
  o <- orphanhood(bolivia, "age", "alive", "dead", M = 28.8, graph = FALSE)
  # the youngest groups, where the counts are those published, match closely
  expect_equal(round(o$table$S10[1:3], 3), c(0.906, 0.840, 0.747))
  # the survivorship ratios do not depend on the ten-year proportion at all
  ref <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                    M = 28.8, graph = FALSE)
  expect_equal(o$table$lx_ratio, ref$table$lx_ratio)
})

test_that("the weighting factors vary with M as the table says", {
  # W rises with the mean age at maternity in every age group
  Ws <- vapply(c(22, 25, 28, 30), function(m) {
    orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10", M = m,
               graph = FALSE)$table$W
  }, numeric(7))
  expect_true(all(apply(Ws, 1, function(r) all(diff(r) > 0))))

  # a whole-number M picks the tabulated column without interpolation
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10", M = 25,
                  graph = FALSE)
  expect_equal(round(o$table$W, 3),
               c(0.673, 0.704, 0.708, 0.701, 0.630, 0.535, 0.334))
})

test_that("inputs are validated", {
  expect_error(orphanhood(bolivia, "age", "alive", "dead", M = 21),
               "outside the tabulated range")
  expect_error(orphanhood(bolivia, "age", "alive", "dead", M = 31),
               "outside the tabulated range")
  expect_error(orphanhood(bolivia, "age", M = 28.8), "either 'prop'")
  expect_error(orphanhood(bolivia, "nope", "alive", "dead", M = 28.8),
               "Column not found")
  rev_b <- bolivia[rev(seq_len(nrow(bolivia))), ]
  expect_error(orphanhood(rev_b, "age", "alive", "dead", M = 28.8),
               "increasing")
  bad <- bolivia; bad$alive[2] <- -1
  expect_error(orphanhood(bad, "age", "alive", "dead", M = 28.8),
               "non-negative")
})

test_that("the object prints and plots", {
  o <- orphanhood(bolivia, "age", "alive", "dead", prop10 = "s10",
                  M = 28.8, survey_year = 1975.5)
  expect_s3_class(o, "dem_orphanhood")
  expect_s3_class(plot(o), "ggplot")
  out <- capture.output(print(o))
  expect_true(any(grepl("maternal orphanhood", out)))
  expect_true(any(grepl("28.80", out)))

  o0 <- orphanhood(bolivia, "age", "alive", "dead", M = 28.8, graph = FALSE)
  expect_null(o0$plot)
  expect_error(plot(o0), "No plot available")
})
