# UN Manual X (1983) worked example: Panama 1976, both sexes, model West.
panama <- data.frame(
  age_group = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
  women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
  ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
  cd    = c(40, 130, 312, 435, 636, 686, 689)
)

test_that("chm_brass reproduces the Manual X Panama multipliers (Table 52)", {
  est <- chm_brass(panama, "women", "ceb", cd_col = "cd", model = "west",
                   survey_year = 1976.7)
  # multipliers k(i), Manual X Table 52 (both sexes)
  expect_equal(round(est$k, 4),
               c(1.0663, 1.0404, 0.9938, 1.0042, 1.0221, 1.0100, 1.0022),
               tolerance = 1e-3)
  # proportions dead D(i), Manual X Table 51
  expect_equal(round(est$D, 4),
               c(0.0718, 0.0494, 0.0656, 0.0715, 0.0946, 0.1077, 0.1306),
               tolerance = 1e-3)
  # child ages to which q(x) refers
  expect_equal(est$x, c(1, 2, 3, 5, 10, 15, 20))
})

test_that("chm_brass reference dates match Manual X Table 55", {
  est <- chm_brass(panama, "women", "ceb", cd_col = "cd", model = "west",
                   survey_year = 1976.7)
  # reference dates decline into the past and bracket the Table 55 values
  expect_true(all(diff(est$ref_date) < 0))
  # q(2) reference date ~1974.3, q(20) ~1961.8
  expect_equal(round(est$ref_date[2], 1), 1974.3, tolerance = 0.2)
  expect_equal(round(est$ref_date[7], 1), 1961.8, tolerance = 0.3)
})

test_that("chm_brass accepts children surviving and validates input", {
  p2 <- panama; p2$cs <- p2$ceb - p2$cd
  est <- chm_brass(p2, "women", "ceb", cs_col = "cs", model = "west")
  est_cd <- chm_brass(panama, "women", "ceb", cd_col = "cd", model = "west")
  expect_equal(est$qx, est_cd$qx)
  expect_error(chm_brass(panama[1:5, ], "women", "ceb", cd_col = "cd"), "7 rows")
  expect_error(chm_brass(panama, "women", "ceb"), "children dead|children surviving")
})

test_that("chm_brass supports all four Coale-Demeny regions", {
  for (m in c("west", "north", "south", "east")) {
    est <- chm_brass(panama, "women", "ceb", cd_col = "cd", model = m)
    expect_true(all(est$qx > 0 & est$qx < 1))
  }
})
