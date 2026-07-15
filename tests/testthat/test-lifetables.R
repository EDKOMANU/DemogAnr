test_that("lifetable_nqx respects interval widths and gives plausible e0", {
  d <- data.frame(
    age = c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85),
    nqx = c(0.05, 0.01, 0.005, 0.002, 0.003, 0.004, 0.005, 0.007, 0.010, 0.015,
            0.020, 0.030, 0.050, 0.070, 0.100, 0.150, 0.200, 0.300, 1.000)
  )
  lt <- lifetable_nqx(d, age = "age", nqx = "nqx")
  # lx is a non-increasing survival curve starting at the radix
  expect_equal(lt$lx[1], 100000)
  expect_true(all(diff(lt$lx) <= 0))
  # deaths must sum to the radix
  expect_equal(sum(lt$dx), 100000)
  # closed-interval person-years must account for interval width n
  expect_equal(lt$Lx[3], 4 * 0 + (lt$Age[4] - lt$Age[3]) * (lt$lx[3] - lt$dx[3] / 2))
  # e0 for this schedule must be a plausible human life expectancy
  expect_gt(lt$ex[1], 55)
  expect_lt(lt$ex[1], 80)
  # ex at 80 should be small but positive
  expect_gt(lt$ex[18], 0)
})

test_that("lifetable uses correct interval widths for q0 and 4q1", {
  data(gphc2010, envir = environment())
  res <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
  lt <- res$lifetable
  m0 <- gphc2010$Deaths[1] / gphc2010$Pop[1]
  a0 <- 0.045 + 2.684 * m0
  # q0 must come from the 1-year conversion, not the 5-year one
  expect_equal(lt$nqx[1], (1 * m0) / (1 + (1 - a0) * m0), tolerance = 1e-10)
  expect_lt(lt$nqx[1], 5 * m0) # sanity: far below the old inflated value
  # last group is open-ended
  expect_equal(lt$nqx[nrow(lt)], 1)
  expect_equal(sum(lt$dx), 100000, tolerance = 1e-6)
  # e0 plausible for Ghana 2010 census data
  expect_gt(res$metrics$LifeExpectancyAtBirth, 40)
  expect_lt(res$metrics$LifeExpectancyAtBirth, 90)
  # female nax constants also work
  resf <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", sex = "female")
  expect_false(identical(res$lifetable$nqx[1], resf$lifetable$nqx[1]))
})
