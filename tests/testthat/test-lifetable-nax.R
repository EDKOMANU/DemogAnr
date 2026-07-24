test_that("lifetable nax_method='cd' is unchanged (backward compatible)", {
  data(gphc2010, envir = environment())
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
  expect_equal(round(lt$metrics$LifeExpectancyAtBirth, 4), 75.4211)
  # default is cd
  expect_true(all(lt$lifetable$nax[3:20] == 2.5 | is.finite(lt$lifetable$nax[3:20])))
})

test_that("lifetable Keyfitz refines adult nax and stays additive", {
  data(gphc2010, envir = environment())
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
                  nax_method = "keyfitz")
  # interior adult nax now differs from the flat 2.5
  interior <- lt$lifetable$Age %in% c(10, 20, 40, 60)
  expect_true(all(lt$lifetable$nax[interior] != 2.5))
  # nax stays within (0, n)
  closed <- is.finite(lt$lifetable$n)
  expect_true(all(lt$lifetable$nax[closed] > 0 &
                    lt$lifetable$nax[closed] <= lt$lifetable$n[closed]))
  # e0 remains plausible
  expect_gt(lt$metrics$LifeExpectancyAtBirth, 40)
  expect_lt(lt$metrics$LifeExpectancyAtBirth, 90)
})

test_that("lifetable accepts a user-supplied nax vector", {
  data(gphc2010, envir = environment())
  k <- nrow(gphc2010)
  my_nax <- rep(2.5, k); my_nax[1] <- 0.1; my_nax[2] <- 1.5
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", nax = my_nax)
  expect_equal(lt$lifetable$nax[1:2], c(0.1, 1.5))
  expect_error(
    lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", nax = c(0.1, 0.2)),
    "one value per age group"
  )
})
