test_that("the packaged rates cover all nine families", {
  data(model_lt, envir = environment())
  expect_equal(nlevels(model_lt$family), 9)
  expect_equal(nrow(model_lt), 9 * 2 * 45 * 18)
  expect_true(all(model_lt$nMx > 0))
  expect_false(anyNA(model_lt))
  expect_equal(range(model_lt$e0), c(54, 76))
})

test_that("a requested level returns a table close to it", {
  for (L in c(54, 60, 70, 76)) {
    m <- model_lifetable("CD_West", "male", e0 = L, graph = FALSE)
    expect_equal(m$model$e0_requested, L)
    # the tables close at 80+, which lifts the realised e0 slightly; the
    # documented drift is under 1.5 years across the whole range
    expect_gt(m$metrics$LifeExpectancyAtBirth, L)
    expect_lt(m$metrics$LifeExpectancyAtBirth - L, 1.5)
  }
})

test_that("interpolation between tabulated levels is smooth and monotone", {
  e <- vapply(seq(60, 62, 0.25), function(L) {
    model_lifetable("CD_West", "male", e0 = L,
                    graph = FALSE)$metrics$LifeExpectancyAtBirth
  }, numeric(1))
  expect_true(all(diff(e) > 0))
  # a level that is tabulated exactly must not be interpolated at all
  exact <- model_lifetable("CD_West", "male", e0 = 60.5, graph = FALSE)
  data(model_lt, envir = environment())
  want <- model_lt$nMx[model_lt$family == "CD_West" & model_lt$sex == "male" &
                       model_lt$e0 == 60.5]
  expect_equal(exact$lifetable$nMx, want)
})

test_that("matching on q1 or q5 reproduces the target exactly", {
  for (q in c(0.02, 0.05, 0.09, 0.13)) {
    m <- model_lifetable("CD_West", "male", q5 = q, graph = FALSE)
    got <- 1 - m$lifetable$lx[m$lifetable$Age == 5] / m$lifetable$lx[1]
    expect_equal(got, q, tolerance = 1e-8)
    expect_equal(m$model$matched_on, "q5")
  }
  m1 <- model_lifetable("CD_West", "female", q1 = 0.06, graph = FALSE)
  got1 <- 1 - m1$lifetable$lx[m1$lifetable$Age == 1] / m1$lifetable$lx[1]
  expect_equal(got1, 0.06, tolerance = 1e-8)
})

test_that("chm_brass() output completes into a life table", {
  panama <- data.frame(
    women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
    ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
    cd    = c(40, 130, 312, 435, 636, 686, 689))
  est <- chm_brass(panama, "women", "ceb", "cd", model = "west")
  q5 <- est$qx[est$x == 5]

  # the short family name from chm_brass() is accepted directly
  f <- model_lifetable("west", "male", q5 = q5, graph = FALSE)
  expect_equal(f$model$family, "CD_West")
  expect_equal(1 - f$lifetable$lx[f$lifetable$Age == 5] / 1e5, q5,
               tolerance = 1e-8)

  # the same child mortality implies different adult mortality by family
  e0 <- vapply(c("CD_West", "CD_North", "CD_South", "CD_East"), function(fam) {
    model_lifetable(fam, "male", q5 = q5,
                    graph = FALSE)$metrics$LifeExpectancyAtBirth
  }, numeric(1))
  expect_true(all(e0 > 55 & e0 < 75))
  expect_gt(diff(range(e0)), 1)     # the choice of family matters
})

test_that("the level and the family are validated", {
  expect_error(model_lifetable("CD_West", "male", e0 = 40), "outside the tabulated")
  expect_error(model_lifetable("CD_West", "male", e0 = 90), "outside the tabulated")
  expect_error(model_lifetable("CD_West", "male", q5 = 0.40), "outside what")
  expect_error(model_lifetable("CD_West", "male"), "exactly one")
  expect_error(model_lifetable("CD_West", "male", e0 = 60, q5 = 0.1), "exactly one")
  expect_error(model_lifetable("Nonexistent", "male", e0 = 60), "Unknown family")
  # the error names brass_lifetable as the way round the tabulated range
  expect_error(model_lifetable("CD_West", "male", q5 = 0.40), "brass_lifetable")
})

test_that("the result behaves as a dem_lifetable", {
  m <- model_lifetable("UN_South_Asian", "female", e0 = 65)
  expect_s3_class(m, "dem_model_lt")
  expect_s3_class(m, "dem_lifetable")
  expect_s3_class(plot(m), "ggplot")
  out <- capture.output(print(m))
  expect_true(any(grepl("UN_South_Asian", out)))

  lt <- m$lifetable
  expect_true(all(diff(lt$lx) < 0))
  expect_true(all(lt$nqx >= 0 & lt$nqx <= 1))
  expect_equal(lt$Tx, rev(cumsum(rev(lt$Lx))))
  expect_equal(sum(lt$dx), 1e5, tolerance = 1e-6)
})

test_that("brass_logit(complete = TRUE) attaches the life table", {
  obs <- data.frame(age = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
                    qx = c(0.10, 0.15, NA, 0.19, 0.21, NA, 0.26, 0.30,
                           NA, 0.38, 0.44))
  plain <- brass_logit(obs, "qx", "age", standard = "African")
  full  <- brass_logit(obs, "qx", "age", standard = "African", complete = TRUE)

  expect_null(plain$lifetable)
  expect_s3_class(full$lifetable, "dem_brass_lt")
  expect_equal(full$coefficients, plain$coefficients)

  # identical to doing it in two steps
  two_step <- brass_lifetable(plain, graph = FALSE)
  expect_equal(full$lifetable$lifetable, two_step$lifetable)
})
