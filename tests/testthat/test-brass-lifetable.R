test_that("alpha = 0, beta = 1 reproduces the standard", {
  data(standards, envir = environment())
  s <- standards[order(standards$Age), ]

  for (st in c("African", "GeneralUN")) {
    g <- brass_lifetable(alpha = 0, beta = 1, standard = st, graph = FALSE)
    # The standards table stores logits to four decimal places, so the
    # round trip is exact only to that precision, not to machine epsilon.
    expect_equal(g$lifetable$lx / 1e5, s[[st]], tolerance = 1e-4)
    expect_equal(g$brass$alpha, 0)
    expect_equal(g$brass$standard, st)
  }
})

test_that("the logit transform itself is exact", {
  data(standards, envir = environment())
  s <- standards[order(standards$Age), ]
  # invert the stored logit directly, and compare with what the function built
  direct <- 1 - 1 / (1 + exp(-2 * s$African_logit))
  direct[s$Age == 0] <- 1
  g <- brass_lifetable(0, 1, "African", graph = FALSE)
  expect_equal(g$lifetable$lx / 1e5, direct, tolerance = 1e-12)
})

test_that("the generated table is internally consistent", {
  lt <- brass_lifetable(-0.4, 1, "African", graph = FALSE)$lifetable
  k <- nrow(lt)

  expect_true(all(diff(lt$lx) < 0))                   # l(x) strictly declines
  expect_true(all(lt$nqx >= 0 & lt$nqx <= 1))
  expect_equal(lt$nqx[k], 1)                          # everyone dies eventually
  expect_equal(lt$dx[-k], -diff(lt$lx))               # dx = lx - lx(next)
  expect_equal(lt$dx[k], lt$lx[k])
  expect_equal(sum(lt$dx), 1e5)                       # deaths sum to the radix
  expect_equal(lt$Tx, rev(cumsum(rev(lt$Lx))))
  expect_equal(lt$ex, lt$Tx / lt$lx)
  expect_true(all(lt$Lx > 0))
})

test_that("alpha sets the level and beta the age pattern", {
  e0 <- function(a, b) {
    brass_lifetable(a, b, "African", graph = FALSE)$metrics$LifeExpectancyAtBirth
  }
  # alpha below zero is lighter mortality than the standard
  expect_gt(e0(-0.4, 1), e0(0, 1))
  expect_lt(e0(0.4, 1), e0(0, 1))

  q5 <- function(b) {
    t <- brass_lifetable(0, b, "African", graph = FALSE)$lifetable
    1 - t$lx[t$Age == 5] / 1e5
  }
  # beta below one raises child mortality relative to adult mortality
  expect_gt(q5(0.8), q5(1.0))
  expect_lt(q5(1.2), q5(1.0))
})

test_that("a brass_logit() fit can be expanded directly", {
  panama <- data.frame(
    women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
    ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
    cd    = c(40, 130, 312, 435, 636, 686, 689))
  est <- chm_brass(panama, "women", "ceb", "cd", model = "west")
  fit <- brass_logit(data.frame(age = est$x, qx = est$qx),
                     qx_col = "qx", age_col = "age", standard = "African")

  expect_s3_class(fit, "dem_brass_fit")
  expect_equal(fit$standard, "African")

  lt <- brass_lifetable(fit)
  # alpha, beta and the standard all come across from the fit
  expect_equal(lt$brass$alpha, fit$alpha)
  expect_equal(lt$brass$beta, fit$beta)
  expect_equal(lt$brass$standard, "African")

  # the fit covered ages 1-20 only; the table spans the whole standard
  expect_equal(nrow(fit$data), 7)
  expect_equal(max(lt$lifetable$Age), 95)

  # Panama's life expectancy in the mid-1970s was around 68 years
  expect_gt(lt$metrics$LifeExpectancyAtBirth, 60)
  expect_lt(lt$metrics$LifeExpectancyAtBirth, 75)

  # passing alpha and beta by hand gives the identical table
  by_hand <- brass_lifetable(fit$alpha, fit$beta, "African", graph = FALSE)
  expect_equal(by_hand$lifetable, lt$lifetable)
})

test_that("the open-interval assumption barely moves e0", {
  e <- function(m) {
    brass_lifetable(-0.4, 1, "African", open_mx = m,
                    graph = FALSE)$metrics$LifeExpectancyAtBirth
  }
  # survivorship at age 95 is about 0.0006 of the radix, so tripling the
  # assumed death rate there should change e0 by well under a hundredth
  expect_lt(abs(e(0.2) - e(0.6)), 0.01)
})

test_that("radix scales the table without changing the rates", {
  a <- brass_lifetable(-0.3, 1.1, "African", radix = 1e5, graph = FALSE)$lifetable
  b <- brass_lifetable(-0.3, 1.1, "African", radix = 1, graph = FALSE)$lifetable
  expect_equal(b$lx[1], 1)
  expect_equal(a$lx / 1e5, b$lx)
  expect_equal(a$nqx, b$nqx)
  expect_equal(a$ex, b$ex)
})

test_that("invalid parameters are rejected", {
  expect_error(brass_lifetable(0, -1, "African"), "positive")
  expect_error(brass_lifetable(0, 0, "African"), "positive")
  expect_error(brass_lifetable(NA, 1, "African"), "finite")
  expect_error(brass_lifetable(0, 1, "Nonexistent"), "'arg' should be one of")
})

test_that("the result behaves as a dem_lifetable", {
  g <- brass_lifetable(-0.4, 1, "African")
  expect_s3_class(g, "dem_brass_lt")
  expect_s3_class(g, "dem_lifetable")
  expect_s3_class(plot(g), "ggplot")            # inherited plot method

  out <- capture.output(print(g))
  expect_true(any(grepl("Brass relational life table", out)))
  expect_true(any(grepl("alpha", out)))
  expect_true(any(grepl("e0", out)))

  g0 <- brass_lifetable(-0.4, 1, "African", graph = FALSE)
  expect_null(g0$plot)
  expect_error(plot(g0), "No plot available")
})

test_that("a user-supplied standard is honoured", {
  data(standards, envir = environment())
  half <- standards[standards$Age <= 50, ]
  g <- brass_lifetable(0, 1, "African", standards_data = half, graph = FALSE)
  expect_equal(max(g$lifetable$Age), 50)
  expect_equal(g$lifetable$nqx[nrow(g$lifetable)], 1)
})
