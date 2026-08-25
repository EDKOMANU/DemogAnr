# The death distribution methods are validated against a stable population
# built by fine numerical integration, so the true completeness is known
# exactly rather than assumed.

# A stable population with a Siler-type force of mortality, integrated at
# 0.01-year steps and then grouped, so that discretisation of the age groups
# is not itself the thing under test.
synth_pop <- function(c_true = 1, r = 0.025, t = 10, coverage2 = 1) {
  h <- 0.01
  a <- seq(0, 110 - h, by = h)
  mu <- 0.006 * exp(-1.1 * a) + 0.0004 + 0.00003 * exp(0.093 * a)
  l <- exp(-cumsum(mu) * h)
  Nd <- 1e5 * exp(-r * a) * l
  Dd <- Nd * mu

  grp <- seq(0, 85, 5); k <- length(grp)
  agg <- function(dens) vapply(seq_len(k), function(i) {
    lo <- grp[i]; hi <- if (i < k) grp[i + 1] else Inf
    sum(dens[a >= lo & a < hi]) * h
  }, numeric(1))

  Nbar <- agg(Nd); Dann <- agg(Dd)
  data.frame(age = grp,
             n1 = Nbar * exp(-r * t / 2),
             n2 = Nbar * exp(r * t / 2) * coverage2,
             dth = Dann * t * c_true)
}

test_that("ggb() recovers a known registration completeness", {
  for (ct in c(1.0, 0.9, 0.7, 0.5, 0.3)) {
    g <- ggb(synth_pop(ct), "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)
    expect_equal(g$completeness, ct, tolerance = 0.01)
    expect_gt(g$r_squared, 0.999)
    # censuses cover equally well here, so the coverage ratio is 1
    expect_equal(g$coverage_ratio, 1, tolerance = 0.01)
  }
})

test_that("seg() recovers a known registration completeness", {
  for (ct in c(1.0, 0.9, 0.7, 0.5, 0.3)) {
    s <- seg(synth_pop(ct), "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)
    expect_equal(s$completeness, ct, tolerance = 0.02)
    expect_lt(s$completeness_sd, 0.05)
  }
})

test_that("seg() solves the open interval jointly with completeness", {
  # The 1/M rule uses observed deaths, which are themselves deflated by the
  # completeness, so the mean age must not drift with c.
  ages <- vapply(c(1.0, 0.7, 0.4), function(ct) {
    seg(synth_pop(ct), "age", "n1", "n2", "dth", 2000, 2010,
        graph = FALSE)$open_mean_age
  }, numeric(1))
  expect_equal(diff(range(ages)), 0, tolerance = 0.01)
  expect_gt(ages[1], 85)      # above the open age
  expect_lt(ages[1], 100)

  # taking the nominal midpoint instead biases the estimate down
  naive <- seg(synth_pop(1), "age", "n1", "n2", "dth", 2000, 2010,
               open_mean_age = 87.5, graph = FALSE)$completeness
  expect_lt(naive, 0.98)
})

test_that("ggb_seg() corrects for unequal census coverage", {
  # the second census under-enumerates by 5 per cent
  d <- synth_pop(c_true = 0.7, coverage2 = 0.95)
  g <- ggb(d, "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)
  h <- ggb_seg(d, "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)

  expect_equal(g$coverage_ratio, 1 / 0.95, tolerance = 0.02)
  expect_equal(h$completeness, 0.7, tolerance = 0.03)
  expect_s3_class(h$ggb, "dem_ggb")
  expect_s3_class(h$seg, "dem_seg")
})

test_that("annual and intercensal death counts agree", {
  d <- synth_pop(0.8)
  a <- ggb(d, "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)
  d2 <- d; d2$dth <- d$dth / 10
  b <- ggb(d2, "age", "n1", "n2", "dth", 2000, 2010, deaths_annual = TRUE,
           graph = FALSE)
  expect_equal(a$completeness, b$completeness)
})

test_that("inputs are validated", {
  d <- synth_pop(0.8)
  expect_error(ggb(d, "nope", "n1", "n2", "dth", 2000, 2010), "Column not found")
  expect_error(ggb(d, "age", "n1", "n2", "dth", 2010, 2000), "later than")
  expect_error(ggb(d, "age", "n1", "n2", "dth", 2000, 2010,
                   age_range = c(60, 62)), "Fewer than three")
  bad <- d; bad$n1[3] <- -1
  expect_error(ggb(bad, "age", "n1", "n2", "dth", 2000, 2010), "must be positive")
  rev_d <- d[rev(seq_len(nrow(d))), ]
  expect_error(ggb(rev_d, "age", "n1", "n2", "dth", 2000, 2010), "increasing")
})

test_that("the objects print and plot", {
  d <- synth_pop(0.7)
  g <- ggb(d, "age", "n1", "n2", "dth", 2000, 2010)
  s <- seg(d, "age", "n1", "n2", "dth", 2000, 2010)
  expect_s3_class(plot(g), "ggplot")
  expect_s3_class(plot(s), "ggplot")
  expect_true(any(grepl("Completeness", capture.output(print(g)))))
  expect_true(any(grepl("Completeness", capture.output(print(s)))))
  expect_true(any(grepl("growth balance",
                        capture.output(print(ggb_seg(d, "age", "n1", "n2",
                                                     "dth", 2000, 2010,
                                                     graph = FALSE))))))
  g0 <- ggb(d, "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)
  expect_error(plot(g0), "No plot available")
})
