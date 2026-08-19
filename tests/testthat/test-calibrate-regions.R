test_that("calibrate_regions reproduces proportional allocation with no deviation", {
  data(region2000, envir = environment())
  base_total <- sum(region2000$base_pop)
  national <- data.frame(year = c(2000, 2010, 2020),
                         total = round(base_total * c(1, 1.28, 1.65)))

  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           graph = FALSE)

  expect_s3_class(cal, "dem_calibration")
  # Every raking factor is exactly 1 when nothing deviates.
  expect_equal(cal$rake$rake_factor, rep(1, nrow(cal$rake)))
  # Calibrated = national total times the base share, region by region, year by
  # year: the ordinary constant-share allocation.
  expected <- cal$national$total[match(cal$table$year, cal$national$year)] *
    cal$table$share
  expect_equal(cal$table$calibrated, expected)
  # Base year returns the base population unchanged.
  b <- cal$table[cal$table$year == 2000, ]
  expect_equal(b$calibrated, b$base_pop)
})

test_that("raked regions sum exactly to the national control total", {
  data(region2000, envir = environment())
  base_total <- sum(region2000$base_pop)
  national <- data.frame(year = seq(2000, 2025, 5),
                         total = round(base_total * exp(0.025 * seq(0, 25, 5))))

  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           deviation = 0.10, num_samples = 200, graph = FALSE)

  for (y in national$year) {
    got <- sum(cal$table$calibrated[cal$table$year == y])
    expect_equal(got, national$total[national$year == y], tolerance = 1e-8)
  }
  expect_true(all(abs(cal$check$difference) < 1e-6))
  # Deviation actually did something: the raking factors are not all 1.
  expect_false(isTRUE(all.equal(cal$rake$rake_factor,
                                rep(1, nrow(cal$rake)))))
})

test_that("the deviation fan opens with the horizon and the base year is exact", {
  data(region2000, envir = environment())
  base_total <- sum(region2000$base_pop)
  national <- data.frame(year = c(2000, 2010, 2020),
                         total = round(base_total * c(1, 1.28, 1.65)))

  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           deviation = 0.10, num_samples = 400, graph = FALSE)
  t0 <- cal$table[cal$table$year == 2000, ]
  expect_equal(t0$lower, t0$upper)                     # no spread in base year
  expect_equal(t0$calibrated, t0$base_pop)

  width <- function(y) {
    d <- cal$table[cal$table$year == y, ]
    mean((d$upper - d$lower) / d$calibrated)
  }
  expect_gt(width(2010), 0)
  expect_gt(width(2020), width(2010))
})

test_that("results are reproducible for a fixed seed and differ across seeds", {
  data(region2000, envir = environment())
  national <- data.frame(year = c(2000, 2015),
                         total = c(sum(region2000$base_pop), 30000000))
  args <- list(region2000, region_col = "Region", pop_col = "base_pop",
               national = national, year_col = "year", total_col = "total",
               deviation = 0.08, num_samples = 100, graph = FALSE)

  a <- do.call(calibrate_regions, args)
  b <- do.call(calibrate_regions, args)
  expect_equal(a$table, b$table)

  c2 <- do.call(calibrate_regions, c(args, list(random_seed = 7)))
  expect_false(isTRUE(all.equal(a$table$median, c2$table$median)))
})

test_that("rake = FALSE leaves the allocation unadjusted but reports the gap", {
  data(region2000, envir = environment())
  national <- data.frame(year = c(2000, 2020),
                         total = c(sum(region2000$base_pop), 35000000))

  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           deviation = 0.12, num_samples = 200,
                           rake = FALSE, graph = FALSE)
  expect_equal(cal$table$calibrated, cal$table$unraked)
  expect_equal(cal$rake$rake_factor, rep(1, nrow(cal$rake)))
  # Without raking the regions need not add to the control total.
  expect_true(any(abs(cal$check$difference) > 1e-6))
})

test_that("calibration by age group rakes within every age group and year", {
  data(gphc2010, envir = environment())
  nat_age <- data.frame(
    Age   = rep(gphc2010$Age, 3),
    year  = rep(c(2010, 2015, 2020), each = nrow(gphc2010)),
    total = round(rep(gphc2010$Pop, 3) *
                    rep(c(1, 1.13, 1.28), each = nrow(gphc2010))))
  north <- round(gphc2010$Pop * 0.4)
  reg_age <- data.frame(
    Region = rep(c("North", "South"), each = nrow(gphc2010)),
    Age    = rep(gphc2010$Age, 2),
    pop    = c(north, gphc2010$Pop - north))

  cal <- calibrate_regions(reg_age, region_col = "Region", pop_col = "pop",
                           age_col = "Age", national = nat_age,
                           year_col = "year", total_col = "total",
                           deviation = 0.06, num_samples = 100, graph = FALSE)

  expect_true("age" %in% names(cal$table))
  expect_true(all(abs(cal$check$difference) < 1e-6))
  expect_equal(nrow(cal$check), nrow(gphc2010) * 3)
  # And the age groups add to the all-age national total in each year.
  for (y in c(2010, 2015, 2020)) {
    got <- sum(cal$table$calibrated[cal$table$year == y])
    expect_equal(got, sum(nat_age$total[nat_age$year == y]), tolerance = 1e-6)
  }
})

test_that("per-region and user-supplied deviations are honoured", {
  data(region2000, envir = environment())
  national <- data.frame(year = c(2000, 2020),
                         total = c(sum(region2000$base_pop), 35000000))

  sds <- stats::setNames(rep(0.02, nrow(region2000)), region2000$Region)
  sds["Greater Accra"] <- 0.20
  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           deviation = sds, num_samples = 500, graph = FALSE)
  d <- cal$table[cal$table$year == 2020, ]
  rel <- (d$upper - d$lower) / d$calibrated
  expect_gt(rel[d$region == "Greater Accra"], max(rel[d$region != "Greater Accra"]))

  # A user sampler in the same form as project_population()
  cal_fn <- calibrate_regions(region2000, region_col = "Region",
                              pop_col = "base_pop", national = national,
                              year_col = "year", total_col = "total",
                              deviation = function(mean, n)
                                stats::rnorm(n, mean, mean * 0.05),
                              num_samples = 200, graph = FALSE)
  expect_true(all(abs(cal_fn$check$difference) < 1e-6))
  expect_true(all(cal_fn$table$calibrated >= 0))

  # deviation = 0 is the deterministic case, identical to deviation = NULL
  cal0 <- calibrate_regions(region2000, region_col = "Region",
                            pop_col = "base_pop", national = national,
                            year_col = "year", total_col = "total",
                            deviation = 0, graph = FALSE)
  calN <- calibrate_regions(region2000, region_col = "Region",
                            pop_col = "base_pop", national = national,
                            year_col = "year", total_col = "total",
                            graph = FALSE)
  expect_equal(cal0$table, calN$table)
})

test_that("calibrate_regions validates its inputs", {
  data(region2000, envir = environment())
  national <- data.frame(year = c(2000, 2020),
                         total = c(sum(region2000$base_pop), 35000000))
  base_args <- list(region2000, region_col = "Region", pop_col = "base_pop",
                    national = national, year_col = "year",
                    total_col = "total", graph = FALSE)

  expect_error(do.call(calibrate_regions,
                       modifyList(base_args, list(pop_col = "nope"))),
               "not found in 'data'")
  expect_error(do.call(calibrate_regions,
                       modifyList(base_args, list(total_col = "nope"))),
               "not found in 'national'")
  expect_error(do.call(calibrate_regions,
                       modifyList(base_args, list(probs = c(0.9, 0.1)))),
               "increasing quantiles")
  expect_error(do.call(calibrate_regions,
                       modifyList(base_args, list(deviation = c(0.1, 0.2)))),
               "one value per region")
  expect_error(do.call(calibrate_regions,
                       modifyList(base_args, list(deviation = -0.1))),
               "non-negative")
  # An age group present in 'data' but absent from the control totals
  d_age <- region2000
  d_age$Age <- rep(c("0-14", "15+"), length.out = nrow(d_age))
  nat_age <- data.frame(Age = "0-14", year = c(2000, 2020),
                        total = c(10000000, 14000000))
  expect_error(calibrate_regions(d_age, region_col = "Region",
                                 pop_col = "base_pop", age_col = "Age",
                                 national = nat_age, year_col = "year",
                                 total_col = "total", graph = FALSE),
               "no control total")
})

test_that("the print and plot methods work", {
  data(region2000, envir = environment())
  national <- data.frame(year = c(2000, 2020),
                         total = c(sum(region2000$base_pop), 35000000))
  cal <- calibrate_regions(region2000, region_col = "Region",
                           pop_col = "base_pop", national = national,
                           year_col = "year", total_col = "total",
                           deviation = 0.05, num_samples = 50, graph = TRUE)
  expect_s3_class(plot(cal), "ggplot")
  expect_output(print(cal), "Regional calibration")
  expect_output(print(cal), "Raking factors")

  cal_ng <- calibrate_regions(region2000, region_col = "Region",
                              pop_col = "base_pop", national = national,
                              year_col = "year", total_col = "total",
                              graph = FALSE)
  expect_error(plot(cal_ng), "No plot available")
})
