test_that("math_project linear and exponential match closed forms", {
  d <- data.frame(
    Region = c("A", "A"), District = c("A1", "A2"),
    p0 = c(1000, 2000), p1 = c(1100, 2400)
  )
  out <- math_project(d, sub_group_col = "District", parent_group_col = "Region",
                      pop_cols = c("p0", "p1"), time_vals = c(2010, 2020),
                      target_times = 2030, method = "linear")
  expect_equal(out$proj_2030, c(1100 + 10 * 10, 2400 + 40 * 10))
  out2 <- math_project(d, sub_group_col = "District", parent_group_col = "Region",
                       pop_cols = c("p0", "p1"), time_vals = c(2010, 2020),
                       target_times = 2030, method = "exponential")
  expect_equal(out2$proj_2030, c(1100 * (1100/1000), 2400 * (2400/2000)))
})

test_that("math_project share methods respect parent control totals", {
  d <- data.frame(
    Region = "A", District = c("A1", "A2"),
    p0 = c(1000, 2000), p1 = c(1200, 2400), tot2030 = c(4000, 4000)
  )
  out <- math_project(d, sub_group_col = "District", parent_group_col = "Region",
                      pop_cols = c("p0", "p1"), time_vals = c(2010, 2020),
                      target_times = 2030, method = "constant_share",
                      parent_totals_cols = "tot2030")
  expect_equal(sum(out$proj_2030), 4000)
  expect_equal(out$proj_2030, c(1200/3600, 2400/3600) * 4000)
})

test_that("project_population returns quantile summaries and is reproducible", {
  data(region2000, envir = environment())
  args <- list(base_year = 2000, future_year = 2003,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    num_samples = 500, graph = FALSE)
  proj  <- do.call(project_population, c(list(region2000[1:2, ]), args))
  expect_true(all(c("region", "subregion", "year", "lower", "median", "mean", "upper") %in% names(proj)))
  expect_equal(nrow(proj), 2 * 4) # base year + 3 projected years, per region
  expect_true(all(proj$lower <= proj$median & proj$median <= proj$upper))
  # uncertainty grows over the horizon (structural parameter uncertainty)
  p1 <- proj[proj$subregion == proj$subregion[1], ]
  expect_gt(p1$upper[4] - p1$lower[4], p1$upper[2] - p1$lower[2])
  # same seed -> identical result
  proj2 <- do.call(project_population, c(list(region2000[1:2, ]), args))
  expect_equal(proj, proj2)
})

test_that("project_population follows the balancing equation g = b - d + m", {
  # With cv = 0 there is no uncertainty, so the median grows at exactly exp(g).
  one <- data.frame(Country = "X", Region = "R", base_pop = 1000,
                    cbr = 0.030, cdr = 0.010, nmr = 0.005)
  pr <- project_population(one, base_year = 2000, future_year = 2005,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    cv = 0, num_samples = 50, graph = FALSE)
  g <- 0.030 - 0.010 + 0.005
  expect_equal(pr$median[pr$year == 2005], 1000 * exp(5 * g), tolerance = 1e-8)
  # with cv = 0 the band collapses onto the median
  expect_equal(pr$lower, pr$upper, tolerance = 1e-8)
})

test_that("project_population accepts user-defined component distributions", {
  data(region2000, envir = environment())
  common <- list(base_year = 2000, future_year = 2010,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    num_samples = 800, graph = FALSE)
  tight <- do.call(project_population, c(list(region2000[1:3, ]),
    c(common, list(death_dist = function(mean, n) rlnorm(n, log(mean), 0.05)))))
  wide  <- do.call(project_population, c(list(region2000[1:3, ]),
    c(common, list(death_dist = function(mean, n) rlnorm(n, log(mean), 0.40)))))
  w_t <- tight$upper[tight$year == 2010] - tight$lower[tight$year == 2010]
  w_w <- wide$upper[wide$year == 2010]  - wide$lower[wide$year == 2010]
  # a wider mortality distribution yields a wider projection band
  expect_true(all(w_w > w_t))
})

test_that("project_population(): a target rate is reached at the horizon", {
  data(region2000, envir = environment())
  one <- region2000[1, , drop = FALSE]
  args <- list(one, base_year = 2000, future_year = 2030,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    cv = 0, num_samples = 20, graph = FALSE)
  tr <- do.call(project_population, c(args, list(birth_target = 0.020)))
  # growth in the final year must equal b_target - d + m
  g_last <- log(tr$median[tr$year == 2030] / tr$median[tr$year == 2029])
  expect_equal(g_last, 0.020 - one$cdr + one$nmr, tolerance = 1e-6)
})

test_that("project_population(): declining fertility lowers the projection", {
  data(region2000, envir = environment())
  args <- list(region2000[1:3, ], base_year = 2000, future_year = 2030,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    cv = 0, num_samples = 20, graph = FALSE)
  flat  <- do.call(project_population, args)
  trend <- do.call(project_population, c(args, list(birth_target = 0.020)))
  expect_true(all(trend$median[trend$year == 2030] < flat$median[flat$year == 2030]))
})

test_that("project_population(): no target reproduces the constant-rate model", {
  data(region2000, envir = environment())
  args <- list(region2000[1:2, ], base_year = 2000, future_year = 2010,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    num_samples = 200, graph = FALSE)
  a <- do.call(project_population, args)
  b <- do.call(project_population, c(args, list(birth_target = NULL,
                                                death_target = NULL,
                                                migration_target = NULL)))
  expect_equal(a, b)
})

test_that("project_population(): targets may be given per region by column", {
  data(region2000, envir = environment())
  d <- region2000[1:3, ]
  d$cbr_2030 <- c(0.020, 0.022, 0.018)
  args <- list(d, base_year = 2000, future_year = 2030,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    cv = 0, num_samples = 20, graph = FALSE)
  bycol <- do.call(project_population, c(args, list(birth_target = "cbr_2030")))
  for (i in seq_len(3)) {
    r <- bycol[bycol$subregion == d$Region[i], ]
    g <- log(r$median[r$year == 2030] / r$median[r$year == 2029])
    expect_equal(g, d$cbr_2030[i] - d$cdr[i] + d$nmr[i], tolerance = 1e-6)
  }
})

test_that("project_population(): target_cv widens the band at the horizon", {
  data(region2000, envir = environment())
  args <- list(region2000[1, , drop = FALSE], base_year = 2000, future_year = 2030,
    region_var = "Country", subregion_var = "Region", base_pop_var = "base_pop",
    birth_rate_var = "cbr", death_rate_var = "cdr", net_migration_var = "nmr",
    birth_target = 0.024, cv = 0.05, num_samples = 800, graph = FALSE)
  tight <- do.call(project_population, c(args, list(target_cv = 0.05)))
  wide  <- do.call(project_population, c(args, list(target_cv = 0.30)))
  w_t <- tight$upper[tight$year == 2030] - tight$lower[tight$year == 2030]
  w_w <- wide$upper[wide$year == 2030]  - wide$lower[wide$year == 2030]
  expect_gt(w_w, w_t)
})
