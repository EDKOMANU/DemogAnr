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
  proj <- project_population(region2000[1:2, ],
    base_year = 2000, future_year = 2003,
    TFR_var = "TFR", base_pop_var = "base_pop",
    region_var = "Country", death_rate_var = "death_rate",
    net_migration_var = "net_migration", subregion_var = "Region",
    num_samples = 500
  )
  expect_true(all(c("region", "subregion", "year", "lower", "median", "mean", "upper") %in% names(proj)))
  expect_equal(nrow(proj), 2 * 4) # base year + 3 projected years, per region
  expect_true(all(proj$lower <= proj$median & proj$median <= proj$upper))
  # uncertainty grows over the horizon
  p1 <- proj[proj$subregion == proj$subregion[1], ]
  expect_gt(p1$upper[4] - p1$lower[4], p1$upper[2] - p1$lower[2])
  proj2 <- project_population(region2000[1:2, ],
    base_year = 2000, future_year = 2003,
    TFR_var = "TFR", base_pop_var = "base_pop",
    region_var = "Country", death_rate_var = "death_rate",
    net_migration_var = "net_migration", subregion_var = "Region",
    num_samples = 500
  )
  expect_equal(proj, proj2)
})
