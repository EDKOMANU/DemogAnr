# The renamed API, and the aliases kept so that earlier code keeps working.

test_that("dem.chm() is the canonical name and dm.chm() still works", {
  data(ghmort2021, envir = environment())
  males <- subset(ghmort2021, Sex == "Male")

  new <- dem.chm(males, age_col = "Age", live_births = "Population",
                 deaths_col = "Deaths", type = "infant")
  old <- suppressWarnings(
    dm.chm(males, age_col = "Age", live_births = "Population",
           deaths_col = "Deaths", type = "infant"))

  expect_equal(as.numeric(old), as.numeric(new))
  expect_s3_class(new, "dem_chm")
  expect_equal(attr(new, "type"), "infant")
})

test_that("dm.chm() warns once, not on every call", {
  data(ghmort2021, envir = environment())
  males <- subset(ghmort2021, Sex == "Male")
  # the alias may already have warned in this session, so clear the record
  if (exists(".dep_warned", envir = asNamespace("DemogAnr"))) {
    reg <- get(".dep_warned", envir = asNamespace("DemogAnr"))
    suppressWarnings(rm(list = ls(reg), envir = reg))
  }
  call_it <- function() {
    dm.chm(males, age_col = "Age", live_births = "Population",
           deaths_col = "Deaths", type = "infant")
  }
  expect_warning(call_it(), "deprecated")
  expect_silent(call_it())
})

test_that("dem.mmr() reports the rate and the ratio under unambiguous names", {
  d <- data.frame(deaths = 120, births = 100000, women = 5e6)
  m <- dem.mmr(d, "deaths", "births", "women")

  # ratio: per 100,000 live births; rate: per 100,000 women
  expect_equal(m$results$MMRatio, 120)
  expect_equal(m$results$MMRate, 2.4)

  # the old names keep their original values, so existing code is unaffected
  expect_equal(m$results$MMR, m$results$MMRate)
  expect_equal(m$results$MMR_Ratio, m$results$MMRatio)
})

test_that("dem.mmr() prints without dumping a long data frame", {
  d <- data.frame(deaths = rep(10, 50), births = rep(1000, 50),
                  women = rep(50000, 50))
  out <- capture.output(print(dem.mmr(d, "deaths", "births", "women")))
  expect_true(any(grepl("MMRate", out)))
  expect_true(any(grepl("MMRatio", out)))
  expect_true(any(grepl("\\(50 rows\\)", out)))
  expect_lt(length(out), 25)
})

test_that("grouped_pop replaces the dataset formerly called 'data'", {
  data(grouped_pop, envir = environment())
  expect_s3_class(grouped_pop, "data.frame")
  expect_true("age_col" %in% names(grouped_pop))
  expect_equal(nrow(grouped_pop), 16)

  # the old name is still shipped and identical, for code written against it
  data(data, envir = environment())
  expect_equal(as.data.frame(get("data", envir = environment())),
               as.data.frame(grouped_pop))

  # and it drives karup_king() as documented
  single <- karup_king(grouped_pop, age_col = "age_col", pops = "2021")
  expect_equal(sum(single[["2021"]]), sum(grouped_pop[["2021"]]),
               tolerance = 0.001)
})
