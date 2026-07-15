test_that("dem.cdr computes CDR and ASDR correctly", {
  d <- data.frame(
    age = c("0-4", "5-9", "10-14"),
    population = c(10000, 12000, 11000),
    deaths = c(100, 50, 30)
  )
  res <- dem.cdr(d, type = "all", age_col = "age",
                 population_col = "population", deaths_col = "deaths")
  expect_s3_class(res, "dem_cdr")
  expect_equal(res$results$CDR, (180 / 33000) * 1000)
  expect_equal(res$results$ASDR, c(10, 50/12, 30/11))
  # 'all' expansion must also work when type is already a vector (R >= 4.2)
  res2 <- dem.cdr(d, type = c("CDR", "ASDR"), age_col = "age",
                  population_col = "population", deaths_col = "deaths")
  expect_equal(res2$results$CDR, res$results$CDR)
})

test_that("dem.fert computes CBR, GFR, ASFR and TFR correctly", {
  d <- data.frame(
    age = c("15-19", "20-24", "25-29"),
    population = c(10000, 12000, 11000),
    women = c(5000, 6000, 5500),
    live_births = c(200, 400, 300)
  )
  res <- dem.fert(d, type = "all", age_col = "age", population_col = "population",
                  women_col = "women", births_col = "live_births")
  expect_equal(res$results$CBR, (900 / 33000) * 1000)
  expect_equal(res$results$GFR, (900 / 16500) * 1000)
  asfr <- c(200/5000, 400/6000, 300/5500) * 1000
  expect_equal(res$results$ASFR, asfr)
  expect_equal(res$results$TFR, sum(asfr) * 5 / 1000)
  # custom age interval
  res1 <- dem.fert(d, type = "TFR", age_col = "age", population_col = "population",
                   women_col = "women", births_col = "live_births", age_interval = 1)
  expect_equal(res1$results$TFR, sum(asfr) / 1000)
})

test_that("dem.mmr computes rate and ratio correctly", {
  d <- data.frame(
    maternal_deaths = c(5, 10, 15),
    live_births = c(5000, 7000, 6000),
    women = c(20000, 22000, 32000)
  )
  res <- dem.mmr(d, deaths_col = "maternal_deaths", births_col = "live_births",
                 pop_women = "women")
  expect_equal(res$results$MMR, 30 / 74000 * 1e5)
  expect_equal(res$results$MMR_Ratio, 30 / 18000 * 1e5)
})

test_that("dm.chm computes mortality per 1,000 and filters ages", {
  d <- data.frame(
    age = 0:10,
    population = c(1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500),
    deaths = c(30, 20, 15, 10, 5, 2, 1, 1, 0, 0, 0)
  )
  u5 <- dm.chm(d, age_col = "age", live_births = "population",
               deaths_col = "deaths", type = "under5")
  expect_equal(as.numeric(u5), 80 / 4500 * 1000)
  inf <- dm.chm(d, age_col = "age", live_births = "population",
                deaths_col = "deaths", type = "infant")
  expect_equal(as.numeric(inf), 30)
  # neonatal counts only deaths in the first month (age_in_months < 1)
  dn <- data.frame(age = c(0, 0, 0), months = c(0, 1, 2),
                   births = c(100, 100, 100), deaths = c(5, 3, 1))
  neo <- dm.chm(dn, age_col = "age", live_births = "births", deaths_col = "deaths",
                type = "neonatal", age_in_months = "months")
  expect_equal(as.numeric(neo), 5 / 100 * 1000)
})
