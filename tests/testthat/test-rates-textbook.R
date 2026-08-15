# Published worked examples for the basic rate functions.

test_that("dem.fert() reproduces the Manual X Bangladesh 1974 schedule", {
  bd <- data.frame(
    age         = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
    women       = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
    live_births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125),
    population  = 71315944)
  f <- dem.fert(bd, type = "all", age_col = "age", population_col = "population",
                women_col = "women", births_col = "live_births", graph = FALSE)
  # published TFR of the reported schedule (Manual X, Table 12)
  expect_equal(round(f$results$TFR, 2), 4.83)
  # age-specific rates (Manual X, Table 10, col 4), per 1000
  expect_equal(round(f$results$ASFR / 1000, 4),
               c(0.1063, 0.2296, 0.2154, 0.1825, 0.1339, 0.0644, 0.0336),
               tolerance = 2e-3)
})

test_that("dem.cdr() reproduces Preston Box 1.2 (Sweden, 1988)", {
  sweden88 <- data.frame(population = 8438477, deaths = 96756)
  d <- dem.cdr(sweden88, type = "CDR", population_col = "population",
               deaths_col = "deaths", graph = FALSE)
  expect_equal(round(d$results$CDR, 2), 11.47)   # published 0.01147
})
