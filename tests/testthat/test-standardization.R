d <- data.frame(
  age  = c("0-14","15-44","45-64","65+"),
  m1   = c(0.002, 0.001, 0.008, 0.060),
  m2   = c(0.003, 0.0015, 0.010, 0.065),
  n1   = c(300000, 400000, 200000, 100000),
  n2   = c(200000, 350000, 250000, 200000)
)

test_that("direct standardization equals the weighted average of rates", {
  s <- standardize(d, "age", method = "direct", rate_col = "m1",
                   std_pop_col = "n2", per = 1000)
  w <- d$n2 / sum(d$n2)
  expect_equal(s$standardized_rate, sum(w * d$m1) * 1000)
})

test_that("indirect standardization returns a correct SMR", {
  dd <- d; dd$deaths <- dd$n1 * dd$m1  # observed deaths in pop 1
  s <- standardize(dd, "age", method = "indirect", pop_col = "n1",
                   deaths_col = "deaths", std_rate_col = "m2",
                   std_pop_col = "n2", per = 1000)
  expect_equal(s$SMR, sum(dd$deaths) / sum(dd$n1 * dd$m2))
  expect_true(is.finite(s$standardized_rate))
})

test_that("Kitagawa decomposition is exactly additive", {
  kd <- decompose_rates(d, "age", "m1", "m2", "n1", "n2", per = 1000)
  expect_equal(kd$rate_component + kd$composition_component, kd$total,
               tolerance = 1e-9)
  # total equals the difference in crude rates
  cdr1 <- sum(d$n1 / sum(d$n1) * d$m1) * 1000
  cdr2 <- sum(d$n2 / sum(d$n2) * d$m2) * 1000
  expect_equal(kd$total, cdr1 - cdr2, tolerance = 1e-9)
})

test_that("Arriaga LE decomposition sums to the e0 difference", {
  data(gphc2010, envir = environment())
  m <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", sex = "male")
  f <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", sex = "female")
  ad <- decompose_LE(m$lifetable, f$lifetable)
  expect_equal(sum(ad$table$contribution), ad$total, tolerance = 1e-6)
  expect_equal(ad$direct + ad$indirect, ad$total, tolerance = 1e-6)
  expect_equal(ad$total, ad$e0_2 - ad$e0_1, tolerance = 1e-9)
})
