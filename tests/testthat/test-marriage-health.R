# SMAM (Hajnal) and Sullivan health expectancy

test_that("smam() reproduces Preston Box 4.4 (Turkish males, 1990)", {
  turkey <- data.frame(
    age    = seq(15, 50, 5),
    men    = c(3165061, 2581153, 2435765, 2096899, 1784121, 1418784,
               1111113, 980115),
    single = c(3030203, 1853222, 629077, 180767, 77134, 43412, 28627, 22527))
  s <- smam(turkey, age = "age", single = "single", total = "men")
  expect_equal(round(s$SMAM, 1), 25.0)              # published SMAM = 25.0 years
  # published pi(50) = 0.0245, averaged from proportions rounded to 3 dp;
  # computing from the raw counts gives 0.0244
  expect_lt(abs(s$prop_never - 0.0245), 5e-4)
})

test_that("smam(): step-function marriage recovers the exact age", {
  # everyone marries at exactly 25 -> single in 15-19 and 20-24 only
  m25 <- data.frame(age = seq(15, 50, 5),
                    s = c(1, 1, 0, 0, 0, 0, 0, 0))
  expect_equal(smam(m25, age = "age", prop_single = "s")$SMAM, 25)
  # everyone marries at exactly 20
  m20 <- data.frame(age = seq(15, 50, 5),
                    s = c(1, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(smam(m20, age = "age", prop_single = "s")$SMAM, 20)
})

test_that("smam(): proportions from counts match proportions directly", {
  d <- data.frame(age = seq(15, 50, 5),
                  single = c(940, 600, 290, 140, 90, 70, 60, 55),
                  total  = rep(1000, 8))
  a <- smam(d, age = "age", single = "single", total = "total")$SMAM
  d$p <- d$single / d$total
  b <- smam(d, age = "age", prop_single = "p")$SMAM
  expect_equal(a, b)
  expect_gt(a, 20); expect_lt(a, 30)     # a plausible SMAM
})

test_that("sullivan_hle(): HLE + DLE = ex at every age", {
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
                  graph = FALSE)$lifetable
  lt$disab <- pmin(0.6, 0.02 + 0.006 * (seq_len(nrow(lt)) - 1))
  s <- sullivan_hle(lt, age = "Age", nLx = "Lx", lx = "lx",
                    prevalence = "disab", graph = FALSE)
  expect_equal(s$table$HLE + s$table$DLE, s$table$ex, tolerance = 1e-9)
  expect_true(all(s$table$prop_healthy >= 0 & s$table$prop_healthy <= 1))
  expect_lt(s$HLE0, s$e0)                 # some life is spent unhealthy
})

test_that("sullivan_hle(): boundary prevalences behave", {
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
                  graph = FALSE)$lifetable
  lt$none <- 0; lt$all <- 1
  healthy <- sullivan_hle(lt, age = "Age", nLx = "Lx", lx = "lx",
                          prevalence = "none", graph = FALSE)
  sick <- sullivan_hle(lt, age = "Age", nLx = "Lx", lx = "lx",
                       prevalence = "all", graph = FALSE)
  expect_equal(healthy$HLE0, healthy$e0, tolerance = 1e-9)   # all healthy
  expect_equal(sick$HLE0, 0, tolerance = 1e-9)               # none healthy
})
