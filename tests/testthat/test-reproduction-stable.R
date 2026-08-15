# Preston et al. (2001), Boxes 5.5 and 7.1

us1991 <- data.frame(
  age   = seq(10, 45, 5),
  fbir  = c(5816, 253979, 532712, 596823, 431694, 162005, 25531, 829),
  women = c(8620, 8371, 9419, 10325, 11125, 10344, 9496, 7188) * 1000,
  nLx   = c(494603, 493804, 492552, 491138, 489356, 486941, 483577, 478475)
)

egypt <- data.frame(
  age = seq(15, 45, 5),
  Lx  = c(4.66740, 4.63097, 4.58518, 4.53206, 4.46912, 4.39135, 4.28969),
  ma  = c(0.00567, 0.06627, 0.11204, 0.07889, 0.05075, 0.01590, 0.00610)
)

test_that("reproduction() reproduces Preston Box 5.5 (US 1991)", {
  r <- reproduction(us1991, age = "age", births = "fbir", women = "women",
                    nLx = "nLx", radix = 1e5, fraction_female = 1, graph = FALSE)
  expect_equal(round(r$GRR, 3), 1.013)   # published: 1.013
  expect_equal(round(r$NRR, 3), 0.995)   # published: 0.995
})

test_that("reproduction(): no mortality makes NRR equal GRR", {
  d <- data.frame(age = seq(15, 45, 5),
                  asfr = c(0.02, 0.12, 0.11, 0.06, 0.02, 0.005, 0.001),
                  nLx = rep(5, 7))          # Lx = n = 5 everywhere (all survive)
  r <- reproduction(d, age = "age", asfr = "asfr", nLx = "nLx", radix = 1,
                    fraction_female = 1, graph = FALSE)
  expect_equal(r$NRR, r$GRR, tolerance = 1e-12)
})

test_that("stable_population() reproduces Preston Box 7.1 (Egypt 1997)", {
  s <- stable_population(egypt, age = "age", asfr = "ma", nLx = "Lx",
                         radix = 1, fraction_female = 1, graph = FALSE)
  expect_equal(round(s$NRR, 2), 1.53)      # published: 1.53
  expect_equal(round(s$r, 5), 0.01424)     # published: 0.01424
  expect_lt(s$iterations, 8)
})

test_that("stable_population(): Lotka identity and internal consistency hold", {
  s <- stable_population(egypt, age = "age", asfr = "ma", nLx = "Lx",
                         radix = 1, fraction_female = 1, graph = FALSE)
  mid <- s$table$age + s$table$n / 2
  Lx  <- s$table$nLx                      # radix = 1 here
  matf <- s$table$asfr
  # y(r) = sum e^{-r*mid} Lx m = 1 at the solution
  expect_equal(sum(exp(-s$r * mid) * Lx * matf), 1, tolerance = 1e-6)
  # b = 1 / sum e^{-r*mid} Lx ; d = b - r ; c(x) sums to 1
  expect_equal(s$b, 1 / sum(exp(-s$r * mid) * Lx), tolerance = 1e-9)
  expect_equal(s$d, s$b - s$r, tolerance = 1e-12)
  expect_equal(sum(s$table$cx), 1, tolerance = 1e-9)
})

test_that("stable_population(): NRR = 1 gives r ~ 0 and c(x) = Lx / sum(Lx)", {
  z <- egypt
  s0 <- stable_population(z, age = "age", asfr = "ma", nLx = "Lx",
                          radix = 1, fraction_female = 1, graph = FALSE)
  z$ma <- z$ma / s0$NRR                    # rescale fertility so NRR = 1
  s <- stable_population(z, age = "age", asfr = "ma", nLx = "Lx",
                         radix = 1, fraction_female = 1, graph = FALSE)
  expect_equal(s$NRR, 1, tolerance = 1e-9)
  expect_equal(s$r, 0, tolerance = 1e-6)
  expect_equal(s$table$cx, z$Lx / sum(z$Lx), tolerance = 1e-5)
})
