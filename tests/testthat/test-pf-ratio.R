# United Nations (1983) Manual X, Chapter II: Brass P/F ratio worked example
# Bangladesh, 1974. Manual X rounds f(i) to 4 dp in its intermediate tables, so
# the published F/PF/K values are matched here to that rounding (tol ~2e-3).

bd <- data.frame(
  age    = seq(15, 45, 5),
  women  = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
  ceb    = c(1160919, 4901382, 9085852, 9910256, 10384001, 9164329, 6905673),
  births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125)
)

test_that("pf_ratio() reproduces Manual X Bangladesh 1974 (parities & rates)", {
  pf <- pf_ratio(bd, age = "age", women = "women", ceb = "ceb",
                 births = "births", graph = FALSE)
  # reported average parities (Table 10, col 3)
  expect_equal(pf$table$P,
               c(0.385, 1.847, 3.485, 4.917, 5.861, 6.194, 6.084),
               tolerance = 2e-3)
  # period fertility rates (Table 10, col 4)
  expect_equal(pf$table$f,
               c(0.1063, 0.2296, 0.2154, 0.1825, 0.1339, 0.0644, 0.0336),
               tolerance = 2e-3)
})

test_that("pf_ratio() reproduces the parity equivalents F and P/F ratios", {
  pf <- pf_ratio(bd, age = "age", women = "women", ceb = "ceb",
                 births = "births", graph = FALSE)
  # estimated parity equivalents F(i) (Table 11, col 4)
  expect_equal(pf$table$F,
               c(0.237, 1.209, 2.338, 3.323, 4.094, 4.503, 4.789),
               tolerance = 2e-3)
  # P/F ratios (Table 11, col 5)
  expect_equal(pf$table$PF,
               c(1.624, 1.528, 1.491, 1.480, 1.432, 1.376, 1.270),
               tolerance = 2e-3)
})

test_that("pf_ratio() reproduces the adjustment factor and adjusted TFR", {
  pf <- pf_ratio(bd, age = "age", women = "women", ceb = "ceb",
                 births = "births", graph = FALSE)
  expect_equal(pf$K, 1.500, tolerance = 2e-3)        # (1.528+1.491+1.480)/3
  expect_equal(pf$TFR, 4.83, tolerance = 2e-3)       # reported TFR
  expect_equal(pf$TFR_adjusted, 7.24, tolerance = 3e-3) # adjusted TFR
})

test_that("pf_ratio() validates its inputs", {
  bad <- bd[1:5, ]
  expect_error(pf_ratio(bad, age = "age", women = "women", ceb = "ceb",
                        births = "births", graph = FALSE), "seven age groups")
})
