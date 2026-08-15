# Cohort-component projection (Preston Ch.6), cross-checked against the stable
# population model: a projection driven by fixed rates must converge to the
# stable structure and grow at the intrinsic rate.

make_schedule <- function() {
  mort <- data.frame(
    Age = seq(0, 85, 5),
    nMx = c(0.020, 0.001, 0.0008, 0.001, 0.0015, 0.002, 0.0025, 0.003,
            0.004, 0.006, 0.009, 0.014, 0.022, 0.035, 0.055, 0.090, 0.150, 0.250))
  lt <- lifetable(mort, age = "Age", nMx = "nMx", graph = FALSE)$lifetable
  data.frame(Age = lt$Age, Lx = lt$Lx,
             asfr = c(0, 0, 0, 0.05, 0.15, 0.20, 0.15, 0.10, 0.05, rep(0, 9)),
             pop  = c(rep(90, 9), 70, 55, 45, 35, 25, 16, 9, 4, 2) * 1000)
}

test_that("cohort_component() reproduces Preston Box 6.1 (Sweden 1993)", {
  sweden <- data.frame(
    Age = c(seq(0, 80, 5), 85),
    N93 = c(293395, 248369, 240012, 261346, 285209, 314388, 281290, 286923,
            304108, 324946, 247613, 211351, 215140, 221764, 223506, 183654,
            141990, 112424),
    Lx  = c(497487, 497138, 496901, 496531, 495902, 495168, 494213, 492760,
            490447, 486613, 480665, 471786, 457852, 436153, 402775, 350358,
            271512, 291707),
    Fx  = c(0, 0, 0, 0.0120, 0.0908, 0.1499, 0.1125, 0.0441, 0.0074, 0.0003,
            rep(0, 8)))
  cc <- cohort_component(sweden, age = "Age", population = "N93", nLx = "Lx",
                         asfr = "Fx", steps = 2, graph = FALSE)
  # Published totals: 4,397,428 (1993) -> 4,449,570 (1998) -> 4,478,712 (2003)
  expect_equal(round(unname(cc$totals[1])), 4397428)
  expect_equal(round(unname(cc$totals[2])), 4449570, tolerance = 2)
  expect_equal(round(unname(cc$totals[3])), 4478712, tolerance = 2)
  # youngest age group in 1998 (published 293,574)
  expect_equal(round(cc$projection$y5[1]), 293574, tolerance = 2)
})

test_that("cohort_component(): long-run growth equals the dominant eigenvalue", {
  dat <- make_schedule()
  cc <- cohort_component(dat, age = "Age", population = "pop", nLx = "Lx",
                         asfr = "asfr", steps = 200, graph = FALSE)
  tot <- cc$totals
  ratio <- tot[length(tot)] / tot[length(tot) - 1]
  expect_equal(unname(ratio), cc$lambda, tolerance = 1e-6)
  expect_true(all(cc$projection[, -1] >= 0))          # non-negative populations
})

test_that("cohort_component(): converges to the stable population", {
  dat <- make_schedule()
  cc <- cohort_component(dat, age = "Age", population = "pop", nLx = "Lx",
                         asfr = "asfr", steps = 200, graph = FALSE)
  sp <- stable_population(dat, age = "Age", asfr = "asfr", nLx = "Lx",
                          radix = 1e5, graph = FALSE)
  # growth ratio matches exp(r*n) up to the discrete/continuous gap
  expect_equal(cc$lambda, exp(sp$r * cc$n), tolerance = 5e-3)
  # long-run age structure matches the stable c(x)
  last <- cc$projection[[ncol(cc$projection)]]
  expect_equal(last / sum(last), sp$table$cx, tolerance = 1e-3)
})

test_that("cohort_component(): one-step accounting is correct", {
  dat <- make_schedule()
  cc <- cohort_component(dat, age = "Age", population = "pop", nLx = "Lx",
                         asfr = "asfr", steps = 1, graph = FALSE)
  L <- cc$matrix
  expect_equal(as.numeric(L %*% dat$pop), cc$projection$y5, tolerance = 1e-6)
  # survivors of the last closed group + open group land in the open group
  k <- nrow(dat)
  open <- dat$Lx[k] / (dat$Lx[k - 1] + dat$Lx[k])
  expect_equal(cc$projection$y5[k], open * (dat$pop[k - 1] + dat$pop[k]),
               tolerance = 1e-6)
})

test_that("cohort_component(): non-uniform age groups are rejected", {
  dat <- make_schedule()
  bad <- dat; bad$Age[2] <- 3        # break the uniform spacing
  expect_error(cohort_component(bad, age = "Age", population = "pop",
                                nLx = "Lx", asfr = "asfr", graph = FALSE),
               "uniform")
})
