# Preston et al. (2001), Chapter 4: multiple-decrement & cause-deleted tables
# Built on the real 2010 Ghana census life table (gphc2010), with the observed
# all-cause deaths split into a rising "circulatory" share and the rest.

make_causes <- function() {
  data(gphc2010, envir = environment())
  g <- gphc2010
  g$circ  <- round(g$Deaths * pmin(0.55, 0.04 + 0.012 * (seq_len(nrow(g)) - 1)))
  g$other <- g$Deaths - g$circ
  g
}

test_that("multiple_decrement() reproduces Preston Box 4.1 (US females 1991)", {
  us91 <- data.frame(
    Age  = c(0, 1, seq(5, 85, 5)),
    Dall = c(15758, 3169, 1634, 1573, 3955, 4948, 6491, 9428, 12027, 15543,
             19264, 25384, 37211, 59431, 88087, 114693, 143554, 164986, 320578),
    Dneo = c(63, 275, 268, 217, 318, 467, 856, 1924, 3532, 5958,
             8434, 11673, 17078, 25263, 33534, 36695, 36571, 30220, 32739),
    lx   = c(100000, 99217, 99050, 98959, 98870, 98637, 98379, 98070, 97653,
             97083, 96289, 95008, 93018, 89882, 85249, 78711, 69618, 57486, 41756),
    nqx  = c(0.00783, 0.00168, 0.00092, 0.00090, 0.00236, 0.00262, 0.00314,
             0.00425, 0.00584, 0.00818, 0.01330, 0.02095, 0.03371, 0.05155,
             0.07669, 0.11552, 0.17427, 0.27363, 1.00000))
  us91$Dother <- us91$Dall - us91$Dneo
  m <- multiple_decrement(us91, age = "Age", causes = c("Dneo", "Dother"),
                          deaths = "Dall", lx = "lx", nqx = "nqx", graph = FALSE)
  # Published: 21,205 of 100,000 newborns eventually die of neoplasms (0.212)
  expect_equal(round(sum(m$table$d_Dneo)), 21205, tolerance = 1)
  expect_equal(round(unname(m$cause_distribution["Dneo"]), 3), 0.212)
  # the published lx/nqx are rounded, so dx sums to the radix only to ~1e-5
  expect_equal(sum(m$cause_distribution), 1, tolerance = 1e-4)
})

test_that("cause_deleted_lt() reproduces Preston Box 4.2 (US females 1991)", {
  us91 <- data.frame(
    Age  = c(0, 1, seq(5, 85, 5)),
    Dall = c(15758, 3169, 1634, 1573, 3955, 4948, 6491, 9428, 12027, 15543,
             19264, 25384, 37211, 59431, 88087, 114693, 143554, 164986, 320578),
    Dneo = c(63, 275, 268, 217, 318, 467, 856, 1924, 3532, 5958,
             8434, 11673, 17078, 25263, 33534, 36695, 36571, 30220, 32739),
    lx   = c(100000, 99217, 99050, 98959, 98870, 98637, 98379, 98070, 97653,
             97083, 96289, 95008, 93018, 89882, 85249, 78711, 69618, 57486, 41756),
    nqx  = c(0.00783, 0.00168, 0.00092, 0.00090, 0.00236, 0.00262, 0.00314,
             0.00425, 0.00584, 0.00818, 0.01330, 0.02095, 0.03371, 0.05155,
             0.07669, 0.11552, 0.17427, 0.27363, 1.00000),
    nax  = c(0.152, 1.605, 2.275, 2.843, 2.657, 2.547, 2.550, 2.616, 2.677,
             2.685, 2.681, 2.655, 2.647, 2.646, 2.631, 2.628, 2.618, 2.570, 6.539))
  d <- cause_deleted_lt(us91, age = "Age", cause = "Dneo", deaths = "Dall",
                        lx = "lx", nqx = "nqx", nax_all = "nax", graph = FALSE)
  expect_equal(round(d$e0, 2), 78.92)              # published all-cause e0
  expect_equal(d$e0_deleted, 82.46, tolerance = 0.05)  # published 82.46
  expect_equal(d$gain, 3.54, tolerance = 0.05)         # published gain
  # the alternative assumption lands in the same place
  d2 <- cause_deleted_lt(us91, age = "Age", cause = "Dneo", deaths = "Dall",
                         lx = "lx", nqx = "nqx", nax_all = "nax",
                         method = "rates", graph = FALSE)
  expect_equal(d2$e0_deleted, d$e0_deleted, tolerance = 0.05)
})

test_that("multiple_decrement(): cohort cause distribution sums to 1", {
  g <- make_causes()
  m <- multiple_decrement(g, age = "Age", causes = c("circ", "other"),
                          pop = "Pop", graph = FALSE)
  expect_equal(sum(m$cause_distribution), 1, tolerance = 1e-9)
  expect_named(m$cause_distribution, c("circ", "other"))
  expect_gt(m$cause_distribution["other"], m$cause_distribution["circ"])
  # base e0 is the real 2010 Ghana value (~75.4), not an artefact
  expect_gt(m$e0, 70); expect_lt(m$e0, 80)
})

test_that("multiple_decrement(): cause-specific decrements add to all-cause dx", {
  g <- make_causes()
  m <- multiple_decrement(g, age = "Age", causes = c("circ", "other"),
                          pop = "Pop", graph = FALSE)
  expect_equal(m$table$d_circ + m$table$d_other, m$table$dx, tolerance = 0.11)
})

test_that("multiple_decrement() base e0 matches a plain life table", {
  g <- make_causes()
  m <- multiple_decrement(g, age = "Age", causes = c("circ", "other"),
                          pop = "Pop", graph = FALSE)
  lt <- lifetable(g, age = "Age", pop = "Pop", Dx = "Deaths", graph = FALSE)
  expect_equal(m$e0, lt$lifetable$ex[1], tolerance = 1e-9)
})

test_that("cause_deleted_lt(): deleting a cause raises e0 by a plausible amount", {
  g <- make_causes()
  del <- cause_deleted_lt(g, age = "Age", cause = "circ",
                          causes = c("circ", "other"), pop = "Pop", graph = FALSE)
  expect_gt(del$gain, 0)
  expect_lt(del$gain, 15)                              # a few years, not decades
  expect_equal(del$e0_deleted, del$e0 + del$gain, tolerance = 1e-9)
  # deleting the larger cause ("other") gains more than the smaller ("circ")
  del2 <- cause_deleted_lt(g, age = "Age", cause = "other",
                           causes = c("circ", "other"), pop = "Pop", graph = FALSE)
  expect_gt(del2$gain, del$gain)
})

test_that("cause_deleted_lt(): deleting a cause with no deaths changes nothing", {
  g <- make_causes(); g$none <- 0
  del <- cause_deleted_lt(g, age = "Age", cause = "none",
                          deaths = "Deaths", pop = "Pop", graph = FALSE)
  expect_equal(del$gain, 0, tolerance = 1e-9)
  expect_equal(del$table$ex_deleted, del$table$ex, tolerance = 1e-9)
})
