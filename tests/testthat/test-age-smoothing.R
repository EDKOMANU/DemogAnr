# Five-year smoothing is checked three ways: one group worked by hand from
# each published formula (which catches a wrong neighbour), the distributions
# each formula provably leaves unchanged, and recovery of a known population
# from a known pattern of misreporting.

# Ten-year totals 1200, 1000, 700, 600, 300, split unevenly
toy <- data.frame(age = seq(0, 45, 5),
                  pop = c(650, 550, 520, 480, 380, 320, 290, 310, 160, 140))
sm <- function(method, ...) {
  smooth_ages(toy, "age", "pop", method = method, open_ended = FALSE,
              graph = FALSE, ...)$table$smoothed
}

test_that("each formula reproduces a group worked by hand", {
  # Karup-King-Newton, ages 10-19: 1/2 (1000) + 1/16 (1200 - 700)
  expect_equal(sm("kkn")[3:4], c(531.25, 468.75))
  # Carrier-Farrag, ages 10-19: 5P15 = 1000 / (1 + (1200 / 700)^(1/4))
  cf15 <- 1000 / (1 + (1200 / 700)^0.25)
  expect_equal(sm("carrier_farrag")[3:4], c(1000 - cf15, cf15))
  # Arriaga, middle group: 5P15 = (-1200 + 11 * 1000 + 2 * 700) / 24
  expect_equal(sm("arriaga")[3:4], c(1000 - 11200 / 24, 11200 / 24))
  # Arriaga, youngest group: 5P5 = (8 * 1200 + 5 * 1000 - 700) / 24
  expect_equal(sm("arriaga")[1:2], c(1200 - 13900 / 24, 13900 / 24))
  # Arriaga, oldest group: 5P40 = (-700 + 5 * 600 + 8 * 300) / 24
  expect_equal(sm("arriaga")[9:10], c(4700 / 24, 300 - 4700 / 24))
  # United Nations, ages 10-14: (-650 + 4 * 550 + 10 * 520 + 4 * 480 - 380) / 16
  expect_equal(sm("united_nations", preserve_total = FALSE)[3], 8290 / 16)
  # Strong: the interior ten-year groups smoothed 1-2-1, rescaled to their
  # reported total of 2300, then split with Arriaga's formula
  T2 <- c(1200 + 2 * 1000 + 700, 1000 + 2 * 700 + 600, 700 + 2 * 600 + 300) / 4
  T2 <- T2 * 2300 / sum(T2)
  s15 <- (-1200 + 11 * T2[1] + 2 * T2[2]) / 24
  expect_equal(sm("strong")[3:4], c(T2[1] - s15, s15))
})

test_that("the methods leave alone what they can't reach", {
  # Carrier-Farrag, KKN and strong cannot split the end ten-year groups
  for (m in c("carrier_farrag", "kkn", "strong")) {
    s <- sm(m)
    expect_equal(s[c(1, 2, 9, 10)], toy$pop[c(1, 2, 9, 10)], label = m)
  }
  # the United Nations formula needs two groups on each side
  expect_equal(sm("united_nations")[c(1, 2, 9, 10)], toy$pop[c(1, 2, 9, 10)])
})

# Group totals of a density over ages 0-79, closed groups
from_density <- function(f) {
  sapply(seq(0, 75, 5), function(a) integrate(f, a, a + 5)$value)
}
unchanged <- function(P, method, ...) {
  s <- smooth_ages(data.frame(age = seq(0, 75, 5), pop = P), "age", "pop",
                   method = method, open_ended = FALSE, graph = FALSE, ...)
  max(abs(s$table$smoothed - P) / P)
}

test_that("each formula leaves its own family of distributions unchanged", {
  expo  <- from_density(function(a) 1e5 * exp(-0.035 * a))
  lin   <- from_density(function(a) 9000 - 100 * a)
  quad  <- from_density(function(a) 9000 - 180 * a + 1.2 * a^2)
  cubic <- from_density(function(a) 9000 - 250 * a + 3 * a^2 - 0.012 * a^3)

  # Carrier-Farrag: geometric decline (a constant ratio between groups)
  expect_lt(unchanged(expo, "carrier_farrag"), 1e-10)
  expect_gt(unchanged(lin, "carrier_farrag"), 1e-4)
  # Karup-King-Newton: quadratics, but not cubics
  expect_lt(unchanged(quad, "kkn"), 1e-10)
  expect_gt(unchanged(cubic, "kkn"), 1e-4)
  # Arriaga: straight lines only; it fits its parabola through the midpoints
  # of the ten-year groups rather than through their areas
  expect_lt(unchanged(lin, "arriaga"), 1e-10)
  expect_gt(unchanged(quad, "arriaga"), 1e-4)
  # United Nations: cubics, with or without rescaling
  expect_lt(unchanged(cubic, "united_nations", preserve_total = FALSE), 1e-10)
  expect_lt(unchanged(cubic, "united_nations"), 1e-10)
  # strong: straight lines
  expect_lt(unchanged(lin, "strong"), 1e-10)
  expect_gt(unchanged(quad, "strong"), 1e-4)
})

test_that("the light methods keep every ten-year total, and the others the overall total", {
  ten <- function(v) colSums(matrix(v, nrow = 2))
  for (m in c("carrier_farrag", "kkn", "arriaga")) {
    expect_equal(ten(sm(m)), ten(toy$pop), label = m)
  }
  for (m in c("strong", "united_nations")) {
    expect_equal(sum(sm(m)), sum(toy$pop), label = m)
    raw <- sm(m, preserve_total = FALSE)
    expect_false(isTRUE(all.equal(sum(raw), sum(toy$pop))), label = m)
  }
})

# A known population with known misreporting
true_pop <- from_density(function(a) 1e5 * exp(-0.035 * a))
ages <- seq(0, 75, 5)
misreport <- function(from_digit, step) {
  # 15% of every group starting at an age ending in `from_digit` is reported
  # in the neighbouring group `step` places away
  h <- true_pop
  for (i in which(ages %% 10 == from_digit & seq_along(ages) + step >= 1 &
                  seq_along(ages) + step <= length(ages))) {
    h[i] <- h[i] - 0.15 * true_pop[i]
    h[i + step] <- h[i + step] + 0.15 * true_pop[i]
  }
  data.frame(age = ages, pop = h)
}
max_err <- function(s) {
  i <- 3:(length(ages) - 2)
  max(abs(s$table$smoothed[i] / true_pop[i] - 1))
}

test_that("misreporting within ten-year groups is undone", {
  # people aged 25-29 reported at 20-24, and so on: the ten-year totals stay
  # right, so splitting them again recovers the population
  heaped <- misreport(5, -1)
  expect_gt(max(abs(heaped$pop / true_pop - 1)), 0.14)
  run <- function(m) smooth_ages(heaped, "age", "pop", method = m,
                                 open_ended = FALSE, graph = FALSE)
  # exactly, for Carrier-Farrag, since the population declines geometrically
  cf <- run("carrier_farrag")
  ch <- cf$table$changed
  expect_equal(cf$table$smoothed[ch], true_pop[ch], tolerance = 1e-9)
  for (m in c("kkn", "arriaga")) expect_lt(max_err(run(m)), 0.005, label = m)
  for (m in c("carrier_farrag", "kkn", "arriaga", "strong", "united_nations")) {
    s <- run(m)
    expect_lt(s$score_smoothed, s$score_reported / 3, label = m)
  }
})

test_that("pairing from age 5 undoes misreporting across ages ending in 0", {
  # people aged 15-19 reported at 20-24: the shift crosses the boundary of
  # the ten-year groups 10-19 and 20-29, so with the usual pairing the
  # ten-year totals are wrong and the light methods cannot repair them.
  # Pairing as 15-24 (ARRSM5) puts both sides in the same ten-year group.
  heaped <- misreport(5, 1)
  run <- function(off) smooth_ages(heaped, "age", "pop", method = "arriaga",
                                   offset = off, open_ended = FALSE,
                                   graph = FALSE)
  expect_gt(max_err(run(0)), 0.03)
  expect_lt(max_err(run(5)), 0.01)
})

test_that("with offset = 5 the ten-year groups are 5-14, 15-24, ...", {
  s <- smooth_ages(toy, "age", "pop", method = "kkn", offset = 5,
                   open_ended = FALSE, graph = FALSE)$table
  # 0-4 has no partner; pairs 5-14 ... 35-44, with 45-49 left over
  expect_false(s$changed[1])
  expect_false(s$changed[10])
  T <- toy$pop[c(2, 4, 6, 8)] + toy$pop[c(3, 5, 7, 9)]
  # KKN splits the middle pairs: 15-24 is the second
  expect_equal(s$smoothed[4], T[2] / 2 + (T[1] - T[3]) / 16)
})

test_that("lower and upper confine the smoothing", {
  s <- smooth_ages(toy, "age", "pop", method = "arriaga", lower = 10,
                   upper = 39, open_ended = FALSE, graph = FALSE)$table
  expect_equal(s$changed, c(FALSE, FALSE, rep(TRUE, 6), FALSE, FALSE))
  expect_equal(s$smoothed[c(1, 2, 9, 10)], toy$pop[c(1, 2, 9, 10)])
  # outside groups still act as neighbours: 10-19 uses the middle formula
  expect_equal(s$smoothed[4], 11200 / 24)
  # a ten-year group is smoothed only if both halves are inside the range
  s2 <- smooth_ages(toy, "age", "pop", method = "arriaga", lower = 15,
                    upper = 39, open_ended = FALSE, graph = FALSE)$table
  expect_false(any(s2$changed[3:4]))

  # strong smoothing keeps the total over the ages it smooths
  st <- smooth_ages(toy, "age", "pop", method = "strong", lower = 10,
                    upper = 29, open_ended = FALSE, graph = FALSE)$table
  expect_equal(sum(st$smoothed[3:6]), sum(toy$pop[3:6]))
  expect_equal(st$smoothed[-(3:6)], toy$pop[-(3:6)])
})

ghana1960 <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
          "40-44","45-49","50-54","55-59","60-64","65+"),
  males = c(643041, 515520, 357831, 275542, 268336, 278601, 242515,
            198231, 168937, 122756, 96775, 59307, 63467, 113185),
  females = c(654258, 503070, 323460, 265534, 322576, 306329, 245883,
              179182, 145572, 95590, 81715, 48412, 54572, 100392))

test_that("on the Ghana 1960 census every method reduces the irregularity", {
  for (sex in c("males", "females")) {
    sc <- sapply(c("carrier_farrag", "kkn", "arriaga", "united_nations",
                   "strong"), function(m) {
      s <- smooth_ages(ghana1960, "age", sex, method = m, graph = FALSE)
      expect_equal(s$total_smoothed, sum(ghana1960[[sex]]))
      expect_equal(s$table$smoothed[14], ghana1960[[sex]][14])  # 65+ kept
      c(s$score_reported, s$score_smoothed)
    })
    expect_true(all(sc[2, ] < sc[1, ]), label = sex)
    # strong smoothing moves people between ten-year groups as well, and so
    # goes furthest
    expect_true(all(sc[2, "strong"] < sc[2, c("carrier_farrag", "kkn",
                                               "arriaga")]), label = sex)
  }
})

test_that("labels and numeric ages, open or closed, agree", {
  num <- ghana1960; num$age <- seq(0, 65, 5)
  a <- smooth_ages(ghana1960, "age", "females", graph = FALSE)$table
  b <- smooth_ages(num, "age", "females", graph = FALSE)$table
  expect_equal(a$smoothed, b$smoothed)
  expect_equal(a$interval[c(1, 14)], c("0-4", "65+"))
  # closing the last group lets it take part
  c <- smooth_ages(num, "age", "females", open_ended = FALSE, graph = FALSE)
  expect_true(c$table$changed[13])
})

test_that("inputs are validated", {
  expect_error(smooth_ages(toy, "age", "pop", offset = 3), "0 or 5")
  expect_error(smooth_ages(toy[1:5, ], "age", "pop", open_ended = FALSE),
               "three")
  expect_error(smooth_ages(toy[1:4, ], "age", "pop",
                           method = "united_nations", open_ended = FALSE),
               "five closed")
  expect_error(smooth_ages(toy, "age", "pop", lower = 20, upper = 24,
                           open_ended = FALSE), "No age group could be")
  expect_error(smooth_ages(toy, "age", "pop", lower = 30, upper = 10),
               "lower <= upper")
  zero <- toy; zero$pop[5:6] <- 0
  expect_error(smooth_ages(zero, "age", "pop", method = "carrier_farrag",
                           open_ended = FALSE), "positive")
  shifted <- toy; shifted$age <- shifted$age + 2
  expect_error(smooth_ages(shifted, "age", "pop"), "multiples of 5")
  expect_error(smooth_ages(toy, "age", "pop", open_ended = NA), "TRUE or FALSE")
})

test_that("a negative smoothed count is reported", {
  steep <- data.frame(age = seq(0, 45, 5),
                      pop = c(1e5, 5e4, 2e4, 5e3, 800, 100, 2000, 5, 1, 1))
  expect_warning(smooth_ages(steep, "age", "pop", method = "kkn",
                             open_ended = FALSE, graph = FALSE), "negative")
})

test_that("the five-year object prints and plots", {
  s <- smooth_ages(ghana1960, "age", "females")
  expect_s3_class(s, "dem_age_smooth")
  expect_s3_class(plot(s), "ggplot")
  out <- capture.output(print(s))
  expect_true(any(grepl("Arriaga", out)))
  expect_true(any(grepl("Left as reported: 60-64, 65\\+", out)))
  expect_true(any(grepl("Age ratio score", out)))
  s0 <- smooth_ages(ghana1960, "age", "females", graph = FALSE)
  expect_null(s0$plot)
  expect_error(plot(s0), "No plot available")
})


# --- Single years -------------------------------------------------------------

# The weights the smoother applies at each age, read off by smoothing each
# unit vector in turn (the moving average is linear): row i is age i.
weight_matrix <- function(k, ...) {
  sapply(seq_len(k), function(j) {
    e <- numeric(k); e[j] <- 1
    smooth_single_ages(data.frame(age = 0:(k - 1), pop = e), "age", "pop",
                       open_ended = FALSE, preserve_total = FALSE,
                       graph = FALSE, ...)$table$smoothed
  })
}

test_that("the POP1SM default weights and their trimmed versions are reproduced", {
  W <- 100 * weight_matrix(15)
  # interior: 4, 8, 12, 16, 20, 16, 12, 8, 4 per cent
  expect_equal(W[8, 4:12], c(4, 8, 12, 16, 20, 16, 12, 8, 4))
  # the trimmed columns of the POP1SM table, three, two and one year from
  # the end
  expect_equal(W[4, 1:7], c(6.25, 12.5, 18.75, 25, 18.75, 12.5, 6.25))
  expect_equal(round(W[3, 1:5], 2), c(11.11, 22.22, 33.33, 22.22, 11.11))
  expect_equal(W[2, 1:3], c(25, 50, 25))
  # the youngest age is left as it is, and the oldest mirrors the youngest
  expect_equal(W[1, ], c(100, rep(0, 14)))
  expect_equal(W[15:12, 15:1], W[1:4, ])
  # every age's weights sum to one
  expect_equal(rowSums(W), rep(100, 15))
})

test_that("POP1SM's alternative weights and their trimmed versions are reproduced", {
  # POP1SM states these weights on a scale that sums to 90, not 100
  W <- 90 * weight_matrix(15, weights = c(4.5, rep(9, 9), 4.5))
  expect_equal(W[8, 3:13], c(4.5, rep(9, 9), 4.5))
  expect_equal(W[5, 1:9], rep(10, 9))
  expect_equal(round(W[4, 1:7], 2), rep(12.86, 7))
  expect_equal(W[3, 1:5], rep(18, 5))
  expect_equal(W[2, 1:3], rep(30, 3))
})

test_that("a straight line passes through untouched", {
  # every window is symmetric, even where it narrows, so a linear
  # distribution is reproduced exactly at every age
  d <- data.frame(age = 0:40, pop = 5000 - 90 * (0:40))
  s <- smooth_single_ages(d, "age", "pop", open_ended = FALSE,
                          preserve_total = FALSE, graph = FALSE)
  expect_equal(s$table$smoothed, d$pop)
})

test_that("digit preference is removed from Ghana's single years", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  s <- smooth_single_ages(male, "Age", "Population", graph = FALSE)
  expect_gt(s$whipple_reported, 110)
  expect_lt(abs(s$whipple_smoothed - 100), 2)
  # Myers' blend over a range the data support in full (the series stops
  # at 80+)
  m_before <- myers(male, "Age", "Population", upper = 71, graph = FALSE)$index
  sm <- data.frame(Age = s$table$age, Population = s$table$smoothed)
  m_after <- myers(sm, "Age", "Population", upper = 71, graph = FALSE)$index
  expect_lt(m_after, m_before / 3)
  # the total is kept, and the open group 80+ is carried through
  expect_equal(s$total_smoothed, sum(male$Population))
  expect_equal(tail(s$table$smoothed, 1), tail(male$Population, 1))
  # without rescaling the total moves, but only a little
  raw <- smooth_single_ages(male, "Age", "Population", preserve_total = FALSE,
                            graph = FALSE)
  expect_lt(abs(raw$total_smoothed / raw$total_reported - 1), 0.01)
})

test_that("a wider kernel smooths harder", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  rough <- function(span) {
    v <- smooth_single_ages(male, "Age", "Population", span = span,
                            graph = FALSE)$table$smoothed
    sum(diff(v[11:70], differences = 2)^2)
  }
  expect_lt(rough(6), rough(4))
  expect_lt(rough(4), rough(2))
})

test_that("lower and upper confine the single-year smoothing", {
  d <- data.frame(age = 0:30, pop = rep(c(1200, 900, 1000, 1000, 900), 7)[1:31])
  s <- smooth_single_ages(d, "age", "pop", lower = 10, upper = 19,
                          open_ended = FALSE, graph = FALSE)$table
  expect_equal(s$changed, d$age >= 10 & d$age <= 19)
  expect_equal(s$smoothed[!s$changed], d$pop[!s$changed])
  expect_equal(sum(s$smoothed[s$changed]), sum(d$pop[s$changed]))
})

test_that("single-year inputs are validated", {
  d <- data.frame(age = 0:20, pop = 100)
  expect_error(smooth_single_ages(d[-5, ], "age", "pop"), "single years")
  expect_error(smooth_single_ages(d, "age", "pop", weights = c(1, 2)), "odd")
  expect_error(smooth_single_ages(d, "age", "pop", weights = c(1, 2, 3)),
               "symmetric")
  expect_error(smooth_single_ages(d, "age", "pop", weights = c(-1, 2, -1)),
               "non-negative")
  expect_error(smooth_single_ages(d, "age", "pop", weights = c(1, 0, 1)),
               "positive")
  expect_error(smooth_single_ages(d, "age", "pop", weights = "box"),
               "triangular")
  expect_error(smooth_single_ages(d, "age", "pop", span = 2.5), "whole number")
  expect_error(smooth_single_ages(d, "age", "pop", lower = 50), "No age")
})

test_that("the single-year object prints and plots", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  s <- smooth_single_ages(male, "Age", "Population")
  expect_s3_class(s, "dem_age_smooth_single")
  expect_s3_class(plot(s), "ggplot")
  out <- capture.output(print(s))
  expect_true(any(grepl("4.00  8.00  12.00  16.00  20.00", out)))
  expect_true(any(grepl("Whipple", out)))
  expect_lt(length(out), 12)
  s0 <- smooth_single_ages(male, "Age", "Population", graph = FALSE)
  expect_error(plot(s0), "No plot available")
})
