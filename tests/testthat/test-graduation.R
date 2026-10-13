# The multiplier tables are validated by the identities that define them,
# which a wrong digit anywhere would break: each row of a panel sums to 1/5
# (every group contributes one fifth of itself to each of its single years in
# a flat distribution), the single years of a group add back to that group,
# and the whole population is handed out exactly once.

methods <- c("sprague", "beers_ordinary", "beers_modified", "karup_king")

# Group totals of a function of single years of age, and the five-year frame
grouped <- function(f, start = 0, k = 12) {
  single <- f(start + seq_len(5 * k) - 1)
  list(single = single,
       data = data.frame(age = start + 5 * (seq_len(k) - 1),
                         pop = colSums(matrix(single, nrow = 5))))
}

test_that("every panel row sums to one fifth", {
  for (m in names(.grad_panels)) {
    p <- .grad_panels[[m]]
    expect_equal(rowSums(p$head), rep(0.2, 10), tolerance = 1e-12)
    expect_equal(rowSums(p$mid), rep(0.2, 5), tolerance = 1e-12)
  }
  for (nm in c("first_coef", "middle_coef", "last_coef")) {
    expect_equal(rowSums(as.matrix(getExportedValue("DemogAnr", nm))),
                 rep(0.2, 5), tolerance = 1e-12)
  }
})

test_that("the whole population is handed out exactly once", {
  # every column of the full matrix sums to one, at the ends as well as in
  # the middle, and for every number of groups
  for (m in methods) for (k in c(5, 6, 9, 16)) {
    expect_equal(colSums(.grad_matrix(m, k)), rep(1, k), tolerance = 1e-12,
                 label = paste(m, "with", k, "groups"))
  }
  expect_equal(colSums(.grad_matrix("karup_king", 3)), rep(1, 3),
               tolerance = 1e-12)
})

test_that("each group is reproduced, except by Beers' modified formula", {
  for (m in methods) {
    M <- .grad_matrix(m, 10)
    blocks <- t(sapply(1:10, function(g) colSums(M[(5 * g - 4):(5 * g), ])))
    if (m == "beers_modified") {
      expect_false(isTRUE(all.equal(blocks, diag(10))))
    } else {
      expect_equal(blocks, diag(10), tolerance = 1e-12, label = m)
    }
  }
})

test_that("Sprague and Beers ordinary reproduce a cubic exactly, ends included", {
  cubic <- function(a) 2000 - 30 * a + 0.6 * a^2 - 0.004 * a^3
  for (m in c("sprague", "beers_ordinary")) {
    g <- grouped(cubic)
    out <- graduate(g$data, "age", "pop", method = m, open_ended = FALSE,
                    graph = FALSE)
    expect_equal(out$table$population, g$single, tolerance = 1e-10, label = m)
    # but not a quartic: the formulas are exact to third degree and no further
    q <- grouped(function(a) cubic(a) + 1e-4 * a^4)
    outq <- graduate(q$data, "age", "pop", method = m, open_ended = FALSE,
                     graph = FALSE)
    expect_gt(max(abs(outq$table$population - q$single)), 0.01)
  }
})

test_that("Beers modified is exact to quadratics and Karup-King to straight lines", {
  quad <- function(a) 1500 - 12 * a + 0.1 * a^2
  g <- grouped(quad)
  bm <- graduate(g$data, "age", "pop", method = "beers_modified",
                 open_ended = FALSE, graph = FALSE)
  expect_equal(bm$table$population, g$single, tolerance = 1e-10)

  lin <- grouped(function(a) 1500 - 12 * a)
  kk <- graduate(lin$data, "age", "pop", method = "karup_king",
                 open_ended = FALSE, graph = FALSE)
  expect_equal(kk$table$population, lin$single, tolerance = 1e-10)
  kkq <- graduate(g$data, "age", "pop", method = "karup_king",
                  open_ended = FALSE, graph = FALSE)
  expect_gt(max(abs(kkq$table$population - g$single)), 0.01)
})

test_that("the results do not depend on where the age range starts", {
  cubic <- function(a) 900 + 4 * a - 0.05 * a^2 + 0.0002 * a^3
  g <- grouped(cubic, start = 15, k = 8)
  out <- graduate(g$data, "age", "pop", open_ended = FALSE, graph = FALSE)
  expect_equal(out$table$age, 15:54)
  expect_equal(out$table$population, g$single, tolerance = 1e-10)
})

test_that("graduate(method = 'karup_king') is karup_king() without the rounding", {
  data(grouped_pop, envir = environment())
  old <- karup_king(grouped_pop, age_col = "age_col", pops = "2021")
  new <- graduate(grouped_pop, "age_col", "2021", method = "karup_king",
                  open_ended = FALSE, graph = FALSE)$table
  expect_equal(new$age, as.numeric(old$age))
  expect_lte(max(abs(new$population - old[["2021"]])), 0.5)
})

test_that("Beers ordinary is smoother than Sprague, which is what it is for", {
  data(grouped_pop, envir = environment())
  r5 <- function(m) {
    v <- graduate(grouped_pop, "age_col", "2021", method = m,
                  open_ended = FALSE, graph = FALSE)$table$population
    sum(diff(v, differences = 5)^2)
  }
  expect_lt(r5("beers_ordinary"), 0.5 * r5("sprague"))
})

test_that("Beers modified keeps the total but redistributes between groups", {
  data(grouped_pop, envir = environment())
  g <- graduate(grouped_pop, "age_col", "2021", method = "beers_modified",
                open_ended = FALSE, graph = FALSE)
  expect_equal(sum(g$table$population), sum(grouped_pop[["2021"]]))
  expect_gt(max(abs(g$groups$difference)), 1)
  expect_equal(sum(g$groups$difference), 0, tolerance = 1e-8)
})

test_that("an open-ended group is carried over and kept out of the arithmetic", {
  g <- grouped(function(a) 1000 - 10 * a, k = 10)
  d <- rbind(g$data, data.frame(age = 50, pop = 12345))
  out <- graduate(d, "age", "pop", graph = FALSE)
  closed <- graduate(g$data, "age", "pop", open_ended = FALSE, graph = FALSE)
  expect_equal(nrow(out$table), 51)
  expect_equal(out$table$population[51], 12345)
  expect_equal(out$table$age[51], 50)
  expect_equal(out$table$population[1:50], closed$table$population)
  expect_equal(nrow(out$groups), 10)
})

test_that("age labels and numeric lower bounds give the same result", {
  g <- grouped(function(a) 1000 - 10 * a, k = 8)
  lab <- g$data
  lab$age <- c(paste0(seq(0, 30, 5), "-", seq(4, 34, 5)), "35+")
  a <- graduate(g$data, "age", "pop", graph = FALSE)
  b <- graduate(lab, "age", "pop", graph = FALSE)
  expect_equal(a$table, b$table)
  expect_equal(b$groups$interval[1], "0-4")
})

test_that("a negative single year is reported, not hidden", {
  steep <- data.frame(age = seq(0, 45, 5),
                      pop = c(1e5, 6e4, 3e4, 1e4, 2e3, 200, 20, 5, 1, 1))
  expect_warning(graduate(steep, "age", "pop", open_ended = FALSE,
                          graph = FALSE), "negative")
})

test_that("inputs are validated", {
  g <- grouped(function(a) 1000 - 10 * a, k = 6)$data
  expect_error(graduate(g[1:4, ], "age", "pop", open_ended = FALSE),
               "at least 5")
  expect_error(graduate(g[c(1:3, 5:6), ], "age", "pop", open_ended = FALSE),
               "contiguous 5-year")
  ten <- data.frame(age = seq(0, 50, 10), pop = 100)
  expect_error(graduate(ten, "age", "pop"), "5-year")
  bad <- g; bad$pop[2] <- -1
  expect_error(graduate(bad, "age", "pop"), "non-negative")
  bad$pop[2] <- NA
  expect_error(graduate(bad, "age", "pop"), "missing")
  expect_error(graduate(g[6:1, ], "age", "pop"), "increasing")
  expect_error(graduate(g, "nope", "pop"), "Column not found")
  expect_error(graduate(g, "age", "pop", method = "grabill"), "should be one of")
  expect_error(graduate(data.frame(age = "young", pop = 1), "age", "pop"),
               "Could not read an age")
})

test_that("the object prints and plots", {
  data(grouped_pop, envir = environment())
  g <- graduate(grouped_pop, "age_col", "2021", open_ended = FALSE)
  expect_s3_class(g, "dem_graduation")
  expect_s3_class(plot(g), "ggplot")
  out <- capture.output(print(g))
  expect_true(any(grepl("Sprague", out)))
  expect_true(any(grepl("80 single years, ages 0-79", out)))

  bm <- graduate(grouped_pop, "age_col", "2021", method = "beers_modified",
                 open_ended = FALSE, graph = FALSE)
  expect_true(any(grepl("not reproduced exactly", capture.output(print(bm)))))
  expect_null(bm$plot)
  expect_error(plot(bm), "No plot available")
})
