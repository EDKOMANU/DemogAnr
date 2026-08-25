# Guards against side effects a package must not have: drawing on the user's
# graphics device from print(), and disturbing the user's random stream.

test_that("print() does not open a graphics device or write Rplots.pdf", {
  data(gphc2010, envir = environment())
  wd <- file.path(tempdir(), paste0("dg", as.integer(runif(1, 1e6, 9e6))))
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  old <- setwd(wd); on.exit({ setwd(old); unlink(wd, recursive = TRUE) }, add = TRUE)

  before <- length(grDevices::dev.list())
  res <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
  expect_s3_class(res$plot, "ggplot")          # the plot is still attached
  invisible(capture.output(print(res)))

  expect_equal(length(grDevices::dev.list()), before)
  expect_false(file.exists(file.path(wd, "Rplots.pdf")))
})

test_that("every print() method leaves the graphics device alone", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  wd <- file.path(tempdir(), paste0("dg", as.integer(runif(1, 1e6, 9e6))))
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  old <- setwd(wd); on.exit({ setwd(old); unlink(wd, recursive = TRUE) }, add = TRUE)

  objs <- list(
    whipple(male, "Age", "Population"),
    myers(male, "Age", "Population", lower = 10, upper = 69),
    dem.fert(data.frame(age = c("15-19", "20-24"), pop = c(1, 1),
                        w = c(5000, 6000), b = c(200, 400)),
             type = "ASFR", age_col = "age", population_col = "pop",
             women_col = "w", births_col = "b")
  )
  before <- length(grDevices::dev.list())
  for (o in objs) invisible(capture.output(print(o)))
  expect_equal(length(grDevices::dev.list()), before)
  expect_false(file.exists(file.path(wd, "Rplots.pdf")))
})

test_that("plot() still returns the attached figure", {
  data(gphc2010, envir = environment())
  res <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
  expect_s3_class(plot(res), "ggplot")
  res0 <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
                    graph = FALSE)
  expect_error(plot(res0), "No plot available")
})

test_that("project_population() leaves the caller's random stream untouched", {
  pp <- data.frame(region = "R", subregion = "S", base_pop = 1e5,
                   cbr = 0.03, cdr = 0.01, nmr = 0.001)

  set.seed(1); expected <- runif(3)
  set.seed(1)
  invisible(project_population(pp, future_year = 2030, base_year = 2025,
                               num_samples = 50, graph = FALSE))
  expect_equal(runif(3), expected)
})

test_that("calibrate_regions() leaves the caller's random stream untouched", {
  d <- data.frame(reg = c("A", "B"), pop = c(600, 400))
  nat <- data.frame(yr = c(2020, 2025), tot = c(1000, 1100))

  set.seed(7); expected <- runif(3)
  set.seed(7)
  invisible(calibrate_regions(d, "reg", "pop", nat, "yr", "tot",
                              deviation = 0.05, num_samples = 20,
                              graph = FALSE))
  expect_equal(runif(3), expected)
})

test_that("seeding still makes results reproducible, and random_seed = NULL opts out", {
  pp <- data.frame(region = "R", subregion = "S", base_pop = 1e5,
                   cbr = 0.03, cdr = 0.01, nmr = 0.001)
  a <- project_population(pp, 2030, 2025, num_samples = 200, graph = FALSE)
  b <- project_population(pp, 2030, 2025, num_samples = 200, graph = FALSE)
  expect_equal(a$median, b$median)

  set.seed(11)
  c1 <- project_population(pp, 2030, 2025, num_samples = 200,
                           random_seed = NULL, graph = FALSE)
  c2 <- project_population(pp, 2030, 2025, num_samples = 200,
                           random_seed = NULL, graph = FALSE)
  expect_false(isTRUE(all.equal(c1$median, c2$median)))
})

test_that("an uninitialised random stream is left uninitialised", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  pp <- data.frame(region = "R", subregion = "S", base_pop = 1e5,
                   cbr = 0.03, cdr = 0.01, nmr = 0.001)
  invisible(project_population(pp, 2027, 2025, num_samples = 20, graph = FALSE))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})
