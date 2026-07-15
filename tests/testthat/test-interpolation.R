test_that("interpolation is exact on linear data", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 4, 6, 8, 10)
  out <- interpolation(x, y, x_new = c(2.5, 3.5), method = "linear")
  expect_equal(unname(out), c(5, 7))
  # auto (spline for >3 points) is also exact on linear data
  out2 <- interpolation(x, y, x_new = c(2.5, 3.5), method = "auto")
  expect_equal(unname(out2), c(5, 7))
})

test_that("interpolation quadratic is exact on quadratic data", {
  x <- c(0, 1, 2)
  y <- x^2
  out <- interpolation(x, y, x_new = 1.5, method = "quadratic")
  expect_equal(unname(out), 2.25)
})

test_that("interpolation validates inputs", {
  expect_error(interpolation(c(1, 1, 2), c(1, 2, 3), 1.5), "duplicate")
  expect_error(interpolation(c(1, 2), c(1, 2), 1.5, method = "quadratic"), "3 points")
  expect_warning(interpolation(c(1, 2), c(1, 2), 5, method = "linear"), "outside")
})

test_that("interpolation works on data frames", {
  d <- data.frame(category = c("A", "B"), y1 = c(2, 5), y2 = c(3, 6))
  out <- interpolation(x = c(1, 2), data = d, y_cols = c("y1", "y2"),
                       x_new = 1.5, method = "linear")
  expect_equal(out$y_1.5, c(2.5, 5.5))
})
