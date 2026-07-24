test_that("pyramid returns a ggplot", {
  data(ghpop2021, envir = environment())
  g <- within(ghpop2021, grp <- ifelse(Age >= 80, "80+",
              paste0((Age %/% 5) * 5, "-", (Age %/% 5) * 5 + 4)))
  pyr <- aggregate(Population ~ grp + Sex, data = g, FUN = sum)
  p <- pyramid(pyr, "grp", "Sex", "Population")
  expect_s3_class(p, "ggplot")
  expect_error(pyramid(pyr, "grp", "Sex", "Population", male_label = "M"),
               "match values")
})

test_that("graph = TRUE attaches a ggplot; graph = FALSE does not", {
  data(gphc2010, envir = environment())
  lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths")
  expect_s3_class(lt$plot, "ggplot")
  expect_s3_class(plot(lt), "ggplot")
  lt0 <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", graph = FALSE)
  expect_null(lt0$plot)
  # results are unchanged whether or not a plot is drawn
  expect_equal(lt$metrics$LifeExpectancyAtBirth, lt0$metrics$LifeExpectancyAtBirth)
})

test_that("age-quality and rate functions attach plots", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  expect_s3_class(whipple(male, "Age", "Population")$plot, "ggplot")
  expect_s3_class(myers(male, "Age", "Population", 10, 69)$plot, "ggplot")
  fr <- dem.fert(
    data.frame(age = c("15-19","20-24","25-29"), pop = c(1,1,1),
               w = c(5000,6000,5500), b = c(200,400,300)),
    type = "ASFR", age_col = "age", population_col = "pop",
    women_col = "w", births_col = "b")
  expect_s3_class(fr$plot, "ggplot")
})

test_that("project_population carries a plot and remains usable as a data frame", {
  data(region2000, envir = environment())
  pp <- project_population(region2000[1:3, ], base_year = 2000, future_year = 2004,
    TFR_var = "TFR", base_pop_var = "base_pop", region_var = "Country",
    death_rate_var = "death_rate", net_migration_var = "net_migration",
    subregion_var = "Region", num_samples = 300)
  expect_s3_class(plot(pp), "ggplot")
  expect_true(all(c("year", "median", "lower", "upper") %in% names(pp)))
  expect_true(nrow(subset(pp, year == 2000)) == 3)
})
