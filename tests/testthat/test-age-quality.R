# UN Manual II / Kpedekpo (1982) worked example: Ghana 1960 census.
ghana1960 <- data.frame(
  age = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
          "40-44","45-49","50-54","55-59","60-64","65+"),
  males = c(643041, 515520, 357831, 275542, 268336, 278601, 242515,
            198231, 168937, 122756, 96775, 59307, 63467, 113185),
  females = c(654258, 503070, 323460, 265534, 322576, 306329, 245883,
              179182, 145572, 95590, 81715, 48412, 54572, 100392)
)

test_that("un_age_sex_accuracy reproduces the Ghana 1960 worked example", {
  u <- un_age_sex_accuracy(ghana1960, "age", "males", "females")
  # Kpedekpo Table 3.3: male age-ratio score 7.9, female 11.9
  expect_equal(round(u$ARSM, 1), 7.9)
  expect_equal(round(u$ARSF, 1), 11.9)
  # sex ratio score ~8.7-8.8 (Kpedekpo rounds to 8.7 before x3)
  expect_equal(u$SRS, 8.77, tolerance = 0.05)
  # joint index: 46.2 exact; 45.9 in the book using SRS rounded to 8.7
  expect_equal(u$index, 3 * u$SRS + u$ARSM + u$ARSF)
  expect_equal(u$index, 45.9, tolerance = 0.5)
  expect_equal(u$quality, "highly inaccurate")
})

test_that("age_ratio computes the UN age ratio and excludes boundaries", {
  ar <- age_ratio(ghana1960[, c("age","males")], "age", "males")
  # males 5-9: 100*2*515520/(643041+357831) = 103.01
  expect_equal(round(ar$age_ratio[2], 2), 103.01)
  expect_true(is.na(ar$age_ratio[1]))    # first group
  expect_true(is.na(ar$age_ratio[14]))   # open group
  expect_true(is.na(ar$age_ratio[13]))   # adjacent to open group
})

test_that("sex_ratio computes males per 100 females", {
  sr <- sex_ratio(ghana1960, "age", "males", "females")
  expect_equal(round(sr$sex_ratio[1], 2), round(100*643041/654258, 2))
})

test_that("Whipple's index hits its analytic bounds", {
  uniform <- data.frame(age = 23:62, pop = rep(1000, 40))
  expect_equal(whipple(uniform, "age", "pop")$index, 100)
  heaped <- data.frame(age = 23:62, pop = ifelse((23:62) %% 5 == 0, 1000, 0))
  expect_equal(whipple(heaped, "age", "pop")$index, 500)
  expect_equal(whipple(heaped, "age", "pop")$quality, "very rough")
})

test_that("Myers' index is zero for a uniform distribution and bounded", {
  uniform <- data.frame(age = 0:99, pop = rep(1000, 100))
  expect_equal(myers(uniform, "age", "pop", 10, 89)$index, 0, tolerance = 1e-8)
  # heaping concentrated on digit 0 gives a large index
  heaped <- data.frame(age = 10:98, pop = ifelse((10:98) %% 10 == 0, 1e4, 100))
  mi <- myers(heaped, "age", "pop", 10, 89)$index
  expect_gt(mi, 50); expect_lt(mi, 90)
})

test_that("age-quality indices run on the bundled Ghana 2021 data", {
  data(ghpop2021, envir = environment())
  male <- subset(ghpop2021, Sex == "Male")
  expect_gt(whipple(male, "Age", "Population")$index, 100)
  expect_gte(myers(male, "Age", "Population", 10, 69)$index, 0)
})
