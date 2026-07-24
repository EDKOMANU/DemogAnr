d <- data.frame(
  age  = c("0-14","15-44","45-64","65+"),
  m1   = c(0.002, 0.001, 0.008, 0.060),
  m2   = c(0.003, 0.0015, 0.010, 0.065),
  n1   = c(300000, 400000, 200000, 100000),
  n2   = c(200000, 350000, 250000, 200000)
)

# Preston et al. (2001) Box 2.1: Sweden vs Kazakhstan, females 1992
box21 <- data.frame(
  age = c("0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
          "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79",
          "80-84","85+"),
  pop_Sw = c(0.0136,0.0524,0.0559,0.0548,0.0604,0.0655,0.0709,0.0641,0.0654,
             0.0703,0.0730,0.0552,0.0481,0.0493,0.0512,0.0508,0.0420,0.0321,0.0251),
  pop_K  = c(0.0200,0.0868,0.1011,0.0929,0.0828,0.0716,0.0843,0.0842,0.0704,
             0.0561,0.0327,0.0579,0.0347,0.0430,0.0295,0.0178,0.0172,0.0102,0.0068),
  m_Sw = c(0.00467,0.00008,0.00013,0.00014,0.00023,0.00030,0.00032,0.00050,
           0.00069,0.00117,0.00201,0.00305,0.00461,0.00759,0.01226,0.02026,
           0.03664,0.06815,0.15729),
  m_K  = c(0.02137,0.00162,0.00045,0.00037,0.00078,0.00108,0.00103,0.00132,
           0.00182,0.00288,0.00430,0.00571,0.01082,0.01392,0.02679,0.03998,
           0.05469,0.10159,0.18030)
)

test_that("standardize() reproduces Preston Box 2.1 (Sweden vs Kazakhstan)", {
  s <- standardize(box21, "age", rate1_col = "m_Sw", rate2_col = "m_K",
                   pop1_col = "pop_Sw", pop2_col = "pop_K",
                   labels = c("Sweden", "Kazakhstan"))
  # age-standardized crude death rates, average standard (Box 2.1)
  expect_equal(round(unname(s$standardized[1]), 2), 7.37)
  expect_equal(round(unname(s$standardized[2]), 2), 11.88)
  # Sweden's crude rate is higher despite a lower standardized rate
  expect_gt(s$crude[1], s$crude[2])
  expect_lt(s$standardized[1], s$standardized[2])
})

test_that("ASCDR difference equals the Kitagawa rate component (average standard)", {
  s <- standardize(box21, "age", rate1_col = "m_K", rate2_col = "m_Sw",
                   pop1_col = "pop_K", pop2_col = "pop_Sw")
  kd <- decompose_rates(box21, "age", "m_K", "m_Sw", "pop_K", "pop_Sw")
  expect_equal(unname(s$standardized[1] - s$standardized[2]),
               kd$rate_component, tolerance = 1e-6)
})

test_that("Kitagawa decomposition is exactly additive", {
  kd <- decompose_rates(d, "age", "m1", "m2", "n1", "n2", per = 1000)
  expect_equal(kd$rate_component + kd$composition_component, kd$total,
               tolerance = 1e-9)
  # total equals the difference in crude rates
  cdr1 <- sum(d$n1 / sum(d$n1) * d$m1) * 1000
  cdr2 <- sum(d$n2 / sum(d$n2) * d$m2) * 1000
  expect_equal(kd$total, cdr1 - cdr2, tolerance = 1e-9)
})

test_that("Arriaga LE decomposition sums to the e0 difference", {
  data(gphc2010, envir = environment())
  m <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", sex = "male")
  f <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths", sex = "female")
  ad <- decompose_LE(m$lifetable, f$lifetable)
  expect_equal(sum(ad$table$contribution), ad$total, tolerance = 1e-6)
  expect_equal(ad$direct + ad$indirect, ad$total, tolerance = 1e-6)
  expect_equal(ad$total, ad$e0_2 - ad$e0_1, tolerance = 1e-9)
})

test_that("Arriaga decomposition recovers a real sex gap (Ghana 2021)", {
  data(ghmort2021, envir = environment())
  male   <- subset(ghmort2021, Sex == "Male")
  female <- subset(ghmort2021, Sex == "Female")
  m <- lifetable(male,   age = "Age", pop = "Population", Dx = "Deaths",
                 sex = "male",   nax_method = "keyfitz")
  f <- lifetable(female, age = "Age", pop = "Population", Dx = "Deaths",
                 sex = "female", nax_method = "keyfitz")
  ad <- decompose_LE(m$lifetable, f$lifetable)
  # females outlive males by a few years, and it decomposes additively
  expect_gt(ad$total, 2)
  expect_lt(ad$total, 6)
  expect_equal(sum(ad$table$contribution), ad$total, tolerance = 1e-6)
})
