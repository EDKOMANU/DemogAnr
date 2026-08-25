## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 6, fig.height = 4, fig.align = "center"
)
options(width = 90)
library(DemogAnr)

## ----fert-------------------------------------------------------------------------------
bangladesh <- data.frame(
  age         = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"),
  women       = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
  live_births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125),
  population  = 71315944
)
fert <- dem.fert(bangladesh, type = "all", age_col = "age",
                 population_col = "population", women_col = "women",
                 births_col = "live_births", graph = FALSE)
fert$results$TFR

## ----cdr--------------------------------------------------------------------------------
sweden88 <- data.frame(population = 8438477, deaths = 96756)
dem.cdr(sweden88, type = "CDR", population_col = "population",
        deaths_col = "deaths", graph = FALSE)$results$CDR

## ----lt---------------------------------------------------------------------------------
lt <- lifetable(gphc2010, age = "Age", pop = "Pop", Dx = "Deaths",
                graph = FALSE)
head(lt$lifetable)
lt$metrics$LifeExpectancyAtBirth

## ----mdlt-------------------------------------------------------------------------------
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
           2.685, 2.681, 2.655, 2.647, 2.646, 2.631, 2.628, 2.618, 2.570, 6.539)
)
us91$Dother <- us91$Dall - us91$Dneo

md <- multiple_decrement(us91, age = "Age", causes = c("Dneo", "Dother"),
                         deaths = "Dall", lx = "lx", nqx = "nqx", graph = FALSE)
round(md$cause_distribution, 3)

## ----cdlt-------------------------------------------------------------------------------
cd <- cause_deleted_lt(us91, age = "Age", cause = "Dneo", deaths = "Dall",
                       lx = "lx", nqx = "nqx", nax_all = "nax", graph = FALSE)
c(e0 = cd$e0, without_neoplasms = cd$e0_deleted, gain = cd$gain)

## ----repro------------------------------------------------------------------------------
us1991 <- data.frame(
  age   = seq(10, 45, 5),
  fbir  = c(5816, 253979, 532712, 596823, 431694, 162005, 25531, 829),
  women = c(8620, 8371, 9419, 10325, 11125, 10344, 9496, 7188) * 1000,
  nLx   = c(494603, 493804, 492552, 491138, 489356, 486941, 483577, 478475)
)
r <- reproduction(us1991, age = "age", births = "fbir", women = "women",
                  nLx = "nLx", radix = 1e5, fraction_female = 1, graph = FALSE)
c(GRR = r$GRR, NRR = r$NRR)

## ----stable-----------------------------------------------------------------------------
egypt <- data.frame(
  age = seq(15, 45, 5),
  Lx  = c(4.66740, 4.63097, 4.58518, 4.53206, 4.46912, 4.39135, 4.28969),
  ma  = c(0.00567, 0.06627, 0.11204, 0.07889, 0.05075, 0.01590, 0.00610)
)
sp <- stable_population(egypt, age = "age", asfr = "ma", nLx = "Lx",
                        radix = 1, fraction_female = 1, graph = FALSE)
c(r = sp$r, b = sp$b, d = sp$d, NRR = sp$NRR)

## ----pf---------------------------------------------------------------------------------
bd <- data.frame(
  age    = seq(15, 45, 5),
  women  = c(3014706, 2653155, 2607009, 2015663, 1771680, 1479505, 1135129),
  ceb    = c(1160919, 4901382, 9085852, 9910256, 10384001, 9164329, 6905673),
  births = c(320406, 609269, 561494, 367833, 237297, 95357, 38125)
)
pf <- pf_ratio(bd, age = "age", women = "women", ceb = "ceb",
               births = "births", graph = FALSE)
c(K = pf$K, TFR = pf$TFR, TFR_adjusted = pf$TFR_adjusted)

## ----chm--------------------------------------------------------------------------------
panama <- data.frame(
  women = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
  ceb   = c(557, 2633, 4757, 6085, 6722, 6367, 5276),
  cd    = c(40, 130, 312, 435, 636, 686, 689)
)
est <- chm_brass(panama, women_col = "women", ceb_col = "ceb",
                 cd_col = "cd", model = "west", survey_year = 1976.7)
est$qx[est$x == 5]

## ----mlt--------------------------------------------------------------------------------
fitted <- model_lifetable("west", sex = "male", q5 = est$qx[est$x == 5],
                          graph = FALSE)
fitted$metrics$LifeExpectancyAtBirth

## ----mlt-families-----------------------------------------------------------------------
sapply(c("CD_West", "CD_North", "CD_South", "CD_East"), function(f)
  model_lifetable(f, sex = "male", q5 = est$qx[est$x == 5],
                  graph = FALSE)$metrics$LifeExpectancyAtBirth)

## ----brass------------------------------------------------------------------------------
fit <- brass_logit(data.frame(age = est$x, qx = est$qx),
                   qx_col = "qx", age_col = "age", standard = "African",
                   complete = TRUE)
round(fit$coefficients, 4)
fit$lifetable$metrics$LifeExpectancyAtBirth

## ----brass-target-----------------------------------------------------------------------
brass_lifetable(target_e0 = 60, standard = "African",
                graph = FALSE)$brass$alpha

## ----orphanhood-------------------------------------------------------------------------
bolivia <- data.frame(
  age   = seq(15, 50, 5),
  alive = c(5540, 3995, 2886, 1910, 1272, 855, 541, 234),
  dead  = c(448, 541, 769, 852, 945, 1059, 985, 788),
  s10   = c(NA, 0.9060, 0.8401, 0.7474, 0.6313, 0.5175, 0.3996, NA)
)
orph <- orphanhood(bolivia, age = "age", alive = "alive", dead = "dead",
                   prop10 = "s10", M = 28.8, survey_year = 1975.5,
                   graph = FALSE)
orph$table[, c("age", "W", "lx_ratio", "t", "ref_date")]

## ----ddm--------------------------------------------------------------------------------
d <- data.frame(
  age = c(0, seq(5, 85, 5)),
  n1  = c(41255, 36289, 31948, 28121, 24743, 21758, 19116, 16770, 14677,
          12798, 11094, 9527, 8060, 6659, 5298, 3970, 2706, 2614),
  n2  = c(52973, 46596, 41022, 36108, 31770, 27938, 24545, 21533, 18846,
          16433, 14245, 12233, 10349, 8550, 6802, 5098, 3475, 3356),
  dth = c(516, 134, 126, 123, 126, 136, 154, 184, 229, 294, 385, 508,
          668, 865, 1084, 1282, 1380, 2833)
)
g <- ggb(d, age = "age", pop1 = "n1", pop2 = "n2", deaths = "dth",
         t1 = 2000, t2 = 2010, graph = FALSE)
c(completeness = g$completeness, coverage_ratio = g$coverage_ratio)

## ----ddm-compare------------------------------------------------------------------------
c(GGB = g$completeness,
  SEG = seg(d, "age", "n1", "n2", "dth", 2000, 2010, graph = FALSE)$completeness,
  combined = ggb_seg(d, "age", "n1", "n2", "dth", 2000, 2010,
                     graph = FALSE)$completeness)

## ----std--------------------------------------------------------------------------------
male   <- subset(ghmort2021, Sex == "Male")
female <- subset(ghmort2021, Sex == "Female")
comp <- data.frame(
  age  = male$Age,
  m_m  = male$Deaths / male$Population,
  m_f  = female$Deaths / female$Population,
  n_m  = male$Population,
  n_f  = female$Population
)
standardize(comp, "age", rate1_col = "m_m", rate2_col = "m_f",
            pop1_col = "n_m", pop2_col = "n_f",
            labels = c("Males", "Females"))

## ----ccm--------------------------------------------------------------------------------
sweden <- data.frame(
  Age = c(seq(0, 80, 5), 85),
  N93 = c(293395, 248369, 240012, 261346, 285209, 314388, 281290, 286923,
          304108, 324946, 247613, 211351, 215140, 221764, 223506, 183654,
          141990, 112424),
  Lx  = c(497487, 497138, 496901, 496531, 495902, 495168, 494213, 492760,
          490447, 486613, 480665, 471786, 457852, 436153, 402775, 350358,
          271512, 291707),
  Fx  = c(0, 0, 0, 0.0120, 0.0908, 0.1499, 0.1125, 0.0441, 0.0074, 0.0003,
          rep(0, 8))
)
cc <- cohort_component(sweden, age = "Age", population = "N93", nLx = "Lx",
                       asfr = "Fx", steps = 2, graph = FALSE)
round(cc$totals)

## ----pproj------------------------------------------------------------------------------
pp <- project_population(region2000[1:4, ],
        base_year = 2000, future_year = 2010,
        region_var = "Country", subregion_var = "Region",
        base_pop_var = "base_pop", birth_rate_var = "cbr",
        death_rate_var = "cdr", net_migration_var = "nmr",
        num_samples = 500, graph = FALSE)
subset(pp, year == 2010)

## ----pproj-trend------------------------------------------------------------------------
pp_trend <- project_population(region2000[1:4, ],
        base_year = 2000, future_year = 2030,
        region_var = "Country", subregion_var = "Region",
        base_pop_var = "base_pop", birth_rate_var = "cbr",
        death_rate_var = "cdr", net_migration_var = "nmr",
        birth_target = 0.024, death_target = 0.008,
        cv = 0.05, target_cv = 0.20,
        num_samples = 500, graph = FALSE)
subset(pp_trend, year == 2030)

## ----calreg-----------------------------------------------------------------------------
base_total <- sum(region2000$base_pop)
national <- data.frame(
  year  = seq(2000, 2020, 5),
  total = round(base_total * exp(0.025 * seq(0, 20, 5)))
)
cal <- calibrate_regions(region2000,
        region_col = "Region", pop_col = "base_pop",
        national = national, year_col = "year", total_col = "total",
        deviation = 0.08, num_samples = 500, graph = FALSE)
cal$check

## ----quality----------------------------------------------------------------------------
males <- subset(ghpop2021, Sex == "Male")
w <- whipple(males, age_col = "Age", pop_col = "Population", graph = FALSE)
w$index

## ----smam-------------------------------------------------------------------------------
turkey <- data.frame(
  age    = seq(15, 50, 5),
  men    = c(3165061, 2581153, 2435765, 2096899, 1784121, 1418784,
             1111113, 980115),
  single = c(3030203, 1853222, 629077, 180767, 77134, 43412, 28627, 22527)
)
smam(turkey, age = "age", single = "single", total = "men")$SMAM

## ----pyramid, fig.height = 5------------------------------------------------------------
g <- within(ghpop2021, grp <- ifelse(Age >= 80, "80+",
            paste0((Age %/% 5) * 5, "-", (Age %/% 5) * 5 + 4)))
pyr <- aggregate(Population ~ grp + Sex, data = g, FUN = sum)
pyramid(pyr, age_col = "grp", sex_col = "Sex", count_col = "Population",
        title = "Ghana, 2021")

