# State-of-field comparison for DemogAnr's JOSS paper.
# Run this with a normal internet-connected R install (e.g. via Claude Code
# on your machine). It installs DemoTools (a well-established, actively
# maintained abridged-life-table package used across the demography
# community) and runs it on the SAME published worked example already used
# in DemogAnr's own vignette/tests (Preston, Heuveline & Guillot 2001, Box
# 4.1, US females 1991), then prints DemogAnr's output on the identical
# input directly underneath it.
#
# Goal: capture REAL console output from both packages on identical data,
# so the "State of the field" section in paper.md quotes actual evidence
# instead of asserting a readability difference (this is exactly what the
# R Journal reviewer flagged as missing).
#
# Save everything this script prints — copy the full console transcript
# into a file (e.g. paper/comparison_output.txt) and send it back.

if (!requireNamespace("DemoTools", quietly = TRUE)) {
  install.packages("DemoTools", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("DemogAnr", quietly = TRUE)) {
  stop("Install DemogAnr first (devtools::install_local(...) or load_all()).")
}

library(DemoTools)
library(DemogAnr)

## ---- Shared input: Preston et al. (2001) Box 4.1, US females 1991 ----
Age  <- c(0, 1, seq(5, 85, 5))
lx   <- c(100000, 99217, 99050, 98959, 98870, 98637, 98379, 98070, 97653,
          97083, 96289, 95008, 93018, 89882, 85249, 78711, 69618, 57486, 41756)
nqx  <- c(0.00783, 0.00168, 0.00092, 0.00090, 0.00236, 0.00262, 0.00314,
          0.00425, 0.00584, 0.00818, 0.01330, 0.02095, 0.03371, 0.05155,
          0.07669, 0.11552, 0.17427, 0.27363, 1.00000)
AgeInt <- c(diff(Age), NA)

cat("\n\n==================== DemoTools::lt_abridged() ====================\n\n")
dt_lt <- lt_abridged(nqx = nqx, Age = Age, AgeInt = AgeInt, axmethod = "un",
                      Sex = "f", OAG = TRUE)
print(dt_lt)
cat("\nColumn names as returned:\n")
print(names(dt_lt))
cat("\nTo get life expectancy at birth from this object:\n")
cat("dt_lt$ex[dt_lt$Age == 0]  ->", dt_lt$ex[dt_lt$Age == 0], "\n")

cat("\n\n==================== DemogAnr::lifetable_nqx() ====================\n\n")
us91 <- data.frame(Age = Age, nqx = nqx)
da_lt <- lifetable_nqx(us91, age = "Age", nqx = "nqx")
print(da_lt)
cat("\nTo get life expectancy at birth from this object:\n")
cat("da_lt$metrics$LifeExpectancyAtBirth ->", da_lt$metrics$LifeExpectancyAtBirth, "\n")

cat("\n\n==================== Session info (for the record) ====================\n\n")
print(sessionInfo())
