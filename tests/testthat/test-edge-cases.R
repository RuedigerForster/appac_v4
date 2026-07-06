# The latent input-shape cases that were hardened (commit history): duplicate /
# missing dates in get_daily_value, and a cylinder that collapses to one peak.

test_that(".get_daily_value tolerates duplicate and missing dates", {
  dv  <- data.frame(date = c(1, 2, 2, 4), value = c(10, 20, 21, 40))
  res <- appac:::.get_daily_value(c(1, 2, 3, 4), dv)
  # dup -> first match, missing -> NA
  expect_equal(as.numeric(res), c(10, 20, NA, 40))
})

test_that(".prepare keeps a single-peak cylinder as a matrix", {
  set.seed(1); n <- 60
  df <- data.frame(
    Sample_Name    = "CTL.1",
    Peak_Name      = rep(c("p1", "p2"), each = n),
    Injection_Date = rep(1:n, 2),
    Air_Pressure   = rep(rnorm(n, 1000, 5), 2),
    # p2 all-NA -> dropped
    Raw_Area       = c(rnorm(n, 1e6, 1e3), rep(NA_real_, n)),
    stringsAsFactors = FALSE)
  A  <- suppressWarnings(
    appac:::.prepare(df, P_ref = 1000, covariates = "Air_Pressure"))
  ra <- A@samples[["CTL.1"]]$raw.area
  expect_true(is.matrix(ra))
  expect_equal(ncol(ra), 1L)
})
