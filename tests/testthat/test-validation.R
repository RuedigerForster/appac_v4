# Input validation of appac(): every stop() guard should fire on the bad input.
# These hit before any fitting, so they are fast and always run.

test_that("appac() rejects malformed input", {
  d  <- suppressMessages(check_cols(
    make_synthetic(n_samples = 2, n_peaks = 3, n_days = 120, seed = 5), acn))
  ap <- as.numeric(d[, "Air_Pressure"])
  pr <- (max(ap) + min(ap)) / 2

  # missing data argument
  expect_error(appac(), "missing")

  # unknown / non-canonical column names
  bad <- d; colnames(bad)[1] <- "Wrong_Name"
  expect_error(appac(bad, P_ref = pr), "Unknown column names")

  # P_ref outside the observed pressure range (either side) or non-numeric
  expect_error(appac(d, P_ref = max(ap) + 50), "out of range")
  expect_error(appac(d, P_ref = min(ap) - 50), "out of range")
  expect_error(appac(d, P_ref = "not a number"), "out of range")

  # incomplete grid: nrow no longer a multiple of the number of peaks
  expect_error(appac(d[-1, ], P_ref = pr), "Missing data points")
})

test_that(
  "appac() guards undersized and degenerate input with clear messages", {
  ok <- function(...) suppressMessages(check_cols(make_synthetic(...), acn))

  # < 2 peaks
  expect_error(
    appac(ok(n_samples = 3, n_peaks = 1, n_days = 80, seed = 1), P_ref = 1000),
    "at least 2 peaks")
  # < 3 samples (the cross-sample drift PCA needs three)
  expect_error(
    appac(ok(n_samples = 2, n_peaks = 3, n_days = 80, seed = 1), P_ref = 1000),
    "at least 3 samples")
  # too few injections per sample
  expect_error(
    appac(ok(n_samples = 3, n_peaks = 3, n_days = 15, seed = 1), P_ref = 1000),
    "at least 20 injections")
  # a constant (zero-variance) peak
  z <- make_synthetic(n_samples = 3, n_peaks = 3, n_days = 80, seed = 1)
  z$raw.area[z$peak.name == "pk2"] <- 500
  expect_error(appac(suppressMessages(check_cols(z, acn)), P_ref = 1000),
               "no variation")
})
