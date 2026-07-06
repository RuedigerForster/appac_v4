# show() methods: a fitted Appac must print a compact summary, not dump the
# multi-megabyte area matrices.

test_that("show() prints a compact summary for the S4 classes", {
  skip_on_cran()
  d   <- suppressMessages(check_cols(
    make_synthetic(n_samples = 3, n_peaks = 3, n_days = 120, seed = 1), acn))
  ap  <- as.numeric(d[, "Air_Pressure"])
  fit <- suppressWarnings(appac(d, P_ref = (max(ap) + min(ap)) / 2))

  out <- capture.output(show(fit))
  expect_lt(length(out), 12L)   # compact, not the ~75k-line dump
  expect_match(out[1], "Appac")
  expect_true(any(grepl("samples", out)))
  expect_true(any(grepl("kappa",   out)))

  expect_match(capture.output(show(fit@trend))[1],      "Compensation")
  expect_match(capture.output(show(fit@correction))[1], "Correction")

  expect_invisible(show(fit))   # returns the object invisibly
})

test_that("print() shows the summary plus goodness-of-fit", {
  skip_on_cran()
  d   <- suppressMessages(check_cols(
    make_synthetic(n_samples = 3, n_peaks = 3, n_days = 120, seed = 1), acn))
  ap  <- as.numeric(d[, "Air_Pressure"])
  fit <- suppressWarnings(appac(d, P_ref = (max(ap) + min(ap)) / 2))

  out <- capture.output(print(fit))
  expect_match(out[1], "Appac")                      # reuses the show header
  expect_true(any(grepl("goodness of fit", out)))    # ... plus GOF
  expect_true(any(grepl("reduced chi-square", out)))
  expect_invisible(print(fit))
})
