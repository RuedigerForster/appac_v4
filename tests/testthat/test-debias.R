# debias_ct(): the de-biased second pass.  Sweeps a scale factor around each
# centre, re-runs appac(), fits the per-peak chi-square as a parabola and takes
# its minimum.  Checked on a clean, seeded synthetic so the optimum stays inside
# the sweep and the test is stable.

test_that("debias_ct() returns per-sample centres that re-fit cleanly", {
  skip_on_cran()
  d    <- suppressMessages(check_cols(
    make_synthetic(kappa = -6e-4, noise = 1e-3, n_samples = 3, n_peaks = 4,
                   n_days = 250, seed = 3), acn))
  ap   <- as.numeric(d[, "Air_Pressure"]); pr <- (max(ap) + min(ap)) / 2
  fit1 <- suppressWarnings(appac(d, P_ref = pr))

  ct <- suppressWarnings(
    debias_ct(fit1, data = d, P_ref = pr, npt = 6, quiet = TRUE))

  # structure: a numeric vector of per-peak centres for every sample
  expect_type(ct, "list")
  expect_named(ct, names(fit1@samples))
  expect_identical(lengths(ct), lengths(fit1@trend@center))
  expect_true(
    all(vapply(ct, function(x) all(is.finite(x) & x > 0), logical(1))))

  # de-biased centres stay near the originals (the sweep is +/-1% around them)
  for (s in names(ct))
    expect_equal(unname(ct[[s]]), unname(fit1@trend@center[[s]]),
                 tolerance = 0.05)

  # feeding them back re-fits and does not blow up the goodness-of-fit
  fit2 <- suppressWarnings(appac(d, ct = ct, P_ref = pr))
  expect_s4_class(fit2, "Appac")
  g1 <- stats::median(goodness_of_fit(fit1)[[1]]$reduced.chisq)
  g2 <- stats::median(goodness_of_fit(fit2)[[1]]$reduced.chisq)
  expect_lt(g2, g1 * 1.5)
})
