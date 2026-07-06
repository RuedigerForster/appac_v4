# Correctness against ground truth (known kappa) and the headline objective:
# demonstrate an >= 8x (1/8) reduction in the standard uncertainty u.
# Tolerances are tight, calibrated just above the measured achievable error;
# the synthetic is seeded, so this is demanding but not flaky.

test_that("appac recovers a planted kappa to <3% (clean)", {
  skip_on_cran()
  for (kappa_true in c(-5e-4, -1e-3, 3e-4)) {   # magnitudes and both signs
    fit <- fit_synthetic(kappa = kappa_true, noise = 1e-3)
    expect_equal(kappa_of(fit), kappa_true, tolerance = 0.03)
  }
})

test_that(
  "kappa is recovered despite slow drift (drift-then-kappa de-confounding)", {
  skip_on_cran()
  fit <- fit_synthetic(kappa = -6e-4, noise = 1e-3, drift = 0.01)
  expect_equal(kappa_of(fit), -6e-4, tolerance = 0.03)
})

test_that(
  "kappa is recovered despite heavy-tailed noise + outliers (robust fit)", {
  skip_on_cran()
  fit <- fit_synthetic(kappa = -6e-4, noise = 1e-3, heavy = TRUE, contam = 0.02)
  expect_equal(kappa_of(fit), -6e-4, tolerance = 0.10)
})

test_that("the correction delivers a 1/8 (>= 8x) reduction in uncertainty u", {
  skip_on_cran()
  # white noise ~ 10% of the total uncertainty (the rest correctable pressure):
  # the noise floor caps the reduction at a credible ~10x, not an artificial
  # 100x.
  fit <- fit_synthetic(kappa = -8e-4, noise = 1e-3)
  red <- u_reduction(fit)   # raw RSD / corrected RSD, per peak
  expect_gte(stats::median(red), 8)             # >= 8x  ==  u shrinks to <= 1/8

  # and the residual is at the noise floor: reduced chi-square ~ 1
  rc <- goodness_of_fit(fit)[[1]]$reduced.chisq
  expect_gt(stats::median(rc), 0.8)
  expect_lt(stats::median(rc), 1.3)
})
