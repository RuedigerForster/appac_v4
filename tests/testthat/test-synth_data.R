# Tests for the bundled `Synth_data` stress fixture: verify its shape and its
# PLANTED ground truth (recoverable straight from the data, no fit needed), then
# that the appac() pipeline runs on it and behaves as documented -- the
# correction still reduces uncertainty, and the fitted kappa keeps the right
# sign and rough magnitude despite the unmodelled breakpoints biasing a
# whole-series fit.  `acn`, `kappa_of()` and `u_reduction()` come from the
# helper files.

test_that("Synth_data has the documented shape and columns", {
  expect_s3_class(Synth_data, "data.frame")
  expect_identical(dim(Synth_data), c(4500L, 6L))
  expect_identical(colnames(Synth_data),
                   c("sample.name", "injection.date", "peak.name",
                     "retention.time", "raw.area", "air.pressure"))
  expect_identical(sort(unique(Synth_data$sample.name)), c("S1", "S2", "S3"))
  expect_length(unique(Synth_data$peak.name), 10L)
  expect_length(unique(Synth_data$injection.date), 150L)          # 3 x 50 runs
  expect_false(anyNA(Synth_data$raw.area))
  expect_true(all(Synth_data$raw.area > 0))
})

test_that("the attached ground truth is complete and self-consistent", {
  tr <- attr(Synth_data, "truth")
  expect_type(tr, "list")
  expect_true(all(
    c("kappa", "P_ref", "breakpoints", "center_shift", "center_shift_bp2",
      "variance_expansion", "reference") %in% names(tr)))
  expect_length(tr$breakpoints, 2L)
  expect_s3_class(tr$breakpoints, "IDate")
  expect_lt(tr$kappa, 0)   # planted sensitivity is negative
  expect_identical(dim(tr$reference), c(30L, 3L))   # 3 samples x 10 peaks
  # within-sample dynamic range is the planted factor of 5
  rng <- tapply(tr$reference$ref.area, tr$reference$sample.name,
                function(x) max(x) / min(x))
  expect_equal(as.numeric(rng), rep(tr$dynamic_range, 3), tolerance = 1e-6)
})

test_that("the two planted breakpoints are present in the data", {
  # composition is FIXED; only the per-episode centre shifts multiplicatively.
  # Split each (sample, peak) series into its three 50-run episodes and recover
  # the centre ratios (robustly, against the heavy tails + contamination).
  D  <- Synth_data[order(Synth_data$injection.date), ]
  tr <- attr(Synth_data, "truth")
  r_bp1 <- r_bp2 <- spread <- numeric(0)
  for (s in unique(D$sample.name)) for (p in unique(D$peak.name)) {
    a  <- D$raw.area[D$sample.name == s & D$peak.name == p]   # date-ordered
    e1 <- a[1:50]; e2 <- a[51:100]; e3 <- a[101:150]
    r_bp1  <- c(r_bp1,  stats::median(e2) / stats::median(e1))
    r_bp2  <- c(r_bp2,  stats::median(e3) / stats::median(e2))
    spread <- c(spread, stats::mad(e3)    / stats::mad(e2))
  }
  # bp1: +5%
  expect_equal(stats::median(r_bp1), 1 + tr$center_shift,     tolerance = 0.01)
  # bp2: -2%
  expect_equal(stats::median(r_bp2), 1 + tr$center_shift_bp2, tolerance = 0.01)
  expect_gt(stats::median(spread), 1.05)   # bp2: variance up
})

test_that("appac() runs on Synth_data and the correction reduces uncertainty", {
  skip_on_cran()
  d   <- suppressMessages(check_cols(Synth_data, acn))
  tr  <- attr(Synth_data, "truth")
  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  expect_s4_class(fit, "Appac")
  expect_gt(stats::median(u_reduction(fit)), 1.2)   # correction still helps
})

test_that(paste(
  "fitted kappa keeps the right sign and rough magnitude",
  "(biased by the breaks)"), {
  skip_on_cran()
  d   <- suppressMessages(check_cols(Synth_data, acn))
  tr  <- attr(Synth_data, "truth")
  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  kap <- kappa_of(fit)
  expect_lt(kap, 0)                                               # correct sign
  # a whole-series fit is BIASED by the unmodelled centre steps, so the
  # tolerance
  # is deliberately loose -- this documents robustness, not clean recovery.
  expect_equal(kap, tr$kappa, tolerance = 0.5)
})

test_that("goodness_of_fit() returns per-sample per-peak diagnostics", {
  skip_on_cran()
  d   <- suppressMessages(check_cols(Synth_data, acn))
  tr  <- attr(Synth_data, "truth")
  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  gof <- goodness_of_fit(fit)
  expect_length(gof, 3L)
  expect_true(
    all(c("reduced.chisq", "chisq", "dof", "p.value") %in% colnames(gof[[1]])))
  expect_identical(nrow(gof[[1]]), 10L)
})

test_that(
  "the level and variance change-point detectors run and are deterministic", {
  skip_on_cran()
  skip_if_not_installed("strucchange")
  d   <- suppressMessages(check_cols(Synth_data, acn))
  tr  <- attr(Synth_data, "truth")
  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  # composition-preserving (common-mode) LEVEL breaks are hard for the PC2
  # detector, and the VARIANCE step is near the detection floor, so assert only
  # that the machinery runs and returns an IDate vector -- not exact recovery
  # (see ?Synth_data).
  lvl <- get_changepoints(fit@samples)
  vrc <- get_variance_changepoints(fit@samples)
  expect_s3_class(lvl, "IDate")
  expect_s3_class(vrc, "IDate")
  # strucchange is deterministic -- no MCMC, no seed
  expect_identical(get_changepoints(fit@samples), lvl)
  expect_identical(get_variance_changepoints(fit@samples), vrc)
})
