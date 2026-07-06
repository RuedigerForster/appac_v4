# Missing-value imputation.  Two distinct paths, both exercised at ~20% loss:
#   B) per-CELL NA in the areas -> filled by .impute_na() in .prepare()
#   (svdImpute
#      / EM low-rank reconstruction) before the loess + PCA ever see them;
#   A) whole missing INJECTIONS -> samples fall on different dates, so the
#      cross-sample daily-factor matrix has holes that .na_pca() reconstructs.
# `acn`, `make_synthetic()`, `kappa_of()` come from the helper files.

## --- B: the .impute_na() engine -------------------------------------------

test_that(paste(
  ".impute_na() fills gaps by low-rank reconstruction and leaves the",
  "rest intact"), {
  set.seed(1)
  U <- matrix(stats::rnorm(40 * 2), 40, 2)
  V <- matrix(stats::rnorm(6 * 2),  6, 2)
  # rank-2 + offset + noise
  M <- U %*% t(V) * 100 + 1000 + matrix(stats::rnorm(240, 0, 2), 40, 6)
  Mna <- M
  na  <- sample(length(M), round(0.20 * length(M)))   # 20% cells missing
  Mna[na] <- NA

  imp <- appac:::.impute_na(Mna, ncp = 2L)
  expect_false(anyNA(imp))                              # every gap filled
  # observed cells untouched
  expect_identical(imp[!is.na(Mna)], M[!is.na(Mna)])
  # good reconstruction
  expect_lt(sqrt(mean((imp[na] - M[na])^2)) / stats::sd(M), 0.10)
})

test_that(".impute_na() is a no-op on complete input", {
  m <- matrix(1:20 + 0, 5, 4)
  expect_identical(appac:::.impute_na(m), m)
})

## --- B: per-cell NA carried end-to-end through appac() ---------------------

test_that(
  "appac() imputes 20% per-cell NA and still recovers kappa on clean data", {
  skip_on_cran()
  # a well-conditioned set (no breakpoints, light noise) isolates imputation
  # quality from the Synth_data torture conditions.
  base <- make_synthetic(kappa = -8e-4, noise = 1e-3, n_samples = 5,
                         n_peaks = 6, n_days = 400, seed = 1)
  set.seed(42)
  base$raw.area[sample(nrow(base), round(0.20 * nrow(base)))] <- NA
  d   <- suppressMessages(check_cols(base, acn))
  ap  <- as.numeric(d[, "Air_Pressure"]); pr <- (max(ap) + min(ap)) / 2
  fit <- suppressWarnings(appac(d, P_ref = pr))
  # imputation preserves the fit
  expect_equal(kappa_of(fit), -8e-4, tolerance = 0.10)
  expect_false(
    any(vapply(fit@samples, function(s) anyNA(s$corrected.area), logical(1))))
})

test_that(
  "appac() is robust to 20% per-cell NA on the Synth_data torture set", {
  skip_on_cran()
  # here the point is robustness, not recovery: the breakpoints already bias the
  # whole-series fit, and heavy imputation on top makes kappa unreliable -- so
  # we assert only that it runs and returns a complete (NA-free) correction.
  tr <- attr(Synth_data, "truth")
  set.seed(42)
  D <- Synth_data
  D$raw.area[sample(nrow(D), round(0.20 * nrow(D)))] <- NA
  d   <- suppressMessages(check_cols(D, acn))
  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  expect_s4_class(fit, "Appac")
  expect_false(
    any(vapply(fit@samples, function(s) anyNA(s$corrected.area), logical(1))))
})

## --- A: missing whole injections -> cross-sample .na_pca() imputation -------

test_that(
  "appac() handles 20% missing injections (staggered dates -> .na_pca)", {
  skip_on_cran()
  tr <- attr(Synth_data, "truth")
  set.seed(7)
  D <- Synth_data
  drop <- logical(nrow(D))
  for (s in unique(D$sample.name)) {   # drop 20% of each sample's dates
    ds   <- unique(D$injection.date[D$sample.name == s])
    gone <- sample(ds, round(0.20 * length(ds)))
    drop <- drop | (D$sample.name == s & D$injection.date %in% gone)
  }
  D <- D[!drop, ]
  d <- suppressMessages(check_cols(D, acn))

  # samples now fall on DIFFERENT dates, so the daily-factor matrix has holes
  # that only .na_pca() can fill -- a successful fit is proof the path works.
  per_sample <- tapply(d$Injection_Date, d$Sample_Name,
                       function(x) length(unique(x)))
  union_dates <- length(unique(d$Injection_Date))
  expect_true(union_dates > max(per_sample))   # dates are genuinely misaligned

  fit <- suppressWarnings(appac(d, P_ref = tr$P_ref))
  expect_s4_class(fit, "Appac")
  kap <- kappa_of(fit)
  expect_lt(kap, 0)   # correct sign survives the gaps
  expect_equal(kap, tr$kappa, tolerance = 0.3)
})
