test_that(
  "appac() returns an Appac and recovers kappa near the validated value", {
  skip_on_cran()
  fit   <- get_fit()
  kappa <- unname(unlist(fit@correction@coefficients))

  expect_s4_class(fit, "Appac")
  expect_length(kappa, 1L)
  expect_lt(kappa, 0)   # pressure sensitivity is negative
  expect_equal(kappa, -7.1e-4, tolerance = 0.15)   # ~ shrubbery-validated value
  expect_gte(length(fit@correction@fit.data), 1L)  # kappa-fit input is stored
})

test_that("the correction reduces run-to-run scatter (RSD) for every peak", {
  skip_on_cran()
  fit <- get_fit()
  rsd <- function(x) stats::sd(x) / mean(x)
  raw <- fit@samples[[1]]$raw.area
  cor <- fit@samples[[1]]$corrected.area
  expect_true(all(apply(cor, 2, rsd) < apply(raw, 2, rsd)))
})
