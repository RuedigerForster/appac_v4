test_that("the plot functions return ggplot / patchwork objects", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- get_fit()
  expect_s3_class(plot_area_pressure(fit, 1, 1), "ggplot")
  expect_s3_class(plot_area_date(fit, 1, 1, show_changepoints = FALSE),
                  "ggplot")
  expect_s3_class(plot_area_pressure_fit(fit), "ggplot")

  skip_if_not_installed("patchwork")
  expect_s3_class(plot_residuals(fit, 1, 1), "patchwork")
})

test_that("plot selectors validate the sample / peak", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- get_fit()
  expect_error(plot_area_pressure(fit, sample = 99), "out of range")
  expect_error(plot_area_pressure(fit, peak = "nope"), "unknown peak")
})
