test_that("goodness_of_fit() returns per-sample reduced chi-square", {
  skip_on_cran()
  fit <- get_fit()
  g   <- goodness_of_fit(fit)

  expect_type(g, "list")
  expect_named(g, names(fit@samples))
  expect_true(
    all(c("reduced.chisq", "chisq", "dof", "p.value") %in% colnames(g[[1]])))
  expect_true(all(g[[1]]$reduced.chisq > 0))
})
