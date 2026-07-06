test_that("check_cols() is independent of the appac_colnames order", {
  shuf <- acn[c("area_col", "date_col", "sample_col", "pressure_col",
                "peak_col")]
  d1 <- suppressMessages(check_cols(PLOT_FID, acn))
  d2 <- suppressMessages(check_cols(PLOT_FID, shuf))
  expect_identical(d1, d2)
  expect_identical(colnames(d1),
    c("Sample_Name", "Peak_Name", "Injection_Date", "Air_Pressure", "Raw_Area"))
})

test_that("check_cols() errors when make.names would merge distinct peaks", {
  bad <- data.frame(sample.name = "S1", peak.name = c("a-b", "a.b"),
                    injection.date = 1:2, air.pressure = 1000,
                    raw.area = c(1, 2), stringsAsFactors = FALSE)
  expect_error(suppressMessages(check_cols(bad, acn)), "collide")
})

test_that("check_cols() errors on a missing required column", {
  bad <- data.frame(sample.name = "S1", peak.name = "p1",
                    injection.date = 1, air.pressure = 1000,
                    stringsAsFactors = FALSE)
  expect_error(suppressMessages(check_cols(bad, acn)))
})
