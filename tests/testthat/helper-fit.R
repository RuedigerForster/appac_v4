# Shared fixtures for the test suite.  PLOT_FID is the package's bundled,
# sanitised demo dataset (lazy-loaded with the package).

acn <- list(sample_col = "sample.name", peak_col = "peak.name",
            date_col = "injection.date", pressure_col = "air.pressure",
            area_col = "raw.area")

# column-checked data + reference pressure
prep_data <- function() {
  d  <- suppressMessages(check_cols(PLOT_FID, acn))
  ap <- as.numeric(d[, "Air_Pressure"])
  list(data  = d,
       P_ref = mean(range(ap, na.rm = TRUE)))
}

# Fit appac() ONCE and cache it: the fit is the slow part, and most tests share
# the same result.  Tests that use this call skip_on_cran() first.
.cache <- new.env(parent = emptyenv())
get_fit <- function() {
  if (is.null(.cache$fit)) {
    pd <- prep_data()
    .cache$data  <- pd$data
    .cache$P_ref <- pd$P_ref
    .cache$fit   <- suppressWarnings(appac(pd$data, P_ref = pd$P_ref))
  }
  .cache$fit
}
