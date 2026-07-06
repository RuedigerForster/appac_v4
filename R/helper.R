## svdImpute / EM: fill the NA cells of a (caller-scaled) matrix by iterated
## rank-`ncp` SVD reconstruction -- gaps seeded at 0, iterated until the fill
## stabilises (rmse change < tol).  Returns the completed matrix in the SAME
## scaling the caller supplied; complete matrices are returned untouched.
## Shared
## by .na_pca() (unit-variance, uncentred daily factors) and .impute_na()
## (centred + unit-variance area matrices).
.svd_impute <- function(xs, ncp, maxiter = 200L, tol = 1e-7) {
  na <- is.na(xs)
  if (!any(na)) return(xs)
  xs[na] <- 0
  prev <- Inf
  for (it in seq_len(maxiter)) {
    sv    <- svd(xs)
    k     <- seq_len(ncp)
    recon <- sv$u[, k, drop = FALSE] %*% diag(sv$d[k], length(k)) %*%
      t(sv$v[, k, drop = FALSE])
    rmse  <- sqrt(mean((xs[na] - recon[na])^2))
    xs[na] <- recon[na]
    if (abs(prev - rmse) < tol) break
    prev <- rmse
  }
  xs
}

.summarize_factors <- function(X) {
  X %>%
    group_by(date) %>%
    summarize(across(where(is.numeric), mean))
}

## Evaluate `expr` silently and RETURN ITS VALUE: sinks stdout and stderr to the
## null device (so even functions that print directly, e.g. glmrob's
## "eliminating columns" diagnostic, stay quiet) and suppresses
## warnings/messages.
.quiet <- function(expr) {
  con <- file(if (.Platform$OS.type == "windows") "NUL" else "/dev/null",
              open = "wt")
  sink(con); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

## Map per-date `daily_value` (a data.frame or matrix whose first column is the
## date) onto the requested `date` vector.  Uses match() -- the first matching
## row per date, NA if the date is absent -- so duplicate or missing dates yield
## NA rows rather than the ragged index that which() produced (which crashed
## with "incorrect number of dimensions" / silently mis-joined on
## duplicates).
.get_daily_value <- function(date, daily_value) {
  dv_date <- if (is.data.frame(daily_value))
    daily_value$date else daily_value[, "date"]
  idx <- match(as.integer(date), as.integer(dv_date))
  daily_value[idx, -1]
}
## goodness_of_fit() lives in its own file R/goodness_of_fit.R (was duplicated
## here).
