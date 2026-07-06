## changepoints.R - structural-break detection (strucchange, Zeileis et al.)
##
## Two complementary detectors, both a change in the MEAN of a per-sample daily
## signal fitted by strucchange (deterministic; no MCMC):
##   * get_changepoints()          - LEVEL breaks: the signal is the
##   daily-averaged
##     PC2 score of PCA(1) on the (per-peak standardised) areas.  PC1 is the
##     dominant smooth drift, PC2 the second mode that carries abrupt level
##     shifts.
##   * get_variance_changepoints() - VARIANCE breaks: the signal is the daily
##     mean of the per-injection noise energy (sum of squared residuals after
##     removing the leading PC).  A change in its mean is a change in the
##     measurement variance.
## An OLS-MOSUM empirical fluctuation test (strucchange::efp / sctest) gates
## significance; strucchange::breakpoints() dates the breaks (Bai-Perron dynamic
## programming, BIC).  Long daily series are coarsened to `max_grid` points
## first,
## so dating stays fast.  Cross-cylinder breakpoints = the union of the
## per-cylinder
## breaks, merged when closer than `merge_within` days.

## Detect changes in the mean of a per-date signal with strucchange.  Returns
## the
## break dates for ONE cylinder (empty if the series is too short, the MOSUM
## test
## is not significant at `alpha`, or no break improves the BIC).
.detect_mean_changes <- function(date, y, alpha = 0.05, h = 0.15,
                                 max_grid = 250L, min_seg = 30L) {
  ord <- order(date); date <- as.integer(date)[ord]; y <- y[ord]
  keep <- is.finite(y); date <- date[keep]; y <- y[keep]
  n <- length(y)
  if (n < 2L * min_seg) return(data.table::as.IDate(integer(0)))

  ## coarsen a long daily axis onto <= max_grid bins (keep the last date of
  ## each)
  if (n > max_grid) {
    g    <- ceiling(seq_len(n) / (n / max_grid))
    date <- as.integer(tapply(date, g, function(z) z[length(z)]))
    y    <- as.numeric(tapply(y, g, mean))
    n    <- length(y)
  }
  y <- as.numeric(scale(y))
  if (!all(is.finite(y))) return(data.table::as.IDate(integer(0)))

  ## MOSUM significance gate: is there a mean shift at all?
  p <- tryCatch(
    strucchange::sctest(
      strucchange::efp(y ~ 1, type = "OLS-MOSUM", h = h))$p.value,
    error = function(e) NA_real_)
  if (!is.finite(p) || p > alpha) return(data.table::as.IDate(integer(0)))

  ## date the breaks (BIC-optimal number and location)
  bp <- tryCatch(strucchange::breakpoints(y ~ 1, h = h)$breakpoints,
                 error = function(e) NA_integer_)
  if (length(bp) == 0L || all(is.na(bp)))
    return(data.table::as.IDate(integer(0)))
  data.table::as.IDate(sort(unique(date[bp])))
}

## Standardise each peak (unit-variance, zero-variance guarded to 0) and run
## PCA(1).  Returns the prcomp object, or NULL if there is nothing to decompose
## (fewer than 2 peaks) or the series is too short (< 2 * min_seg injections).
.standardise_pca <- function(Y, min_seg) {
  Y <- as.matrix(Y)
  if (ncol(Y) < 2L || nrow(Y) < 2L * min_seg) return(NULL)
  y <- apply(Y, 2, scale); y[!is.finite(y)] <- 0        # standardise each peak
  prcomp(y, center = FALSE, scale. = FALSE)
}

## Daily-average a per-injection signal onto its date axis, then look for
## changes
## in its mean (the shared epilogue of both detectors).
.detect_daily <- function(d, sig, alpha, h, max_grid, min_seg) {
  dr <- aggregate(
    sig ~ date,
    data = data.frame(date = as.integer(d), sig = sig), FUN = mean)
  .detect_mean_changes(dr$date, dr$sig, alpha, h, max_grid, min_seg)
}

## Per-cylinder LEVEL signal: daily-averaged PC2 score.
.get_chp <- function(d, Y, alpha = 0.05, h = 0.15, max_grid = 250L,
                     min_seg = 30L) {
  pca <- .standardise_pca(Y, min_seg)
  if (is.null(pca) || ncol(pca$x) < 2L) return(data.table::as.IDate(integer(0)))
  .detect_daily(d, pca$x[, 2], alpha, h, max_grid, min_seg)
}

## Per-cylinder VARIANCE signal: daily mean of the per-injection noise energy
## (sum of squared residuals after dropping the leading `drop` PC(s)).
.get_var_chp <- function(d, Y, drop = 1L, alpha = 0.05, h = 0.15,
                         max_grid = 250L, min_seg = 30L) {
  pca <- .standardise_pca(Y, min_seg)
  if (is.null(pca) || ncol(pca$x) < drop + 2L)
    return(data.table::as.IDate(integer(0)))
  keep <- (drop + 1L):ncol(pca$x)   # drop leading PC(s) = signal
  ## per-injection noise energy = ||residual row||^2 = row sum of squares of the
  ## KEPT PC scores (orthonormal rotation), no residual-matrix reconstruction.
  energy <- rowSums(pca$x[, keep, drop = FALSE]^2)
  .detect_daily(d, energy, alpha, h, max_grid, min_seg)
}

## Union the per-cylinder breaks and merge near-coincident ones.
.merge_breaks <- function(per_sample, merge_within) {
  bp <- sort(unique(do.call(c, per_sample)))
  if (length(bp) == 0L) return(data.table::as.IDate(integer(0)))
  data.table::as.IDate(bp[c(TRUE, diff(as.integer(bp)) > merge_within)])
}

#' Detect episode level breakpoints across cylinders
#'
#' Detects abrupt level shifts ("episode" boundaries) in the area series.  Per
#' cylinder the change signal is the second principal-component score of a
#' one-component PCA on the per-peak standardised areas; a structural-break
#' model
#' (\code{\link[strucchange]{breakpoints}}, gated by an OLS-MOSUM fluctuation
#' test,
#' \code{\link[strucchange]{efp}} / \code{\link[strucchange]{sctest}}) is fitted
#' to
#' the daily-averaged PC2 series and the BIC-optimal breakpoints are kept.  The
#' cross-cylinder breakpoints are the union of the per-cylinder breaks, merged
#' when
#' closer than \code{merge_within} days.  Detection is deterministic.
#'
#' @param samples Named list, each element with \code{$date} (vector) and
#'   \code{$raw.area} (peak matrix) -- i.e. the \code{samples} slot of an
#'   \code{"Appac"} object.
#' @param alpha Significance level of the OLS-MOSUM test that gates whether a
#'   cylinder has any break.
#' @param h Minimum segment width, as a fraction of the (coarsened) series
#' length,
#'   for the MOSUM bandwidth and the breakpoint search.
#' @param max_grid Long daily series are coarsened onto at most this many bins
#'   before dating, so the search stays fast.
#' @param min_seg Minimum daily-series length (days) required to attempt
#' detection.
#' @param merge_within Merge breakpoints from different cylinders that fall
#' within
#'   this many days of each other.
#' @return An \code{\link[data.table]{IDate}} vector of breakpoint dates.
#' @seealso \code{\link{get_variance_changepoints}}, \code{\link{appac}},
#'   \code{\link{plot_area_date}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' get_changepoints(fit@samples)   # episode breakpoint dates
#' }
#' @export
get_changepoints <- function(samples, alpha = 0.05, h = 0.15, max_grid = 250L,
                             min_seg = 30L, merge_within = 3) {
  chp <- lapply(samples, function(s)
    .get_chp(s$date, s$raw.area, alpha, h, max_grid, min_seg))
  .merge_breaks(chp, merge_within)
}

#' Detect variance breakpoints across cylinders
#'
#' Detects abrupt changes in the measurement \emph{variance} (precision), the
#' second-moment counterpart of \code{\link{get_changepoints}}.  Per cylinder
#' the
#' signal is the daily mean of the per-injection noise energy -- the sum of
#' squared
#' residuals after removing the leading \code{drop} principal component(s) (the
#' level / pressure / drift) -- and a change in its mean is a change in the
#' noise
#' variance.  The same strucchange machinery (OLS-MOSUM gate +
#' \code{\link[strucchange]{breakpoints}}) dates the breaks; results are merged
#' and
#' returned as for \code{\link{get_changepoints}}.
#'
#' A pure mean (level) shift leaves the variance unchanged, so this detector
#' does
#' \emph{not} flag the level breaks that \code{\link{get_changepoints}} finds --
#' the
#' two are complementary.  Small variance steps sit near the detection floor: a
#' change is only found when it is significant at \code{alpha}.
#'
#' @param samples Named list, each element with \code{$date} and
#' \code{$raw.area},
#'   i.e. the \code{samples} slot of an \code{"Appac"} object.
#' @param drop Number of leading principal components to remove (the level /
#'   pressure / drift) before forming the noise energy.
#' @param alpha,h,max_grid,min_seg,merge_within As in
#' \code{\link{get_changepoints}}.
#' @return An \code{\link[data.table]{IDate}} vector of variance-breakpoint
#' dates.
#' @seealso \code{\link{get_changepoints}}, \code{\link{appac}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' get_variance_changepoints(fit@samples)   # precision-change dates
#' }
#' @export
get_variance_changepoints <- function(samples, drop = 1L, alpha = 0.05,
                                       h = 0.15, max_grid = 250L,
                                       min_seg = 30L, merge_within = 3) {
  chp <- lapply(samples, function(s)
    .get_var_chp(s$date, s$raw.area, drop, alpha, h, max_grid, min_seg))
  .merge_breaks(chp, merge_within)
}
