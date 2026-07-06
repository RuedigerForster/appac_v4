#' Goodness of fit of the correction
#'
#' Tests, per peak, the null hypothesis that the corrected areas are constant
#' (equal to the centre) up to measurement noise.  Residual structure left
#' behind
#' by the correction (drift, steps) inflates the statistic.
#'
#' The chi-square sums the squared residuals about the centre, each scaled by
#' the
#' measurement-noise variance of its peak.  That variance is estimated from
#' successive injections -- the von Neumann lag-1 estimator
#' \eqn{\sigma^2 = mean(d^2)/2} with \eqn{d} the date-ordered first differences
#' --
#' so it captures the irreducible short-term noise and is independent of any
#' slow
#' residual drift.  Scaling by a noise estimate (rather than by the centre, as a
#' Pearson form would, which assumes a Poisson variance) makes the statistic
#' correct for continuous peak areas and dimensionless.
#'
#' Use the \strong{reduced} chi-square (chi-square / dof): about 1 means the
#' corrected areas are down to the short-term noise floor (a good correction),
#' above 1 flags residual structure.  The p-value is returned too but is of
#' little use here: with thousands of injections the distribution is so tight
#' that any reduced chi-square a touch above 1 gives p near 0.
#'
#' @param Appac A corrected \code{"Appac"} object (e.g. from a second
#'   \code{\link{appac}} pass).
#' @return A named list (per sample) of data frames, one row per peak, with
#'   columns \code{reduced.chisq}, \code{chisq}, \code{dof} and \code{p.value}.
#' @seealso \code{\link{appac}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' goodness_of_fit(fit)[[1]]   # per-peak reduced chi-square, sample 1
#' }
#' @export
goodness_of_fit <- function(Appac) {

  samples <- names(Appac@samples)

  out <- lapply(seq_along(samples), function(i) {
    X   <- Appac@samples[[i]]$corrected.area
    ctr <- Appac@trend@center[[i]]
    ord <- order(as.integer(Appac@samples[[i]]$date))   # noise needs time order
    sigma2 <- vapply(seq_len(ncol(X)),
                     function(j) mean(diff(X[ord, j])^2, na.rm = TRUE) / 2,
                     numeric(1))
    chisq <- colSums(sweep(sweep(X, 2, ctr, "-")^2, 2, sigma2, "/"),
                     na.rm = TRUE)
    dof   <- length(Appac@samples[[i]]$date) - 1L
    peaks <- colnames(X)
    if (is.null(peaks)) peaks <- paste0("peak", seq_along(chisq))
    data.frame(reduced.chisq = as.numeric(chisq) / dof,
               chisq         = as.numeric(chisq),
               dof           = dof,
               p.value       = stats::pchisq(chisq, dof, lower.tail = FALSE),
               row.names     = peaks)
  })
  names(out) <- samples

  return(out)
}
