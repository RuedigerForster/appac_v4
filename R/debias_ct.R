#' De-bias the per-peak centres
#'
#' Refines ("de-biases") the per-peak centres by minimising the chi-square of
#' the
#' corrected-area residuals about the centre.  Sweeps a scale factor \code{cf}
#' around the current centre (the whole-series mean), re-runs
#' \code{\link{appac}}
#' at each \code{cf}, fits the per-peak chi-square as a parabola in \code{cf}
#' and
#' takes its (closed-form) minimum.  Trades a little variance for reduced bias.
#'
#' @param Appac A fitted \code{"Appac"} object (from a first \code{\link{appac}}
#'   pass).
#' @param data The same column-checked data passed to that \code{appac()} call.
#' @param P_ref The same reference pressure (hPa).
#' @param npt Number of sweep points for \code{cf} in \code{[0.99, 1.01]}.
#' @param quiet Suppress the progress dots.
#' @return A list of de-biased centres (one numeric vector per sample), to feed
#'   back as \code{appac(..., ct = <this>)}.
#' @seealso \code{\link{appac}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit1 <- appac(dat, P_ref = P_ref)
#' ct   <- debias_ct(fit1, data = dat, P_ref = P_ref, quiet = TRUE)
#' fit  <- appac(dat, ct = ct, P_ref = P_ref)   # de-biased second pass
#' }
#' @export
debias_ct <- function(Appac, data, P_ref, npt = 20, quiet = FALSE) {

  samples <- names(Appac@samples)
  cfl <- seq(from = 0.99, to = 1.01, length.out = npt)
  ct  <- Appac@trend@center

  ## sample chi^2(cf): re-run appac() at each centre scaling cf, drawing a
  ## 21-wide progress bar (<=====    >) that redraws in place after each sweep.
  chisqu <- list()
  for (i in seq_along(cfl)) {
    sct <- lapply(ct, function(x) cfl[i] * x)
    XX  <- appac(data = data, ct = sct, P_ref = P_ref)
    chisqu[[i]] <- lapply(seq_along(samples), function(j)
      colSums(
        sweep(XX@samples[[j]]$corrected.area, 2, XX@trend@center[[j]], "-")^2)
    )
    if (!quiet) {
      f <- round(i / npt * 21)
      cat(sprintf("\rChi-square minimisation: <%s%s>",
                  strrep("=", f), strrep(" ", 21L - f)))
      utils::flush.console()
    }
  }
  if (!quiet) cat("\n")

  ## per peak, fit chi^2(cf) = b[1] + b[2]*cf + b[3]*cf^2 and take the minimum.
  ## The minimum is the parabola vertex cf* = -b[2] / (2*b[3]) in closed form;
  ## b[3] > 0 confirms it is a minimum.  (This replaces a uniroot() on the
  ## derivative, which threw "values at end points not of opposite sign"
  ## whenever the optimum fell outside the [0.99, 1.01] sweep — exactly the
  ## large-bias case the function exists for.  The closed form has no bracket
  ## to escape.)
  bias <- list()
  for (smp in seq_along(samples)) {
    dd <- sapply(chisqu, function(x) x[[smp]])      # rows = peaks, cols = cf
    if (all(is.na(dd))) { bias[[smp]] <- numeric(0); next }
    b <- sapply(seq_len(nrow(dd)),
                function(k) coef(lm(dd[k, ] ~ cfl + I(cfl^2))))
    vertex <- -b[2, ] / (2 * b[3, ])
    bad <- !(b[3, ] > 0) | !is.finite(vertex)   # non-convex / degenerate fit
    if (any(bad)) {
      warning(sprintf(paste0("debias_ct(): %d peak(s) in '%s' have no ",
                             "chi-square minimum; centre left unbiased."),
                      sum(bad), samples[smp]))
      vertex[bad] <- 1                          # leave those centres unchanged
    }
    oor <- !bad & (vertex < min(cfl) | vertex > max(cfl))
    if (any(oor))
      warning(sprintf(paste0("debias_ct(): %d peak(s) in '%s' optimise ",
                             "outside the sweep [%.3f, %.3f] (extrapolated) ",
                             "- widen the cf range."),
                      sum(oor), samples[smp], min(cfl), max(cfl)))
    bias[[smp]] <- vertex
  }
  names(bias) <- samples

  ct <- lapply(seq_along(samples), function(i)
    if (length(bias[[i]]) == 0) ct[[i]] else ct[[i]] * bias[[i]]
  )
  names(ct) <- samples
  return(ct)
}
