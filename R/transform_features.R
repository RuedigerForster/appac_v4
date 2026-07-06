## Dispersion-normalise the areas and split them into three components by PCA:
## correlated (PC1 — the dominant common drift + pressure signal), uncorrelated
## (PC2+ — abrupt / per-peak structure) and noise.  Each peak is referenced to
## its
## whole-series mean — the centre debias_ct() de-biases around (see below).
##
## Arguments:
##   d  — per-sample date vectors.
##   Y  — per-sample raw-area matrices (cols = peaks).
##   ct — optional centres per sample (NULL = whole-series mean per peak).
##   sc — optional scales per sample (NULL = per-peak sd).
## Value:
##   list(correlated, uncorrelated, noise); each a per-sample list with a
##   leading
##   date column and attr "mean" (= ct) / "scale".
.transform_features <- function(d, Y, ct = NULL, sc = NULL) {

  samples <- names(d)

  ##*****************************************************

  if (is.null(sc)) {
    sc <- lapply(Y, function(x) apply(x, 2, sd))
  }
  if (is.null(ct)) {
    ## Reference each cylinder to its WHOLE-SERIES mean.  debias_ct()
    ## de-biases the centre by sweeping +/-1% around it, so the centre must
    ## start
    ## near the chi-square optimum — the overall mean does, a day-1/ep0 mean
    ## does
    ## not (it sits >1% off, which puts the optimum outside the sweep and
    ## crashes
    ## the uniroot).  The starting reference does not change the final de-biased
    ## centre anyway, since debias_ct() re-optimises it.
    ct <- lapply(Y, function(x) apply(x, 2, mean))
    names(ct) <- names(Y)
  }

  ## linearized sc ~ 0 + ct is used for dispersion normalization
  x <- unlist(ct)
  x2 <- x^2
  y <- unlist(sc)
  ex <- coef(lm(y ~ 0 + x + x2))

  ##*****************************************************

  Y_scld <- lapply(Y, function(x) scale(x))

  ##*****************************************************

  ## dispersion normalization
  sc_new <- lapply(ct, function(x)
    x * ex[1] / (x * ex[1] + x^2 * ex[2])
  )
  Y_scld <- lapply(seq_along(samples), function(i)
    sweep(Y_scld[[i]], 2, sc_new[[i]], "*")
  )

  ##*****************************************************

  ## 1st Principal Component Analysis
  ## divides the data into correlated (through the sample = PC1) and
  ## uncorrelated through the sample. although correlated through the
  ## instrument = PC2) arrays and completely uncorrelated noise (= PC3 and
  ## higher)
  # browser()
  pca <- lapply(Y_scld, function(x)
    prcomp(x, scale = FALSE, center = FALSE)
  )

  ## PC1 includes everything correlated
  Y_cor <- lapply(pca, function(z)
    z$x[, 1] %*% t(z$rotation[, 1])
  )

  ## PC2, includes the bias, omits the noise
  Y_unc <- lapply(seq_along(pca), function(i)
    pca[[i]]$x[, 2] %*% t(pca[[i]]$rotation[, 2])
  )

  ## PC3 and higher contain the noise
  Y_noiz <- lapply(seq_along(pca), function(i)
    pca[[i]]$x[, -c(1, 2)] %*% t(pca[[i]]$rotation[, -c(1, 2)])
  )

  ##*****************************************************

  ## reverse scaling
  Y_cor <- lapply(seq_along(Y_cor), function(i)
    sweep(Y_cor[[i]], 2, sc[[i]], "*")
  )
  Y_cor <- lapply(seq_along(Y_cor), function(i)
    sweep(Y_cor[[i]], 2, ct[[i]], "+")
  )
  Y_cor <- lapply(seq_along(Y_cor), function(i)
    cbind(date = d[[i]], Y_cor[[i]])
  )
  names(Y_cor) <- samples
  # browser()
  ##*****************************************************

  Y_unc <- lapply(seq_along(Y_unc), function(i)
    sweep(Y_unc[[i]], 2, sc[[i]], "*")
  )
  Y_unc <- lapply(seq_along(Y_unc), function(i)
    sweep(Y_unc[[i]], 2, ct[[i]], "+")
  )
  Y_unc <- lapply(seq_along(Y_unc), function(i)
    cbind(date = d[[i]], Y_unc[[i]])
  )
  names(Y_unc) <- samples

  ##*****************************************************

  Y_noiz <- lapply(seq_along(Y_noiz), function(i)
    sweep(Y_noiz[[i]], 2, sc[[i]], "*")
  )
  Y_noiz <- lapply(seq_along(Y_noiz), function(i)
    sweep(Y_noiz[[i]], 2, ct[[i]], "+")
  )
  Y_noiz <- lapply(seq_along(Y_noiz), function(i)
    cbind(date = d[[i]], Y_noiz[[i]])
  )
  names(Y_noiz) <- samples

  std_dev <- lapply(Y_noiz, function(x)
    apply(x[, -1], 2, sd)
  )

  ##*****************************************************

  for (i in seq_along(samples)) {
   attr(Y_cor[[i]], "mean") <- ct[[i]]
   attr(Y_cor[[i]], "scale") <- sc_new[[i]]
   attr(Y_unc[[i]], "mean") <- ct[[i]]
   attr(Y_noiz[[i]], "mean") <- ct[[i]]
   attr(Y_noiz[[i]], "standard.deviation") <- std_dev[[i]]
  }

  out <- list(
    correlated = Y_cor,
    uncorrelated = Y_unc,
    noise = Y_noiz
  )
  return(out)
}
