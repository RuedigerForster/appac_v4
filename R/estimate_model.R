## NA-tolerant PCA in base R only (no dependency).  Replaces
## pcaMethods::pca(method="bpca") (Bioconductor).  Like bpca, it imputes the
## missing daily-factor cells by RECONSTRUCTION rather than skipping them (as
## NIPALS does): unit-variance scale (center=FALSE, matching the old
## scale="uv"),
## seed the gaps at the mean, then iterate a rank-ncp SVD reconstruction into
## the
## missing cells until stable (the svdImpute / EM scheme).  Returns the pieces
## the callers use: scores, loadings, per-component variance fraction R2.
.na_pca <- function(x, ncp = 2, maxiter = 200L, tol = 1e-7) {
  s  <- apply(x, 2, function(col) stats::sd(col, na.rm = TRUE))
  # unit-variance, no centering
  xs <- .svd_impute(sweep(x, 2, s, "/"), ncp, maxiter, tol)
  sv <- svd(xs)
  list(scores   = sv$u[, 1:ncp, drop = FALSE] %*% diag(sv$d[1:ncp], ncp),
       loadings = sv$v[, 1:ncp, drop = FALSE],
       R2       = sv$d[1:ncp]^2 / sum(sv$d^2))
}

## Estimate the per-sample slow drift/bias from the uncorrelated (PC2+) daily
## factors: a cross-sample PCA whose first two components give the trend and
## bias.
##
## Arguments:
##   uncorrelated_data — per-sample uncorrelated component from
##   transform_features.
## Value:
##   list(date, trend, bias, scaling) — daily trend/bias series over `date`,
##   plus a
##   per-sample scaling.
.estimate_uncorrelated_model <- function(uncorrelated_data) {

  ##*****************************************************
  ## capture the names
  ##*****************************************************

  samples <- names(uncorrelated_data)

  ##*****************************************************
  ## scaling of measured data
  ##*****************************************************

  y_scld <- lapply(seq_along(samples), function(i)
    sweep(
      subset(uncorrelated_data[[i]], select = -date),
      2,
      attr(uncorrelated_data[[i]], "mean"),
      "/"
    )
  )

  ##*****************************************************

  sc <- lapply(seq_along(samples), function(i)
    apply(y_scld[[i]], 2, sd)
  )
  names(sc) <- samples
  ct <- lapply(seq_along(samples), function(i)
    apply(y_scld[[i]], 2, mean)
  )
  names(ct) <- samples

  y_scld <- lapply(seq_along(samples), function(i)
    scale(y_scld[[i]], center = ct[[i]], scale = sc[[i]])
  )
  names(y_scld) <- samples

  ##*****************************************************
  ## take care of the sign in scaling, it may change
  ##*****************************************************

  sgn <- lapply(y_scld, function(x)
    cov(x)[1, ]
  )
  ## just in case the covariance matrix is not
  ## a matrix of ones: use sign function
  sgn <- lapply(sgn, function(x)
    if(x[1] > 0) sign(x) else -sign(x)
  )
  sc <- lapply(seq_along(samples), function(i)
    sc[[i]] * sgn[[i]]
  )

  ##*****************************************************

  daily_factor <- lapply(seq_along(samples), function(i)
    data.frame(
      date = as.IDate(uncorrelated_data[[i]][, 1]),
      factor = y_scld[[i]][, 1]
    )
  )
  daily_factor <- lapply(seq_along(samples), function(i)
    .summarize_factors(daily_factor[[i]])
  )
  daily_factor <- suppressWarnings(
    purrr::reduce(daily_factor, merge, by = "date", all = TRUE, no.dups = FALSE)
  )
  daily_factor[, "date"] <- as.IDate(daily_factor[, "date"])
  colnames(daily_factor) <- c("date", samples)

  ##*****************************************************

  ## PCA of uncorrelated data
  ## NAs make trouble
  n_na <- sapply(2:ncol(daily_factor), function(i)
    sum(!is.na(daily_factor[, i])) / nrow(daily_factor)
  )
  names(n_na) <- samples
  usable <- c(date = FALSE, n_na > 0.1)
  if (sum(usable) < 3)
    stop("Insufficient number of experiments to determine drift and bias.")
  sample <- as.matrix(daily_factor[, usable])
  ## remove rows where all cols are NAs
  idx <- rowSums(is.na(sample)) != ncol(sample)
  sample <- as.matrix(sample[idx, ])
  pca <- .na_pca(sample, ncp = 2)

  PCs <- data.frame(
    date = daily_factor[idx, 1],
    trend = (pca$scores[, 1] %*% t(pca$loadings[, 1]))[, 1],
    bias = (pca$scores[, 2] %*% t(pca$loadings[, 2]))[, 1]
  )
  attr(PCs$trend, "importance") <- pca$R2[1]
  attr(PCs$bias, "importance") <- pca$R2[2]

  ##*****************************************************

  out <- list(
    date = PCs$date,
    trend = PCs$trend,
    bias = PCs$bias,
    scaling = sc
  )

  return(out)
}


## Estimate the cross-sample correlated daily factor (PC1 of the concurrent
## correlated components) — the common day-to-day drift after kappa correction.
##
## Arguments:
##   correlated_data — per-sample correlated component from transform_features.
## Value:
##   list(daily.factor, rolling.variance, scaling) — the daily factor (date +
##   value),
##   its rolling variance, and a per-sample scaling.
.estimate_correlated_model <- function(correlated_data) {

  ##*****************************************************
  ## capture the names
  ##*****************************************************

  samples <- names(correlated_data)

  ##*****************************************************
  ## scaling of measured data
  ##*****************************************************

  y_scld <- lapply(seq_along(samples), function(i)
    sweep(
      subset(correlated_data[[i]], select = -date),
      2,
      attr(correlated_data[[i]], "mean"),
      "/"
    )
  )

  ##*****************************************************

  sc <- lapply(seq_along(samples), function(i)
    apply(y_scld[[i]], 2, sd)
  )
  names(sc) <- samples
  ct <- lapply(seq_along(samples), function(i)
    apply(y_scld[[i]], 2, mean)
  )
  names(ct) <- samples

  y_scld <- lapply(seq_along(samples), function(i)
    scale(y_scld[[i]], center = ct[[i]], scale = sc[[i]])
  )
  names(y_scld) <- samples

  ##*****************************************************

  daily_factor <- lapply(seq_along(samples), function(i)
    data.frame(
      date = as.IDate(correlated_data[[i]][, 1]),
      factor = y_scld[[i]][, 1]
    )
  )
  daily_factor <- lapply(seq_along(samples), function(i)
    .summarize_factors(daily_factor[[i]])
  )
  daily_factor <- suppressWarnings(
    purrr::reduce(daily_factor, merge, by = "date", all = TRUE, no.dups = FALSE)
  )
  daily_factor[, "date"] <- as.IDate(daily_factor[, "date"])
  colnames(daily_factor) <- c("date", samples)

  ##*****************************************************

  ## PCA of correlated data
  ## NAs make trouble
  n_na <- sapply(2:ncol(daily_factor), function(i)
    sum(!is.na(daily_factor[, i])) / nrow(daily_factor)
  )
  names(n_na) <- samples
  usable <- c(date = FALSE, n_na > 0.7)
  if (sum(usable) < 3) stop("Insufficient number of concurrrent experiments.")
  sample <- daily_factor[, usable]
  idx <- rowSums(is.na(sample)) != ncol(sample)
  sample <- as.matrix(sample[idx, ])
  pca <- .na_pca(sample, ncp = 2)

  pca1 <- data.frame(
    date = daily_factor[idx, 1],
    value = (pca$scores[, 1] %*% t(pca$loadings[, 1]))[, 1]
  )
  attr(pca1, "importance") <- pca$R2[1]

  pca2 = data.frame(
    date = daily_factor[idx, 1],
    value = (pca$scores[, 2] %*% t(pca$loadings[, 2]))[, 1]
  )
  attr(pca2, "importance") <- pca$R2[2]
  roll_var <- data.frame(
    date = pca2[, 1],
    time.series = pca2[, 2],
    variance = kza::rlv(pca2[, 2], 7)
  )

  ##*****************************************************

  out <- list(
    daily.factor = pca1,
    rolling.variance = roll_var,
    scaling = sc
  )

  return(out)
}


