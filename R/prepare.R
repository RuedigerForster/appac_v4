
## Impute the missing cells of a (dates x peaks) area matrix by low-rank
## reconstruction (svdImpute / EM), so the downstream loess + PCA see a complete
## matrix.  Within a sample the peaks share the common drift / pressure / daily
## factor, so a rank-`ncp` reconstruction fills a gap from the OTHER peaks of
## the
## same injection.  Complete matrices are returned untouched; non-missing cells
## are never altered.  Shares the EM loop (.svd_impute) with .na_pca(), but
## centres
## the columns and returns the completed matrix rather than PCA scores.
.impute_na <- function(x, ncp = 2L, maxiter = 200L, tol = 1e-7) {
  na <- is.na(x)
  if (!any(na)) return(x)
  mu <- colMeans(x, na.rm = TRUE)
  s  <- apply(x, 2, function(col) stats::sd(col, na.rm = TRUE))
  s[!is.finite(s) | s == 0] <- 1
  ncp <- min(ncp, ncol(x) - 1L, nrow(x) - 1L)
  if (ncp < 1L) {   # too small to decompose: mean-fill
    for (j in seq_len(ncol(x))) x[na[, j], j] <- mu[j]
    return(x)
  }
  # centre + unit variance
  xs <- .svd_impute(sweep(sweep(x, 2, mu, "-"), 2, s, "/"), ncp, maxiter, tol)
  filled <- sweep(sweep(xs, 2, s, "*"), 2, mu, "+") # undo the scaling
  x[na] <- filled[na]
  x
}

## Pivot the long-format, column-checked data into per-sample matrices and wrap
## them in an "Appac" object ready for .calculate().  Peaks missing in a sample
## are handled; areas are stored per (sample x peak) with their injection dates.
##
## Arguments:
##   df         — long-format data with canonical column names (see check_cols).
##   P_ref      — reference value for the covariate(s) (e.g. reference
##   pressure).
##   covariates — name(s) of the covariate column(s) to carry through.
## Value:
##   An "Appac" S4 object whose @samples[[s]] holds $date, $raw.area (peak
##   matrix)
##   and an empty $corrected.area to be filled by .calculate().
.prepare <- function(df, P_ref, covariates = "pressure") {

  Correction <- methods::new("Correction")
  Compensation <- methods::new("Compensation")

  samples <- list(

  )

  spls <- unique(df$Sample_Name)
  cmps <- unique(df$Peak_Name)

  for (smp in spls) {
    #----------------------------------------------------------
    # some data wrangling
    #----------------------------------------------------------

    df_ <- df %>% dplyr::filter(df$Sample_Name == smp)
    idx <- lapply(cmps, function(x) df_$Peak_Name == x)

    #----------------------------------------------------------
    # Input data:
    #
    # P: ambient pressure vector (averaged)
    # X: date vector (ambiguity checked)
    # Y: data frame of raw areas; colnames are peak names
    # also take care of peaks which are missing in a sample
    #----------------------------------------------------------

    ## covariate: pressure
    Correction@covariates <- covariates
    for (cv in covariates) {
      Correction@coefficients[[cv]] <- NA
      Correction@reference.values[[cv]] <- P_ref
    }

    P <- rowMeans(do.call(cbind, lapply(idx, function(x) {
      df_$Air_Pressure[x]
    })), na.rm = TRUE)

    X <- do.call(cbind, lapply(idx, function(x) df_$Injection_Date[x]))
    if (!all(sapply(2:ncol(X), function(x) identical(X[, 1], X[, x])))) {
      stop("Inconsistent date vectors.")
    } else {
      X <- as.integer(X[, 1])
    }
    unique_X <- sort(unique(X))
    l <- length(unique_X)
    Compensation@date <- unique_X
    Compensation@bias <- rep(0, l)
    Compensation@trend <- rep(0, l)
    Compensation@correlated.features <- rep(0, l)

    Y <- do.call(cbind, lapply(idx, function(x) df_$Raw_Area[x]))
    missing.peak <- colSums(is.na(Y)) > 0.9 * nrow(Y)
    if (any(missing.peak)) {
      cmpl <- cmps[-which(missing.peak)]
      Y <- Y[, -which(missing.peak), drop = FALSE]
    } else {
      cmpl <- cmps
    }
    colnames(Y) <- cmpl
    rsd <- sapply(
      seq_len(ncol(Y)),
      function(x) stats::sd(Y[, x], na.rm = TRUE) / mean(Y[, x], na.rm = TRUE)
    )
    cutoff <- 0.025
    if (any(rsd > cutoff)) {
      warning(
        "The noise in the input data of peak(s): ",
        paste0("'", colnames(Y)[rsd > cutoff], sep = "'"),
        " in sample '", smp,
        "' exceeds the allowed noise level of ",
        sprintf("%.1f", cutoff * 100), "%."
      )
    }

    # NA in the data cause trouble. Peaks with more than `max_na_fraction`
    # NA (default 30%) are eliminated; the rest are imputed below.
    idx <- sapply(
      seq_len(ncol(Y)),
      function(y) sum(is.na(Y[, y])) / nrow(Y)
    )
    idx <- idx <= max_na_fraction
    Y <- Y[, idx, drop = FALSE]
    cmpl <- cmpl[idx]

    ## impute the surviving (<= 30% per peak) NA cells so the loess de-drift and
    ## the PCA in .transform_features() receive a complete matrix; complete
    ## input
    ## is passed through untouched.
    if (anyNA(Y)) Y <- .impute_na(Y, ncp = 2L)

    c <- ncol(Y)
    r <- nrow(Y)


    samples[[smp]] <- list(
      sample.name = smp,
      date = data.table::as.IDate(X),
      # pressure = P,
      raw.area = Y,
      corrected.area = matrix(NA_real_, ncol = c, nrow = r)
    )

    Correction@samples[[smp]] <- list(
      sample.name = smp,
      date = data.table::as.IDate(X),
      pressure = P,
      raw.area = Y,
      corrected.area = matrix(NA, ncol = ncol(Y), nrow = nrow(Y))
    )

  }

  Appac <- methods::new("Appac", samples = samples, trend = Compensation,
                        correction = Correction)
  return(Appac)
}
