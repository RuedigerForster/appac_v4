
## Assemble and apply the full APPAC correction.  Kappa is estimated from a
## DRIFT-REDUCED signal (the areas are de-drifted with a robust loess low-pass
## on
## a COPY, so kappa is not confounded with the slow drift); the
## daily-factor/trend
## drift model is then estimated on the kappa-corrected data, and the combined
## multiplier (kappa + drift) is divided out of the raw areas.
##
## Arguments:
##   Appac — an "Appac" object from .prepare() (raw areas, dates, covariates).
##   ct    — optional fixed centres (passed through to transform_features).
## Value:
##   The Appac object with @samples[[i]]$corrected.area,
##   @correction@coefficients
##   (kappa) and @trend (the drift model) populated.
.calculate <- function(Appac, ct) {

  ## .get_daily_value() is defined once in helper.R (matrix/data.frame
  ## compatible,
  ## match()-based so duplicate/missing dates don't crash); was a local copy
  ## here.

  ##*****************************************************
  ## prepare the data
  ##*****************************************************

  d <-  covar <-  Y <- list()
  samples <- sapply(Appac@samples, function(x) x$sample.name)
  for (sample_name in samples) {
    ## Y: list of raw areas
    Y[[sample_name]] <- Appac@samples[[sample_name]]$raw.area
    ## p: list of pressure vectors (P in mbar), referenced
    for (cv in Appac@correction@covariates) {
      covar[[cv]][[sample_name]] <- unname(unlist(
        unlist(Appac@correction@samples[[sample_name]][cv]) -
        unlist(Appac@correction@reference.values[cv])
      ))
      # Appac@correction@samples[[sample_name]][cv] <-
      #   covar[[cv]][[sample_name]]
    }
    ## d: list of dates vectors
    d[[sample_name]] <- Appac@samples[[sample_name]]$date
  }

  ##*****************************************************
  ## some consistency checks
  ##*****************************************************

  ## these critical errors should never occur
  if (!identical(sapply(d, length), sapply(Y, nrow)))
    stop("Not compliant lengths of dates vectors and area matrices.")
  for (cv in Appac@correction@covariates) {
    if (!identical(sapply(covar[[cv]], length), sapply(Y, nrow)))
      stop("Not compliant lengths of covariate vectors and area matrices.")
  }

  ##*****************************************************
  ## here's the actual calculation
  ##*****************************************************

    ## DRIFT-THEN-KAPPA: de-drift with a robust temporal low-pass (loess,
    ## family = "symmetric") on a COPY of the areas, then fit kappa; original Y
    ## drives the correction below.  Two rejected alternatives:
    ##  - appac's daily-factor model absorbs the residual pressure correlation
    ##    and strips it from the kappa fit (kappa -> ~0);
    ##  - per-episode-mean de-bias over-corrects (kappa -> -1.2e-3) because it
    ##    removes the pressure-correlated level differences BETWEEN episodes,
    ##    not
    ##    just the slow drift.
    ## The continuous low-pass removes only the slow drift, preserving the
    ## within-day pressure response, and recovers the validated kappa.
    Y_dd <- lapply(seq_along(Y), function(i) {
      ym <- Y[[i]]; dt <- as.numeric(d[[i]])
      for (j in seq_len(ncol(ym))) {
        tr <- fitted(loess(ym[, j] ~ dt, span = 0.1, degree = 1,
                           family = "symmetric"))
        ym[, j] <- ym[, j] / (tr / mean(tr))
      }
      ym
    })
    Y_div <- .transform_features(d, Y_dd, ct = ct)
    kappa <- .estimate_kappa(Y_div$correlated, covariates = covar)
    Y_kappa_corrected <- Y
    for (cv in Appac@correction@covariates) {
      Y_kappa_corrected <- lapply(seq_along(samples), function(i)
        Y_kappa_corrected[[i]] / (1 + kappa[[1]][[cv]] * covar[[cv]][[i]])
      )
    }

    # browser()
    Y_div <- .transform_features(d, Y_kappa_corrected, ct = ct)
    # Restore the IDate class of the date column (col 1) for every sample —
    # was hardcoded to a single sample name ("X17k"), which is absent in any
    # other dataset and crashed .calculate() with NULL[,1] <-.
    for (smp in names(Y_div$correlated)) {
      Y_div$correlated[[smp]][, 1]   <- as.IDate(Y_div$correlated[[smp]][, 1])
      Y_div$uncorrelated[[smp]][, 1] <- as.IDate(Y_div$uncorrelated[[smp]][, 1])
    }
    daily_factor_trend <- .estimate_uncorrelated_model(Y_div$uncorrelated)
    daily_factor_features <- .estimate_correlated_model(Y_div$correlated)
    Appac@trend@center <- lapply(Y_div$correlated, function(x)
      attr(x, "mean")
    )
    Appac@correction@coefficients <- kappa[[1]]
    ## binned per-(sample x peak) input to the kappa fit, for
    ## plot_area_pressure_fit
    Appac@correction@fit.data <- kappa[[3]]

    ##*****************************************************
    ## evaluate the results
    ##*****************************************************
    Appac@trend@date <- as.integer(daily_factor_trend$date)
      # as.integer(daily_factor_trend$trend[, "date"])
    Appac@trend@bias <- daily_factor_trend$bias
      # [, "value"]
    Appac@trend@trend <- daily_factor_trend$trend
      # [, "spline"]
    Appac@trend@correlated.features <-
      daily_factor_features$daily.factor[, "value"]
    scaling <- lapply(daily_factor_trend$scaling, function(x)
      x
    )
    names(scaling) <- samples
    Appac@trend@bias.trend.scaling <- scaling

    scaling <- lapply(daily_factor_features$scaling, function(x)
      x
    )
    names(scaling) <- samples
    Appac@trend@correlated.features.scaling <- scaling
    trend <- lapply(seq_along(samples), function(i)
      .get_daily_value(
        date = d[[i]],
        data.frame(date = daily_factor_trend$date,
                   value = daily_factor_trend$trend))
    )

    trend <- lapply(seq_along(samples), function(i)
      trend[[i]] %*% t(scaling[[i]])
    )

    bias <- lapply(seq_along(samples), function(i)
      .get_daily_value(
        date = d[[i]],
        data.frame(date = daily_factor_trend$date,
                   value = daily_factor_trend$bias))
    )
    bias <- lapply(seq_along(samples), function(i)
      bias[[i]] %*% t(scaling[[i]])
    )


    correlated_features <- lapply(seq_along(samples), function(i)
      .get_daily_value(date = d[[i]],
                       daily_value = daily_factor_features$daily.factor)
    )
    correlated_features <- lapply(seq_along(samples), function(i)
      correlated_features[[i]] %*% t(scaling[[i]])
    )

    covariates_features <- list()
    for (cvt in Appac@correction@covariates) {
      covariates_features[[cvt]]  <- lapply(seq_along(samples), function(i)
        (unname
          (unlist(Appac@correction@samples[[i]][cvt]) -
           unlist(Appac@correction@reference.values[cvt])
          )
        ) * unlist(Appac@correction@coefficients[cvt])
      )
      names(covariates_features[[cvt]]) <- samples
    }
    summed_features <- lapply(seq_along(samples), function(i)
      rowSums(sapply(covariates_features, function(x) x[[i]]))
    )

    for (i in seq_along(samples)) {
      qu <- .get_daily_value(
        date = Appac@samples[[i]]$date,
        daily_value = cbind(
          date = Appac@trend@date,
          1 +
            (Appac@trend@bias + Appac@trend@trend) %*%
              t(Appac@trend@bias.trend.scaling[[i]]) +
            (Appac@trend@correlated.features %*%
              t(Appac@trend@correlated.features.scaling[[i]]))
        )
      )
      if (nrow(qu) != length(summed_features[[i]]))
        stop("Inconsistent sizes of 'qu' and 'summed_features'.")
      qu <- sweep(qu, 1, summed_features[[i]], "+")
      Appac@samples[[i]]$corrected.area <- Appac@samples[[i]]$raw.area / qu
      Appac@correction@samples[[i]]$corrected.area <- sweep(
        Appac@samples[[i]]$raw.area, 1, (1 + summed_features[[i]]), "/"
      )
    }

    return(Appac)
}

