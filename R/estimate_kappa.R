## Estimate the per-covariate sensitivity kappa from the inverse regression
## (covariate ~ peak areas) on the correlated (PC1) component.  The areas are
## binned by the covariate (0.1-unit bins) before fitting; for a linear model
## predict(area-ratio 1.1) - predict(0.9) = 0.2 * sum(beta), so kappa =
## 1/sum(beta).
##
## The fit is a heavy-tail-robust GLM (robustbase::glmrob, Mqle): the kappa-fit
## residuals are heavy-tailed, so a robust fit down-weights the tails instead of
## being dragged by them.
##
## Arguments:
##   correlated_data - per-sample correlated component (PC1) from
##   transform_features,
##                     each carrying attr "mean" (the day-1 reference) + a date
##                     column.
##   covariates      - named list of referenced covariate vectors per sample
##                     (e.g. pressure - reference value).
## Value:
##   list(kappa, rsq, fit_data).  kappa and rsq are named-by-covariate lists;
##   fit_data holds the binned per-(sample x peak) input for plotting / audit.
.estimate_kappa <- function(correlated_data, covariates) {

  summarize_values <- function(X, smp) {
    X %>%
      as_tibble() %>%
      group_by(.covar) %>%
      summarise(across(starts_with(smp), mean))
  }

  samples   <- names(correlated_data)
  cvt_names <- names(covariates)

  ## scale each peak by its reference (attr "mean" = the day-1 value)
  y_scld <- lapply(seq_along(samples), function(i)
    sweep(
      subset(correlated_data[[i]], select = -date), 2,
      attr(correlated_data[[i]], "mean"), "/"
    )
  )

  kappa <- rsq <- fit_data <- list()
  for (cvts in cvt_names) {
    for (i in seq_along(samples)) {
      idx <- order(covariates[[cvts]][[i]])
      colnames(y_scld[[i]]) <- sapply(colnames(y_scld[[i]]), function(x)
        paste0(samples[i], ".", x)
      )
      x <- round(covariates[[cvts]][[i]], 1)
      y <- cbind(x, y_scld[[i]])
      colnames(y)[1] <- ".covar"
      y <- y[idx, ]
      y <- summarize_values(y, samples[i])
      if (i == 1) {
        y_cvt <- y
      } else {
        y_cvt <- merge(y_cvt, y, by = ".covar", all = TRUE, no.dups = FALSE)
      }
    }

    ## keep the binned per-(sample x peak) input (pre-impute) for plotting /
    ## audit
    fit_data[[cvts]] <- y_cvt

    ## numeric design; mean-impute missing bins per column (some samples do not
    ## cover every covariate bin after the outer merge above)
    m  <- data.matrix(y_cvt)
    cm <- colMeans(m, na.rm = TRUE)
    for (j in seq_len(ncol(m))) m[is.na(m[, j]), j] <- cm[j]
    y <- m[, 1]                 # response: the covariate (e.g. pressure)
    x <- m[, -1, drop = FALSE]  # predictors: binned peak areas

    ## Inverse regression (covariate ~ peak areas) for kappa.  For a linear
    ## model
    ##   predict(area-ratio 1.1) - predict(0.9) = 0.2 * sum(beta),
    ## so the original  kappa = 0.2 / (x1 - x2)  is just  1 / sum(beta)
    ## (the intercept cancels) - computed directly from the coefficients here.
    ## Robust GLM (robustbase::glmrob, Mqle): the corrected-area distribution is
    ## heavy-tailed, so the robust fit down-weights the tails instead of being
    ## dragged by them.  complete.cases() drops incomplete rows; glmrob
    ## eliminates collinear (aliased) columns, returning NA for them (set to 0
    ## below).  Mqle can warn "did not converge" on a Gaussian response - the
    ## estimate is still usable.
    cc  <- stats::complete.cases(x, y)
    dat <- data.frame(.response = y[cc], x[cc, , drop = FALSE],
                      check.names = FALSE)
    best_model <- .quiet(
      robustbase::glmrob(.response ~ ., data = dat, family = gaussian()))
    cf <- coef(best_model)
    cf[is.na(cf)] <- 0
    betas    <- cf[-1]
    fitted_y <- as.numeric(fitted(best_model))
    y        <- y[cc]

    kappa[[cvts]] <- 1 / sum(betas)
    ## R-squared on the binned data
    rsq[[cvts]] <- as.numeric(1 - sum((y - fitted_y)^2) / sum((y - mean(y))^2))
  }

  list(kappa, rsq, fit_data)
}
