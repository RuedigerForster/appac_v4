## ggplot2 plots for a fitted "Appac" object.  Each takes the object plus a
## sample and peak selector (name or positive index).  Every visual setting is
## configurable in global variables: COLOURS in default_palette, GEOMETRY
## (point/line sizes, alpha, bins, ...) in default_aes, and TEXT/LAYOUT in the
## .*_plot_theme() functions -- all in AllGlobalVariables.R.  The plot functions
## hold no hard-coded aesthetics.  ggplot2 (and, for the residual panel,
## patchwork) are Suggested dependencies; the functions stop clearly if absent.

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for appac plots; ",
         "install it with install.packages(\"ggplot2\").", call. = FALSE)
}

## resolve a sample/peak selector (name or positive index) to a name
.resolve <- function(x, choices, what) {
  if (is.numeric(x)) {
    if (x < 1 || x > length(choices))
      stop(sprintf("%s index %d is out of range (1..%d)",
                   what, x, length(choices)))
    return(choices[x])
  }
  if (!x %in% choices)
    stop(sprintf("unknown %s '%s'; available: %s",
                 what, x, paste(choices, collapse = ", ")))
  x
}

## long-format raw + corrected area for one (sample, peak), carrying date,
## pressure and the reference (centre) as attributes
.plot_data <- function(appac, sample, peak) {
  s    <- .resolve(sample, names(appac@samples), "sample")
  raw  <- appac@samples[[s]]$raw.area
  corr <- appac@samples[[s]]$corrected.area
  pk   <- .resolve(peak, colnames(raw), "peak")
  meta <- data.frame(date     = as.Date(appac@samples[[s]]$date),
                     pressure = appac@correction@samples[[s]]$pressure)
  long <- rbind(
    data.frame(meta, area = as.numeric(raw[,  pk]), series = "raw"),
    data.frame(meta, area = as.numeric(corr[, pk]), series = "corrected")
  )
  long$series <- factor(long$series, levels = c("raw", "corrected"))
  attr(long, "sample")    <- s
  attr(long, "peak")      <- pk
  attr(long, "reference") <- appac@trend@center[[s]][[pk]]
  long
}

.area_colours <- function() c(raw       = default_palette$lowlight_color,
                              corrected = default_palette$highlight_color)

## shared raw/corrected decoration: colour scale (with a readable legend key) +
## main theme + a top legend with readable labels
.area_decor <- function(size) {
  list(
    ggplot2::scale_colour_manual(
      values = .area_colours(), name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(size = default_aes$legend_key_size, alpha = 1))),
    .main_plot_theme(size),
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(
                     size = ggplot2::rel(default_aes$legend_text_rel)))
  )
}

## reference (centre) horizontal line, configurable
.ref_line <- function(yintercept)
  ggplot2::geom_hline(yintercept = yintercept,
                      colour = default_palette$line_color,
                      linewidth = default_aes$reference_linewidth)

.area_points <- function()
  ggplot2::geom_point(size = default_aes$point_size,
                      alpha = default_aes$point_alpha)

#' Plot peak area versus air pressure
#'
#' Scatter of raw and corrected peak area against air pressure for one
#' (sample, peak), with the reference (centre) line.  The pressure dependence
#' visible in the raw area should be flattened in the corrected area.
#'
#' @param appac A fitted \code{"Appac"} object.
#' @param sample Sample selector: a cylinder name or a positive index.
#' @param peak Peak selector: a peak name or a positive index.
#' @param size Base font size passed to the theme.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot_area_date}}, \code{\link{plot_residuals}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   plot_area_pressure(fit, sample = 1, peak = 1)
#' }
#' @export
plot_area_pressure <- function(appac, sample = 1, peak = 1, size = 12) {
  .need_ggplot2()
  d <- .plot_data(appac, sample, peak)
  ggplot2::ggplot(d, ggplot2::aes(x = pressure, y = area, colour = series)) +
    .area_points() +
    .ref_line(attr(d, "reference")) +
    ggplot2::labs(title = "Peak area vs. air pressure",
                  subtitle = sprintf("%s / %s",
                                     attr(d, "sample"), attr(d, "peak")),
                  x = "air pressure", y = "peak area") +
    .area_decor(size)
}

#' Plot peak area versus date
#'
#' Scatter of raw and corrected peak area against injection date for one
#' (sample, peak), with the reference (centre) line and, optionally, the robust
#' episode change-points as dotted vertical lines.
#'
#' @param appac A fitted \code{"Appac"} object.
#' @param sample Sample selector: a cylinder name or a positive index.
#' @param peak Peak selector: a peak name or a positive index.
#' @param show_changepoints Draw the detected change-points (see
#'   \code{\link{get_changepoints}}) as dotted vertical lines.
#' @param size Base font size passed to the theme.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot_area_pressure}}, \code{\link{get_changepoints}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   plot_area_date(fit, sample = 1, peak = 1, show_changepoints = FALSE)
#' }
#' @export
plot_area_date <- function(appac, sample = 1, peak = 1,
                           show_changepoints = TRUE, size = 12) {
  .need_ggplot2()
  d <- .plot_data(appac, sample, peak)
  s <- attr(d, "sample")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = date, y = area, colour = series)) +
    .area_points() +
    .ref_line(attr(d, "reference"))
  if (isTRUE(show_changepoints)) {
    bp <- tryCatch(as.Date(get_changepoints(appac@samples[s])),
                   error = function(e) as.Date(integer(0)))
    if (length(bp))
      p <- p + ggplot2::geom_vline(
        xintercept = bp,
        linetype   = default_aes$changepoint_linetype,
        colour     = default_palette$line_color,
        linewidth  = default_aes$changepoint_linewidth)
  }
  p +
    ggplot2::labs(title = "Peak area vs. date",
                  subtitle = sprintf("%s / %s", s, attr(d, "peak")),
                  x = "date", y = "peak area") +
    .area_decor(size)
}

#' Residual diagnostic panel
#'
#' A 2x2 panel of residual diagnostics for one (sample, peak): histogram with a
#' fitted normal, normal Q-Q plot, residual versus date, and residual versus
#' pressure.  The residual is the corrected area minus the reference (centre).
#'
#' @param appac A fitted \code{"Appac"} object.
#' @param sample Sample selector: a cylinder name or a positive index.
#' @param peak Peak selector: a peak name or a positive index.
#' @param size Base font size passed to the theme.
#' @return A \code{patchwork} object assembling the four panels (requires the
#'   \pkg{patchwork} package).
#' @seealso \code{\link{plot_area_date}}, \code{\link{goodness_of_fit}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE))
#'   plot_residuals(fit, sample = 1, peak = 1)
#' }
#' @export
plot_residuals <- function(appac, sample = 1, peak = 1, size = 12) {
  .need_ggplot2()
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for the residual panel; ",
         "install it with install.packages(\"patchwork\").", call. = FALSE)
  d   <- .plot_data(appac, sample, peak)
  s   <- attr(d, "sample"); pk <- attr(d, "peak"); ref <- attr(d, "reference")
  cor <- d[d$series == "corrected", ]
  r   <- data.frame(date = cor$date, pressure = cor$pressure,
                    residual = cor$area - ref)
  r   <- r[is.finite(r$residual), ]
  hi  <- default_palette$highlight_color; lo <- default_palette$lowlight_color
  ln  <- default_palette$line_color;      fl <- default_palette$fill_color
  pts <- function(x, y)
    ggplot2::geom_point(ggplot2::aes(x = {{x}}, y = {{y}}), colour = hi,
                        size = default_aes$point_size,
                        alpha = default_aes$point_alpha)

  g_hist <- ggplot2::ggplot(r, ggplot2::aes(x = residual)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = default_aes$histogram_bins, fill = fl,
                            colour = NA) +
    ggplot2::stat_function(fun = stats::dnorm,
                           args = list(mean = mean(r$residual),
                                       sd = stats::sd(r$residual)),
                           colour = ln, linewidth = default_aes$fit_linewidth) +
    ggplot2::labs(title = "Histogram", x = "residual", y = NULL) +
    .histogram_plot_theme(size)

  g_qq <- ggplot2::ggplot(r, ggplot2::aes(sample = residual)) +
    ggplot2::stat_qq(colour = lo, size = default_aes$point_size,
                     alpha = default_aes$point_alpha) +
    ggplot2::stat_qq_line(colour = ln, linewidth = default_aes$fit_linewidth) +
    ggplot2::labs(title = "Normal Q-Q", x = "theoretical", y = "sample") +
    .residuals_plot_theme(size)

  g_date <- ggplot2::ggplot(r) + pts(date, residual) +
    ggplot2::geom_hline(yintercept = 0, colour = ln,
                        linewidth = default_aes$reference_linewidth) +
    ggplot2::labs(title = "Residual vs. date", x = "date", y = "residual") +
    .residuals_plot_theme(size)

  g_press <- ggplot2::ggplot(r) + pts(pressure, residual) +
    ggplot2::geom_hline(yintercept = 0, colour = ln,
                        linewidth = default_aes$reference_linewidth) +
    ggplot2::labs(title = "Residual vs. pressure", x = "air pressure",
                  y = "residual") +
    .residuals_plot_theme(size)

  patchwork::wrap_plots(g_hist, g_qq, g_date, g_press, ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("Residual diagnostics - %s / %s", s, pk),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
}

#' Plot the kappa-fit input
#'
#' Shows the data that actually feeds the \code{kappa} fit: the reference-scaled
#' correlated area, binned by pressure deviation, for each (sample, peak) and
#' coloured by peak, with the fitted \code{kappa}-slope line overlaid.  All
#' series share the single common slope.
#'
#' Pass the \strong{de-biased} result -- \code{appac(ct = debias_ct(...))} -- so
#' this shows the chi-square-minimum fit (the correct \code{kappa}), not the
#' first pass.
#'
#' @param appac A fitted \code{"Appac"} object whose \code{correction} carries
#' the
#'   stored kappa-fit input (populated by \code{\link{appac}}).
#' @param covariate Covariate selector: a covariate name or a positive index
#'   (normally the pressure covariate, index 1).
#' @param size Base font size passed to the theme.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{appac}}, \code{\link{debias_ct}}
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
#' fit  <- appac(dat, ct = ct, P_ref = P_ref)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   plot_area_pressure_fit(fit)
#' }
#' @export
plot_area_pressure_fit <- function(appac, covariate = 1, size = 12) {
  .need_ggplot2()
  fd <- appac@correction@fit.data
  if (!length(fd))
    stop("no stored kappa-fit data on this object; was it produced by appac()?",
         call. = FALSE)
  cv <- if (is.numeric(covariate)) names(fd)[covariate] else covariate
  d  <- fd[[cv]]
  dev    <- d[[".covar"]]
  series <- setdiff(colnames(d), ".covar")
  long <- do.call(rbind, lapply(series, function(col)
    data.frame(deviation = dev, area = d[[col]], series = col)))
  long <- long[is.finite(long$area), ]
  ## split "sample.peak" back into sample + peak via the known sample names
  samples <- names(appac@samples)
  long$sample <- vapply(long$series, function(z) {
    hit <- samples[startsWith(z, paste0(samples, "."))]
    if (length(hit)) hit[which.max(nchar(hit))] else NA_character_
  }, character(1))
  long$peak <- substring(long$series, nchar(long$sample) + 2L)

  kappa <- appac@correction@coefficients[[cv]]
  ## kappa-slope line through the data centroid (scaled area = 1 + kappa * dev)
  b0 <- mean(long$area, na.rm = TRUE) -
    kappa * mean(long$deviation, na.rm = TRUE)
  ggplot2::ggplot(long, ggplot2::aes(x = deviation, y = area, colour = peak)) +
    ggplot2::geom_point(size = default_aes$point_size,
                        alpha = default_aes$point_alpha) +
    ggplot2::geom_abline(intercept = b0, slope = kappa,
                         colour = default_palette$line_color,
                         linewidth = default_aes$fit_linewidth) +
    ggplot2::scale_colour_manual(values = okabe_ito_palette, name = "peak") +
    ggplot2::labs(title = paste0("Kappa-fit input: reference-scaled area ",
                                 "vs. pressure deviation"),
                  subtitle = sprintf("%s   kappa = %.3g", cv, kappa),
                  x = sprintf("%s deviation from reference", cv),
                  y = "reference-scaled area") +
    .main_plot_theme(size) +
    ggplot2::theme(legend.position = "right")
}
