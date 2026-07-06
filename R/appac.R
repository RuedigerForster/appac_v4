
#' Fit the APPAC pressure correction
#'
#' Fits the APPAC (Atmospheric Pressure Peak Area Correction) model and returns
#' the
#' corrected peak areas with diagnostics.  The data must already be
#' column-checked with \code{\link{check_cols}}.  \code{appac} drives the whole
#' pipeline: pivot the data, decompose it by principal components, estimate the
#' common pressure sensitivity \code{kappa} with a heavy-tail-robust fit on a
#' drift-reduced signal, estimate the drift / daily-factor model, and divide the
#' assembled correction out of the raw areas.
#'
#' Missing areas are tolerated: a peak with up to 30\% \code{NA} is kept and its
#' gaps are imputed by low-rank reconstruction before the fit; whole missing
#' injections (samples on staggered dates) are handled by the cross-sample
#' reconstruction in the drift model.
#'
#' @param data Long-format data frame with the canonical columns
#'   (Sample_Name, Peak_Name, Injection_Date, Air_Pressure, Raw_Area), e.g. the
#'   output of \code{\link{check_cols}}.
#' @param ct Optional fixed per-sample centres (a list, named by sample, of
#'   per-peak reference values).  \code{NULL} uses each peak's whole-series
#'   mean;
#'   pass the result of \code{\link{debias_ct}} for the de-biased centres.
#' @param P_ref Reference pressure (hPa) at which the correction is unity.
#' @return An object of class \code{"Appac"}: the \code{samples} slot holds, per
#'   cylinder, the \code{corrected.area}; \code{correction@@coefficients} the
#'   fitted \code{kappa}; and \code{trend} the drift / daily-factor model.
#' @seealso \code{\link{check_cols}}, \code{\link{debias_ct}},
#'   \code{\link{goodness_of_fit}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat   <- check_cols(PLOT_FID, acn)
#' ap    <- as.numeric(dat[, "Air_Pressure"])
#' P_ref <- mean(range(ap, na.rm = TRUE))   # mid-range reference pressure
#' \donttest{
#' fit <- appac(dat, P_ref = P_ref)
#' unlist(fit@correction@coefficients)   # the fitted common kappa
#' }
#' @export
appac <- function(
  data,
  ct = NULL,
  P_ref = 1013.25
)

{

  if (missing(data)) {
    stop("Input 'data' is missing")
  }

  if (!all(colnames(data) %in% column_names)) {
    stop(
      paste("Unknown column names.",
        "Please call 'check_cols(data, appac_colnames)'",
        "before calling 'appac()'.")
    )
  }

  #----------------------------------------------------------
  # check if P_ref is within the range of the input pressures
  #----------------------------------------------------------
  if (!is.numeric(P_ref) || P_ref < min(data$Air_Pressure, na.rm = TRUE) ||
        P_ref > max(data$Air_Pressure, na.rm = TRUE)) {
    stop("P_ref: ", P_ref, " is out of range: ",
      min(data$Air_Pressure, na.rm = TRUE), " <  P_ref < ",
      max(data$Air_Pressure, na.rm = TRUE)
    )
  }

  peaks <- unique(data$Peak_Name)
  spls  <- unique(data$Sample_Name)

  #----------------------------------------------------------
  # check if data contains all peaks for every file
  # nrow(data) must be a multiple of length(peaks)
  #----------------------------------------------------------
  if (nrow(data) %% length(peaks) != 0) {
    stop("Missing data points")
  }

  #----------------------------------------------------------
  # minimum-size / non-degeneracy guards, so undersized or degenerate input
  # fails with a clear message instead of a deep numeric error (loess span,
  # PCA subscript, zero-variance scaling, cross-sample drift PCA).
  #----------------------------------------------------------
  if (length(peaks) < 2L) {
    stop("appac needs at least 2 peaks to decompose the areas; got ",
         length(peaks), ".")
  }
  if (length(spls) < 3L) {
    stop("appac needs at least 3 samples/cylinders to separate the common ",
         "drift; got ", length(spls), ".")
  }
  ## each sample must carry the FULL peak set -- the pipeline needs a
  ## rectangular
  ## (sample x peak x date) grid; a sample missing peaks is reported here rather
  ## than crashing later in .prepare() with an opaque dimnames/date error (and
  ## it
  ## makes the injections-per-sample count below correct).
  peaks_per_sample <- tapply(data$Peak_Name, data$Sample_Name,
                             function(p) length(unique(p)))
  if (any(peaks_per_sample < length(peaks))) {
    short <- names(peaks_per_sample)[peaks_per_sample < length(peaks)]
    stop("These samples do not contain all ", length(peaks), " peaks: ",
         paste(sQuote(short), collapse = ", "), ".")
  }
  inj <- table(data$Sample_Name) / length(peaks)   # complete grid -> injections
  if (any(inj < 20)) {
    stop("Each sample needs at least 20 injections; too few in: ",
         paste(sQuote(names(inj)[inj < 20]), collapse = ", "), ".")
  }
  ## every KEPT (sample, peak) series must vary.  A peak that is >
  ## max_na_fraction NA in a sample is dropped downstream (see .prepare), so
  ## flag only series that would survive that drop yet have no variation
  ## (constant / all non-NA equal).
  key   <- paste(data$Sample_Name, data$Peak_Name, sep = "\r")
  grp   <- split(data$Raw_Area, key)   # group once, then two cheap passes
  na_fr <- vapply(grp, function(z) mean(is.na(z)), numeric(1))
  vv    <- vapply(grp, function(z) stats::var(z, na.rm = TRUE), numeric(1))
  bad_i <- na_fr <= max_na_fraction & (!is.finite(vv) | vv == 0)
  if (any(bad_i)) {
    stop("These sample/peak series have no variation (constant): ",
         paste(sub("\r", "/", names(grp)[bad_i]), collapse = ", "), ".")
  }

  Appac <- .prepare(data, P_ref)
  Appac <- .calculate(Appac, ct)

  return(Appac)
}
