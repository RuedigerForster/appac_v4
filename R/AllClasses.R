#' Drift / daily-factor model
#'
#' The temporal-drift component of a fitted \code{\link{Appac-class}} object:
#' the
#' daily-factor and bias/trend decomposition of the area drift.
#'
#' @slot date Integer date axis the model is evaluated on.
#' @slot bias Per-cylinder additive bias.
#' @slot trend Per-cylinder linear trend.
#' @slot correlated.features The common daily-factor signal.
#' @slot correlated.features.scaling Scaling applied to the correlated features.
#' @slot bias.trend.scaling Scaling applied to the bias/trend component.
#' @slot center Per-sample, per-peak centre (reference) values.
#' @slot samples Per-cylinder inputs to the drift model.
#' @name Compensation-class
#' @export
setClass("Compensation",
  slots = c(
    date = "numeric",
    bias = "numeric",
    trend = "numeric",
    correlated.features = "numeric",
    correlated.features.scaling = "list",
    bias.trend.scaling = "list",
    center = "list",
    samples = "list"
  ),
  prototype = c(
    date = NA_integer_,
    bias = NA_real_,
    trend = NA_real_,
    correlated.features = NA_real_,
    bias.trend.scaling = list(),
    correlated.features.scaling = list(),
    center = list(),
    samples = list()
  )
)

#' Pressure correction model
#'
#' The pressure-correction component of a fitted \code{\link{Appac-class}}
#' object.
#'
#' @slot covariates Names of the correction covariates (e.g. air pressure).
#' @slot coefficients Fitted coefficient (\code{kappa}) per covariate.
#' @slot reference.values Reference value (e.g. \code{P_ref}) per covariate.
#' @slot samples Per-cylinder correction inputs.
#' @slot fit.data Stored kappa-fit input, used by
#'   \code{\link{plot_area_pressure_fit}}.
#' @name Correction-class
#' @export
setClass("Correction",
  slots = c(
    covariates = "character",
    coefficients = "list",
    reference.values = "list",
    samples = "list",
    fit.data = "list"
  ),
  prototype = c(
    covariates = character(0),
    coefficients = list(),
    reference.values = list(),
    samples = list(),
    fit.data = list()
  )
)

#' Fitted APPAC model
#'
#' The object returned by \code{\link{appac}}: the corrected areas plus the
#' drift
#' and pressure-correction models.
#'
#' @slot samples Per-cylinder list carrying \code{raw.area},
#'   \code{corrected.area}, dates and pressure.
#' @slot trend The drift / daily-factor model (a
#' \code{\link{Compensation-class}}).
#' @slot correction The pressure correction (a \code{\link{Correction-class}}).
#' @name Appac-class
#' @export
setClass("Appac",
  slots = c(
    samples = "list",
    trend = "Compensation",
    correction = "Correction"
  ),
  prototype = c(
    samples = list()
  )
)
