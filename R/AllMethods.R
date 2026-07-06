## S4 methods for the appac classes (the AllMethods counterpart to
## AllClasses.R).
## Concise show() / print() so a fitted object prints a summary at the console
## instead of dumping every slot (the raw / corrected area matrices run to
## megabytes).

#' @rdname Appac-class
#' @param object the object to display.
#' @export
setMethod("show", "Appac", function(object) {
  s     <- object@samples
  ns    <- length(s)
  peaks <- unique(unlist(lapply(s, function(x) colnames(x$raw.area))))
  dts   <- unlist(lapply(s, function(x) as.integer(x$date)))
  cat("<Appac> atmospheric-pressure peak-area correction\n")
  cat(sprintf("  samples: %d%s\n", ns,
              if (ns)
                paste0(" (", paste(names(s), collapse = ", "), ")") else ""))
  cat(sprintf("  peaks  : %d%s\n", length(peaks),
              if (length(peaks))
                paste0(" (", paste(peaks, collapse = ", "), ")") else ""))
  if (any(is.finite(dts)))
    cat(sprintf("  data   : %d injections, %s to %s\n", length(dts),
                as.character(data.table::as.IDate(min(dts, na.rm = TRUE))),
                as.character(data.table::as.IDate(max(dts, na.rm = TRUE)))))
  co <- object@correction@coefficients
  for (cv in names(co))
    cat(sprintf("  kappa[%s]: %.4g  (P_ref = %g)\n",
                cv, co[[cv]], object@correction@reference.values[[cv]]))
  invisible(object)
})

#' @rdname Appac-class
#' @param x the \code{"Appac"} object to print.
#' @param ... ignored.
#' @export
setMethod("print", "Appac", function(x, ...) {
  show(x)                                            # the compact header
  gof <- tryCatch(goodness_of_fit(x), error = function(e) NULL)
  if (!is.null(gof)) {
    cat("  goodness of fit (reduced chi-square):\n")
    for (nm in names(gof)) {
      rc <- gof[[nm]]$reduced.chisq; rf <- rc[is.finite(rc)]
      if (length(rf))
        cat(sprintf("    %-10s median %.2f  [%.2f, %.2f]  over %d peaks\n",
                    nm, stats::median(rf), min(rf), max(rf), length(rc)))
      else
        cat(sprintf("    %-10s (no finite chi-square) over %d peaks\n",
                    nm, length(rc)))
    }
  }
  invisible(x)
})

#' @rdname Compensation-class
#' @param object the object to display.
#' @export
setMethod("show", "Compensation", function(object) {
  cat("<Compensation> drift / daily-factor model\n")
  cat(sprintf("  dates  : %d\n", length(object@date)))
  cat(sprintf("  samples: %d%s\n", length(object@center),
              if (length(object@center))
                paste0(" (", paste(names(object@center), collapse = ", "),
                       ")") else ""))
  invisible(object)
})

#' @rdname Correction-class
#' @param object the object to display.
#' @export
setMethod("show", "Correction", function(object) {
  cat("<Correction> pressure correction\n")
  cat(sprintf("  covariates: %s\n", paste(object@covariates, collapse = ", ")))
  co <- object@coefficients
  for (cv in names(co))
    cat(sprintf("  kappa[%s] : %.4g  (ref = %g)\n",
                cv, co[[cv]], object@reference.values[[cv]]))
  invisible(object)
})
