
#' Validate and canonicalise the input columns
#'
#' Validates that the input carries the required columns, renames them to the
#' canonical internal names, and cleans the sample and peak names to valid R
#' names with \code{\link[base]{make.names}}.  Call this before
#' \code{\link{appac}}.
#'
#' The mapping is keyed by \emph{role} (\code{sample_col}, \code{peak_col},
#' \code{date_col}, \code{pressure_col}, \code{area_col}), so the order in which
#' \code{appac_colnames} lists them does not matter: each role's column is found
#' by name and relabelled to its canonical name.
#'
#' @param data The raw long-format data frame.
#' @param appac_colnames Named list mapping the roles (\code{sample_col},
#'   \code{peak_col}, \code{date_col}, \code{pressure_col}, \code{area_col}) to
#'   the actual column names in \code{data}.
#' @param verbose If \code{TRUE}, report (via \code{\link[base]{message}}) which
#'   column, peak and sample names were renamed.  Off by default.
#' @return The data frame with columns renamed to the canonical names and
#'   sample/peak names cleaned.  Stops if a role or column is missing, or if
#'   \code{make.names()} would collapse two distinct sample or peak names into
#'   one.
#' @seealso \code{\link{appac}}
#' @examples
#' acn <- list(sample_col = "sample.name", peak_col = "peak.name",
#'             date_col = "injection.date", pressure_col = "air.pressure",
#'             area_col = "raw.area")
#' dat <- check_cols(PLOT_FID, acn)
#' head(dat)
#' @export
check_cols <- function(data, appac_colnames, verbose = FALSE) {

  ## canonical role order -- MUST line up with column_names
  roles <- c("sample_col", "peak_col", "date_col", "pressure_col", "area_col")

  #----------------------------------------------------------
  # validate, keyed by role so the order of appac_colnames is irrelevant
  #----------------------------------------------------------
  missing_role <- !roles %in% names(appac_colnames)
  if (any(missing_role)) {
    stop("appac_colnames is missing the role(s): ",
         paste(roles[missing_role], collapse = ", "))
  }
  # column names, in canonical role order
  wanted <- unname(unlist(appac_colnames[roles]))
  orig_col_names <- unname(vapply(colnames(data),
                                  function(x) gsub("\n", "", x), character(1)))
  missing_col <- !wanted %in% orig_col_names
  if (any(missing_col)) {
    stop("Could not find the column(s): '",
         paste0(wanted[missing_col], collapse = "', '"), "'")
  }

  #----------------------------------------------------------
  # extract + relabel: select each role's column, then assign the canonical
  # names
  # in the SAME (role) order -- independent of appac_colnames' order
  #----------------------------------------------------------
  idx       <- match(wanted, orig_col_names)
  orig_kept <- orig_col_names[idx]
  data      <- data[, idx]
  colnames(data) <- column_names

  #----------------------------------------------------------
  # make the sample and peak names compatible with R naming conventions, i.e.
  # convert special characters to '.' but keep the (upper, lower) case.  Guard
  # against make.names() collapsing two DISTINCT names into one (which would
  # silently merge peaks/samples).
  #----------------------------------------------------------
  clean_or_stop <- function(orig, what) {
    orig  <- unname(vapply(orig, function(x) gsub("\n", "", x), character(1)))
    clean <- make.names(orig)
    if (anyDuplicated(clean)) {
      dup <- clean %in% clean[duplicated(clean)]
      stop("Distinct ", what, " names collide under make.names(): '",
           paste0(orig[dup], collapse = "', '"),
           "' map to the same R name. Please make them distinct.")
    }
    list(orig = orig, clean = clean)
  }
  pk <- clean_or_stop(unique(data$Peak_Name),   "peak")
  sp <- clean_or_stop(unique(data$Sample_Name), "sample")

  for (i in seq_along(pk$orig))
    data$Peak_Name[data$Peak_Name == pk$orig[i]] <- pk$clean[i]
  for (i in seq_along(sp$orig))
    data$Sample_Name[data$Sample_Name == sp$orig[i]] <- sp$clean[i]

  #----------------------------------------------------------
  # report what was renamed (only when verbose)
  #----------------------------------------------------------
  if (verbose) {
    report <- function(orig, new, what) {
      if (!identical(orig, new)) {
        ind <- orig != new
        message(what, ": '", paste0(orig[ind], collapse = "', '"),
                "'\n\t have been replaced by: '",
                paste0(new[ind], collapse = "', '"), "'\n")
      }
    }
    report(orig_kept, column_names, "Column names")
    report(pk$orig, pk$clean, "Peak names")
    report(sp$orig, sp$clean, "Sample names")
  }

  return(invisible(data))
}
