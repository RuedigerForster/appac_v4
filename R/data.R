#' Flame-ionization-detector peak areas for APPAC
#'
#' An example set of gas-chromatography flame-ionization-detector (FID)
#' injections from a real instrument, with expert-annotated peak integration,
#' used to demonstrate the pressure correction.  Long format: one row per
#' injection and peak, from repeated analyses of several control cylinders.
#'
#' @format A data frame with 77365 rows and 6 columns:
#' \describe{
#'   \item{sample.name}{Cylinder (control sample) identifier.}
#'   \item{injection.date}{Date of the injection.}
#'   \item{peak.name}{Component / peak identifier.}
#'   \item{retention.time}{Peak retention time.}
#'   \item{raw.area}{Raw integrated peak area.}
#'   \item{air.pressure}{Ambient air pressure (hPa) at the injection.}
#' }
#' @keywords datasets
"PLOT_FID"

#' Synthetic stress-test data for APPAC
#'
#' A compact, fully synthetic data set with a \emph{known} ground truth, for
#' unit
#' tests and runnable examples.  It is deliberately adversarial: three samples
#' of
#' largely different composition, ten peaks each spanning a factor-of-five
#' dynamic
#' range, measured over three fifty-run episodes separated by two planted
#' breakpoints, with every classical adversary of the pressure fit stacked on
#' top
#' (a small \code{kappa} buried under brown (reddened, AR(1) \eqn{\phi = 0.85})
#' heavy-tailed Student-t noise, ~2\% positive outlier contamination, and a slow
#' per-sample drift confounded with the pressure signal).  The AR(1) coefficient
#' is
#' calibrated to the real \code{\link{PLOT_FID}} noise colour (lag-1
#' autocorrelation ~0.87, spectral slope ~1.4).
#'
#' The forward model is
#' \deqn{area = ref \cdot c_{ep} \cdot (1 + \kappa (P - P_{ref}))
#'   \cdot (1 + drift) \cdot (1 + noise),}
#' with a shared per-run air pressure \eqn{P} (same instrument day, same ambient
#' pressure across samples).  The \strong{composition is fixed} throughout --
#' the same samples, the same relative peak amounts -- and the two breakpoints
#' shift only the per-episode \emph{centre} \eqn{c_{ep}}, multiplicatively and
#' predictably:
#' \describe{
#'   \item{breakpoint 1}{an abrupt \eqn{+5\%} shift of the centre
#'     (episodes 2--3).}
#'   \item{breakpoint 2}{a further \eqn{-2\%} centre step \emph{and} a nominal
#'     \eqn{+200\%} step in the measurement-noise variance (innovation sd
#'     \eqn{\times\sqrt{3}}, episode 3).  The AR(1) memory and the short 50-run
#'     episodes pull the \emph{realized} spread ratio to about \eqn{1.5\times},
#'     so
#'     this break sits \emph{near} the detection floor rather than clearing it
#'     (see \sQuote{Limitations} in \code{\link{appac-package}}).}
#' }
#' Because the shifts are composition-preserving (common-mode) they are, by
#' design, hard for the trend/PC2 change-point detector to recover -- a naive
#' whole-series \code{\link{appac}} fit is likewise \emph{biased} by the
#' unmodelled
#' steps.  That is the point: the set exercises robustness and motivates episode
#' splitting rather than promising clean recovery.
#'
#' \strong{Two moments, two detectors.}  Breakpoint 2 deliberately carries two
#' \emph{separable} signatures that live in different statistical moments of the
#' per-sample PCA (standardise each peak, then decompose):
#' \itemize{
#'   \item the level steps (breakpoint 1, and breakpoint 2's \eqn{-2\%} step)
#'   are a
#'     \emph{first-moment} change -- a shift in the \emph{mean} of the dominant
#'     common-mode component (PC1).  A trend / level change-point model recovers
#'     these (from PC1, not the PC2 that \code{\link{get_changepoints}}
#'     currently
#'     reads; see there).
#'   \item the variance expansion is a \emph{second-moment} change -- a shift in
#'   the
#'     \emph{spread} of an uncorrelated (non-PC1) component, with its mean flat.
#'     A
#'     mean/trend model is blind to it by construction; it needs a
#'     \emph{dispersion}
#'     change-point statistic (e.g. a change-point on the rolling variance of
#'     the
#'     uncorrelated subspace).
#' }
#' Which component carries the variance step is sample-dependent (an artefact of
#' the
#' whole-series PCA ordering: PC2 in one sample, PC3 or higher in another), so a
#' variance detector must scan the uncorrelated subspace rather than a fixed PC.
#' Physically the two moments are distinct instrument events: a calibration /
#' level
#' shift (mean) versus a loss of precision / repeatability (variance).
#'
#' The planted ground truth is attached as \code{attr(Synth_data, "truth")}, a
#' list with \code{kappa}, \code{P_ref}, the \code{breakpoints}
#' (\code{\link[data.table]{IDate}}), the \code{center_shift} /
#' \code{center_shift_bp2}
#' / \code{variance_expansion} magnitudes, \code{n_episodes},
#' \code{runs_per_episode},
#' \code{dynamic_range}, the per-sample per-peak \code{reference} composition,
#' and
#' the \code{noise} parameters.  Regenerate with \code{data-raw/Synth_data.R}.
#'
#' @format A data frame with 4500 rows (3 samples x 10 peaks x 150 runs) and 6
#'   columns, matching \code{\link{PLOT_FID}}:
#' \describe{
#'   \item{sample.name}{Sample identifier (\code{S1}, \code{S2}, \code{S3}).}
#'   \item{injection.date}{Date of the injection
#'     (\code{\link[data.table]{IDate}}).}
#'   \item{peak.name}{Peak identifier (\code{pk01}..\code{pk10}).}
#'   \item{retention.time}{Peak retention time (elution order).}
#'   \item{raw.area}{Raw integrated peak area.}
#'   \item{air.pressure}{Ambient air pressure (hPa), shared per run across
#'     samples.}
#' }
#' @seealso \code{\link{PLOT_FID}}, \code{\link{appac}},
#' \code{\link{get_changepoints}}
#' @examples
#' truth <- attr(Synth_data, "truth")
#' truth$kappa            # the planted common pressure sensitivity
#' truth$breakpoints      # the two planted episode boundaries
#' @keywords datasets
"Synth_data"
