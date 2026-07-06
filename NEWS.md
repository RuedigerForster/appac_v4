# appac 4.0.3

* Added a package vignette (`vignette("appac")`) covering usage on `PLOT_FID` and
  the method: the multiplicative forward model, the PCA decomposition into
  correlated / uncorrelated / noise components, robust estimation of the common
  `kappa`, the NA-tolerant drift/daily-factor imputation, and the change-point
  detectors.

# appac 4.0.2

* Change-point detection moved from 'Rbeast' to 'strucchange': `get_changepoints()`
  now dates episode level breakpoints with a deterministic structural-break model
  (OLS-MOSUM test + BIC-optimal `breakpoints()`), dropping the heavy 'Rbeast'
  dependency and the need for a random seed.
* New `get_variance_changepoints()`: detects precision (variance) breakpoints on the
  noise-energy signal — the second-moment counterpart of `get_changepoints()`.
* New example dataset `Synth_data`: a compact, fully synthetic stress-test set with a
  known ground truth (attached as `attr(., "truth")`) — three samples, ten peaks,
  three episodes split by two planted level/variance breakpoints, with brown (AR(1))
  heavy-tailed noise at a 1% repeatability — for unit tests and examples.
* `appac()` now imputes missing area cells: peaks with up to 30% `NA` are filled by
  low-rank reconstruction (svdImpute / EM) before the fit; whole missing injections
  (staggered dates) are handled by the cross-sample reconstruction.
* `appac()` validates minimum-size and degenerate input (at least 3 samples, 2 peaks
  and 20 injections per sample, and non-constant areas), failing with an explanatory
  error instead of a deep numeric one.
* `show()` and `print()` methods for the `Appac`, `Compensation` and `Correction`
  classes: a compact summary at the console (`print()` also lists per-sample
  goodness-of-fit) instead of dumping the full object.
* `check_cols()` gains a `verbose` argument (default `FALSE`) that reports which
  column, peak and sample names were renamed.
* `debias_ct()` shows a progress bar during the chi-square minimisation sweep.
* Documented the package limitations (see `?appac-package`).

# appac 4.0.1

First CRAN release.

* `appac()` runs the correction pipeline: it decomposes per-cylinder peak areas
  by principal components into a pressure-correlated component and per-peak
  drift, estimates the common pressure-sensitivity coefficient `kappa` with a
  heavy-tail-robust fit on a drift-reduced signal, and removes slow drift plus a
  daily factor. Corrects the response of standard, atmosphere-open detectors
  (FID, and more weakly TCD).
* `check_cols()` validates and canonicalises the input columns (role-keyed, so
  the order of the mapping does not matter).
* `debias_ct()` refines the per-peak centres by closed-form chi-square
  minimisation, for an optional de-biased second pass.
* `goodness_of_fit()` reports, per peak, the reduced chi-square of the corrected
  areas against a noise-floor estimate.
* `get_changepoints()` provides Bayesian episode/breakpoint detection on the
  PC2 drift signal (via 'Rbeast').
* `plot_area_pressure()`, `plot_area_date()`, `plot_residuals()` and
  `plot_area_pressure_fit()` visualise a fitted object (require the suggested
  'ggplot2' / 'patchwork').
* Example dataset `PLOT_FID`: real FID injections from several control cylinders.
* Scope: APPAC is an *a posteriori* correction of already-measured areas. It has
  no forecasting ability and makes no prediction beyond the acquired data.
