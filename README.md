# APPAC

**Atmospheric Pressure Peak Area Correction for Gas Chromatography with Standard Detectors**

APPAC corrects gas-chromatography peak areas for the influence of ambient air
pressure on **standard detectors** — those in open communication with the
ambient atmosphere, whose response can therefore vary with barometric pressure.
The flame ionization detector (FID) is the prime example, and the one
demonstrated here: standard textbooks and current reviews treat the FID response
as independent of ambient pressure, yet a real and correctable dependence
exists. The thermal conductivity detector (TCD) shows the same effect, more
weakly, and APPAC applies to it accordingly. Other open detectors (PID, FPD,
PFPD) have not yet been investigated. Detectors sealed from ambient pressure,
such as a mass spectrometer under vacuum, fall outside this class.

Per-cylinder peak areas are
decomposed by principal components into a pressure-correlated component and
per-peak drift; a single common pressure-sensitivity coefficient **kappa** is
estimated with a heavy-tail-robust fit on a drift-reduced signal, and the slow
drift plus a daily factor are removed. It returns the corrected areas with a
chi-square goodness-of-fit diagnostic, and provides Bayesian change-point
detection (via RBeast) for episode and breakpoint analysis.

## Scope — an *a posteriori* correction, not a prediction

APPAC is an **a posteriori** method: it corrects peak areas that have **already
been measured**, using the ambient pressure recorded at the time of each
injection, to remove a **known** pressure artifact. It has **no forecasting
ability** — it does not predict areas, pressures, or future measurements, and
infers nothing beyond the acquired data it is given.

## Installation

From a local clone of the repository:

```r
# install.packages("remotes")
remotes::install_local(".")
```

or, from the shell:

```sh
R CMD INSTALL .
```

The package depends only on CRAN packages; `ggplot2` and `patchwork` (for the
plots) are optional.

## Quick start

APPAC ships an example dataset, `PLOT_FID` — FID injections from several
control cylinders, recorded on a real instrument with expert-annotated peak
integration.

```r
library(appac)

# map your column names to the canonical roles (order does not matter)
acn <- list(
  sample_col   = "sample.name",
  peak_col     = "peak.name",
  date_col     = "injection.date",
  pressure_col = "air.pressure",
  area_col     = "raw.area"
)

data  <- check_cols(PLOT_FID, acn)
ap    <- as.numeric(data[, "Air_Pressure"])
P_ref <- (max(ap, na.rm = TRUE) - min(ap, na.rm = TRUE)) / 2 + min(ap, na.rm = TRUE)

# two passes: estimate kappa + drift, de-bias the centres, refit
fit1 <- appac(data, P_ref = P_ref)
ct   <- debias_ct(fit1, data = data, P_ref = P_ref, quiet = TRUE)
fit  <- appac(data, ct = ct, P_ref = P_ref)

# fitted common pressure sensitivity
unlist(fit@correction@coefficients)

# per-peak goodness of fit (reduced chi-square ~ 1 = down to the noise floor)
goodness_of_fit(fit)[[1]]
```

On `PLOT_FID` the correction reduces the run-to-run scatter (RSD) of the control
peaks from roughly 0.7 % to about 0.14 %.

## Plots

Four `ggplot2` plots visualise a fitted object:

```r
plot_area_pressure(fit, sample = 1, peak = "n.C4H10")     # area vs pressure, raw + corrected
plot_area_date(fit, sample = 1, peak = "n.C4H10")         # area vs date, with change-points
plot_residuals(fit, sample = 1, peak = "n.C4H10")         # residual diagnostic panel
plot_area_pressure_fit(fit)                               # what actually enters the kappa fit
```

A full worked walkthrough — usage and the method decomposition — is in the package
vignette: `vignette("appac")`.

## References

- Boček, P., Novák, J., Janák, J. (1969). *Effect of pressure on the performance
  of the flame ionization detector.* Journal of Chromatography.
  [doi:10.1016/S0021-9673(00)99223-9](https://doi.org/10.1016/S0021-9673(00)99223-9)
- Ayers, B. O., Clardy, E. K. (1985). *Pressure Compensation for a
  Chromatograph.* US Patent 4,512,181.
- Agilent Technologies (2005). *The importance of area and retention time
  precision in gas chromatography.* Technical Note 5989-3425EN.

## License

GPL-3. See [COPYING](COPYING).
