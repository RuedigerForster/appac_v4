# Synthetic data with a KNOWN kappa, for CORRECTNESS tests (recover what we put
# in) and to demonstrate the headline objective: an >= 8x (1/8) reduction in the
# standard uncertainty u of a single area measurement.
#
# Forward model:  area[i,j,t] = ref[i,j] * (1 + kappa*(P[t]-P_ref)) * (1+drift)
# * (1+noise)
# P[t] is shared across cylinders (same day = same atmospheric pressure).
make_synthetic <- function(kappa = -5e-4, n_samples = 5, n_peaks = 4,
                           n_days = 800, noise = 0.002, drift = 0,
                           heavy = FALSE, contam = 0, P_ref = 1000, seed = 1) {
  set.seed(seed)
  dates    <- as.integer(as.Date("2020-01-01")) + sort(sample(0:1600, n_days))
  pressure <- stats::rnorm(n_days, P_ref, 12)
  peaks    <- paste0("pk", seq_len(n_peaks))
  ref_pk   <- 1e5 * (1 + seq_len(n_peaks))
  span     <- diff(range(dates))
  out <- list()
  for (s in seq_len(n_samples)) for (j in seq_len(n_peaks)) {
    ref  <- ref_pk[j] * stats::runif(1, 0.95, 1.05)
    dr   <- if (drift > 0)
      drift * sin(2 * pi * (dates - min(dates)) / (span / 1.5) + s) else 0
    base <- ref * (1 + kappa * (pressure - P_ref)) * (1 + dr)
    nz   <- if (heavy)
      stats::rt(n_days, df = 3) * noise else stats::rnorm(n_days, 0, noise)
    a    <- base * (1 + nz)
    if (contam > 0) {
      k <- sample(n_days, round(contam * n_days))
      a[k] <- a[k] * stats::runif(length(k), 1.05, 1.15)
    }
    out[[length(out) + 1L]] <- data.frame(
      sample.name = sprintf("S%d", s), peak.name = peaks[j],
      injection.date = dates, air.pressure = pressure, raw.area = a,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# generate -> check_cols -> appac
fit_synthetic <- function(...) {
  d  <- suppressMessages(check_cols(make_synthetic(...), acn))
  ap <- as.numeric(d[, "Air_Pressure"])
  pr <- mean(range(ap, na.rm = TRUE))
  suppressWarnings(appac(d, P_ref = pr))
}

kappa_of <- function(fit) unname(unlist(fit@correction@coefficients))

# per-(sample, peak) reduction in standard uncertainty u: raw RSD / corrected
# RSD
u_reduction <- function(fit) {
  rsd <- function(x) stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  unlist(lapply(seq_along(fit@samples), function(i)
    apply(fit@samples[[i]]$raw.area, 2, rsd) /
    apply(fit@samples[[i]]$corrected.area, 2, rsd)))
}
