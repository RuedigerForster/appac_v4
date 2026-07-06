# data-raw/Synth_data.R
# -----------------------------------------------------------------------------
# Reproducible generator for the bundled `Synth_data` stress-test fixture.
# Run with:  Rscript data-raw/Synth_data.R   (from the package root)
#
# A deliberately ADVERSARIAL synthetic set: 3 samples of largely different
# composition, 10 peaks each spanning a factor-of-5 dynamic range, over three
# 50-run episodes separated by two planted breakpoints:
#   * breakpoint 1  -- an abrupt +5% shift of the composition centre (persists),
#   * breakpoint 2  -- a small -2% centre shift plus a nominal +200% step in the
#                      measurement-noise variance (innovation sd x sqrt(3)).  The
#                      AR(1) memory and short episodes pull the realized spread to
#                      ~1.5x, so this break sits NEAR the detection floor rather
#                      than clearing it (a deliberately marginal variance break).
# On top of the breakpoints every classical adversary of the pressure fit is
# stacked at once: a small kappa buried under brown/reddened heavy-tailed noise,
# ~2% positive outlier contamination, and a slow per-sample drift that is
# confounded with the pressure signal.  It is tuned to remain RECOVERABLE --
# appac() still returns the planted kappa within a loose tolerance -- so it
# doubles as a demanding-but-passing testthat fixture.
#
# The planted ground truth is attached as attr(Synth_data, "truth").
# -----------------------------------------------------------------------------

## --- reproducibility --------------------------------------------------------
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(2026)

## --- planted parameters (the "truth") ---------------------------------------
kappa_true   <- -6.5e-4          # per hPa; small vs. the noise floor
P_ref        <- 1000             # hPa, correction is unity here (= pressure centre)
n_samples    <- 3L
n_peaks      <- 10L
n_episodes   <- 3L
runs_per_ep  <- 50L
n_runs       <- n_episodes * runs_per_ep          # 150 runs (dates) per sample
dyn_range    <- 5                 # within-sample max/min peak-area ratio
center_shift <- 0.05              # breakpoint 1: +5% level shift (episodes 2,3)
center_shift2 <- -0.02            # breakpoint 2: extra -2% level step (episode 3)
var_expand   <- 2.00              # breakpoint 2: +200% VARIANCE (episode 3); sd x sqrt(1+var)

repeatability <- 0.01            # base marginal noise sd = manufacturer repeatability spec (1%)
noise_df     <- 4                 # Student-t d.o.f. -> heavy tails
noise_phi    <- 0.85             # AR(1) coefficient -> brown/reddened noise (real data: acf1~0.87, alpha~1.4)
drift_amp    <- 0.010            # slow drift amplitude (confounds with kappa)
contam_frac  <- 0.02              # fraction of positive outlier spikes
contam_lo    <- 1.05
contam_hi    <- 1.15

## --- shared injection axis (same instrument days for every sample) ----------
## 150 distinct, sorted dates; pressure is a per-day atmospheric value shared
## across the cylinders (same day == same ambient pressure).
dates    <- as.integer(as.Date("2021-01-04")) + sort(sample(0:900L, n_runs))
pressure <- stats::rnorm(n_runs, P_ref, 30 / sqrt(3))  # Gaussian, sd = 30/sqrt(3) ~ 17.3 hPa
episode  <- rep(seq_len(n_episodes), each = runs_per_ep)      # 1..3 per run
bp_dates <- data.table::as.IDate(dates[c(runs_per_ep + 1L, 2L * runs_per_ep + 1L)])

## --- per-sample, per-peak reference composition -----------------------------
## Each sample gets a distinct overall scale AND a distinct within-sample shape
## (a fresh random pattern rescaled to span exactly a factor of `dyn_range`),
## so the three compositions are "largely different".
peaks    <- sprintf("pk%02d", seq_len(n_peaks))
ret_time <- round(cumsum(stats::runif(n_peaks, 0.8, 3.0)), 3)   # elution order
base_scale <- c(5e4, 2e4, 1e5)                                  # different magnitudes

ref <- matrix(NA_real_, nrow = n_samples, ncol = n_peaks,
              dimnames = list(sprintf("S%d", seq_len(n_samples)), peaks))
for (s in seq_len(n_samples)) {
  raw <- stats::runif(n_peaks, 0, 1)
  shape <- (raw - min(raw)) / (max(raw) - min(raw)) * (dyn_range - 1) + 1  # in [1, 5]
  ref[s, ] <- base_scale[s] * shape
}

## --- forward model ----------------------------------------------------------
## area = ref * centre_shift * (1 + kappa*(P - P_ref)) * (1 + drift) * (1 + noise)
## with heavy-tailed noise, episode-3 variance expansion, and outlier spikes.
rows <- vector("list", n_samples * n_peaks)
k <- 0L
for (s in seq_len(n_samples)) {
  phase <- 2 * pi * (s - 1) / n_samples                        # per-sample drift phase
  drift <- drift_amp * sin(2 * pi * (dates - min(dates)) / (diff(range(dates)) / 1.5) + phase)
  cfac   <- c(1, 1 + center_shift, (1 + center_shift) * (1 + center_shift2))
  cshift <- cfac[episode]                                      # bp1 +5%, then bp2 -2%
  sd_run <- ifelse(episode >= 3L, sqrt(1 + var_expand), 1)     # bp2: ep3 innovations x sqrt(1+var_expand)
  for (j in seq_len(n_peaks)) {
    mu    <- ref[s, j] * cshift * (1 + kappa_true * (pressure - P_ref)) * (1 + drift)
    innov <- stats::rt(n_runs, df = noise_df) * sd_run          # heavy-tailed, episode-scaled innovations
    noise <- as.numeric(stats::filter(innov, noise_phi, method = "recursive"))  # AR(1) -> brown/reddened
    noise <- noise - mean(noise)                                # DC offset absorbed by the peak centre
    noise <- noise / stats::sd(noise[episode <= 2L]) * repeatability  # pin base marginal sd = spec (1%)
    area  <- mu * (1 + noise)
    ci    <- sample.int(n_runs, round(contam_frac * n_runs))    # positive outlier spikes
    area[ci] <- area[ci] * stats::runif(length(ci), contam_lo, contam_hi)
    k <- k + 1L
    rows[[k]] <- data.frame(
      sample.name    = sprintf("S%d", s),
      injection.date = data.table::as.IDate(dates),
      peak.name      = peaks[j],
      retention.time = ret_time[j],
      raw.area       = area,
      air.pressure   = pressure,
      stringsAsFactors = FALSE)
  }
}
Synth_data <- do.call(rbind, rows)
## order like a real export: by date, then sample, then elution order
Synth_data <- Synth_data[order(Synth_data$injection.date,
                               Synth_data$sample.name,
                               Synth_data$retention.time), ]
rownames(Synth_data) <- NULL

## --- attach the ground truth ------------------------------------------------
reference <- data.frame(
  sample.name = rep(sprintf("S%d", seq_len(n_samples)), each = n_peaks),
  peak.name   = rep(peaks, times = n_samples),
  ref.area    = as.numeric(t(ref)),                # episode-1 centre (pre-shift)
  stringsAsFactors = FALSE)

attr(Synth_data, "truth") <- list(
  kappa              = kappa_true,
  P_ref              = P_ref,
  breakpoints        = bp_dates,
  center_shift       = center_shift,      # breakpoint 1: +5% level step
  center_shift_bp2   = center_shift2,     # breakpoint 2: extra -2% level step
  variance_expansion = var_expand,        # breakpoint 2: +200% variance (sd x sqrt(1+var_expand))
  n_episodes         = n_episodes,
  runs_per_episode   = runs_per_ep,
  dynamic_range      = dyn_range,
  reference          = reference,
  noise              = list(repeatability = repeatability, df = noise_df, phi = noise_phi,
                            drift_amp = drift_amp, contam_frac = contam_frac))

## --- save -------------------------------------------------------------------
save(Synth_data, file = "data/Synth_data.rda", compress = "xz")
cat(sprintf("Wrote data/Synth_data.rda: %d rows, %d samples x %d peaks x %d runs\n",
            nrow(Synth_data), n_samples, n_peaks, n_runs))
