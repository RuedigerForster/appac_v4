# Global RNG state for the whole test run, so any stochastic step that does not
# set its own seed is reproducible (belt-and-suspenders alongside
# make_synthetic()'s
# internal seed; the change-point detectors are deterministic via strucchange).
# Pinning RNGkind too keeps results stable across R versions (sampler defaults
# changed in R 3.6).  Package FUNCTIONS never call set.seed() themselves -- that
# would clobber the caller's RNG -- so this lives in the test setup only.
suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
set.seed(2026)
