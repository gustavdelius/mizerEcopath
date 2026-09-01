# Project notes

- The package now targets mizer 3.4.0 and mizerExperimental 3.3.0. Because
  mizer 3.4.0 is newer than the current CRAN release, `DESCRIPTION` includes
  the sizespectrum/mizer remote.
- The four packaged `MizerParams` datasets were regenerated as S3 objects by
  `data-raw/upgrade_mizer_3_4.R`. Package data remains `.rda`; runtime and
  user-facing persistence uses `saveParams()` and `readParams()`.
- `bindParams()` must use list-style `[[` access for dynamically named model
  elements and validate its result with `validParams()`.
- With second-order size integration, `rowSums(getDietMatrix())` and
  `getConsumption()` can differ by about 2.7% at an occupied-grid boundary.
  The old double prey-quadrature bug is fixed in mizer 3.3; this remaining
  difference comes from `getDiet()` masking empty predator bins before its
  predator-size quadrature.
