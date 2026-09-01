# Project notes

- Tuning-gadget undo logs and downloads contain `MizerParams` objects and use
  `saveParams()` and `readParams()` so extension metadata is preserved.
- Packaged lazy-load datasets are the deliberate exception to the persistence
  rule: `data-raw/upgrade_mizer_3_4.R` writes them with `usethis::use_data()`.
- mizerEcopath depends on mizerExperimental but does not register extension
  methods or marker classes, so no extension-dispatch migration was needed for
  the mizer 3.4 upgrade.
