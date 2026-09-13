# mizerEcopath 0.4.0

* **The Ecopath-to-mizer core has moved to the new
  [mizerFromEcopath](https://github.com/gustavdelius/mizerFromEcopath)
  package.** mizerEcopath now depends on it, so `library(mizerEcopath)` still
  puts `matchCatch()`, `getConsumption()`, `celtic_params` and the rest on the
  search path and no user code needs to change. What stays here is the
  interactive tuning gadget, the allometric starting model
  (`newAllometricParams()`, `setFeedingLevels()`), the growth-diffusion and
  age-at-size module, and the general mizer utilities that touch no Ecopath
  data (`bindParams()`, `addCatch()`, `plotPreyAvailability()`,
  `plotReproductiveEfficiency()`, `plotDeathX()`).

* `matchEcopath()`, `matchProduction()`, `matchProductionOnce()`,
  `matchExtMortOnce()`, `matchRespirationOnce()` and
  `matchGonadicProportionOnce()` have been removed rather than moved.
  `matchProductionOnce()` rescaled the metabolic rate without compensating the
  encounter rate, so it changed the energy available for growth -- on a model
  from `newVonBertalanffyParams()` it moved somatic growth by up to a factor of
  four. `matchExtMortOnce()` additionally discarded the result of its own
  `steady = TRUE` pipe, and `matchGonadicProportionOnce()` referenced an
  undefined variable.

* `plotYieldVsSizeByGear()` has been removed; it was never exported or called.

* mizerEcopath no longer contains compiled code. The TMB objective used by
  `matchCatch()` lives in mizerFromEcopath.

# mizerEcopath 0.3.2.9000

* mizerEcopath now requires mizer 3.4.0 and mizerExperimental 3.3.0. All
  S4-only model introspection and test expectations have been migrated to the
  S3 `MizerParams` representation introduced by mizer 3.4.
* The packaged `MizerParams` datasets have been upgraded to the S3
  representation. The tuning gadget now saves and restores models with
  `saveParams()` and `readParams()`, preserving mizer's validation, upgrade and
  extension metadata handling.
* Calls and documentation now use `tuneSteadyState()`,
  `projectUntilSettled()` and the current stored-rate accessors instead of
  their superseded aliases.
* Size-integrated Ecopath quantities now use mizer's `sizeIntegral()` helper,
  so future numerical-scheme changes are inherited from mizer rather than
  duplicated locally.
* The `get...()` functions are now consistent with `mizer::second_order_w()`.
  On a model with bin averaging switched on they trapezoidally bin-average the
  weight of their size integral, as mizer's own summary functions do, so that
  `getConsumption()`, `getSomaticProduction()`, `getZB()`, `getDietMatrix()` and
  the rest agree with the mizer quantities they correspond to. Models on mizer's
  default first-order path are unaffected.
* `getReproductiveEfficiency()` is now exactly the offspring biomass produced per
  unit of energy invested into reproduction whatever `mizer::second_order_w()` is
  set to. Previously that identity held only without bin averaging.
* `matchCatch()` now measures production, biomass and yield with the same
  quadrature as those functions. On a bin-averaged model this fixes a mismatch
  between the spectrum the optimisation scaled to `biomass_observed` and the one
  produced by the `mizer::matchBiomasses()` call that follows it.
* The catch panels of `tuneEcopath()` now report the same yield as
  `mizer::getYield()` on a bin-averaged model.
* mizer 3.3 fixed the double prey-bin quadrature in `mizer::getDiet()` that was
  tracked in [sizespectrum/mizer#474](https://github.com/sizespectrum/mizer/issues/474).
