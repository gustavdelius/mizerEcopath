# mizerEcopath 0.3.2.9000

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
* Known limitation: on a bin-averaged model `rowSums(getDietMatrix())` does not
  equal `getConsumption()`, because `mizer::getDiet()` applies the prey-bin
  quadrature twice and comes out a factor `(1 + beta) / 2` too large, where
  `beta` is the size-grid ratio. Reported upstream as
  [sizespectrum/mizer#474](https://github.com/sizespectrum/mizer/issues/474).
