# Project notes: extending mizer in mizerEcopath

## How this package extends mizer

It does **not** customise mizer's dynamics. There is no `setRateFunction()`,
`setComponent()`, `other_params()`, `customFunction()`, S4 subclass or `.onLoad`
hook anywhere in `R/`. Everything is setup-time array surgery on stock slots
(`ext_encounter`, `mu_b`, `ext_diffusion`, `metab`, `interaction`), so
`project()` runs pure mizer. Keep it that way — it is the lightest mechanism the
skill asks for, and it means mizer's own solvers and diagnostics apply
unchanged.

The one `setRateFunction()` in the repo is a commented-out chunk in
`vignettes/Celtic_Sea_Status_Quo_Model.qmd` calling a function that no longer
exists in the tree.

## Quadrature: `flux` and `bin_average` are independent switches

`params@second_order_w` has two entries and they can differ —
`second_order_w(p) <- list(flux = "upwind", bin_average = TRUE)` is accepted.
`second_order_w(p) <- TRUE` flips both, which hides mistakes.

mizer gates on `bin_average` alone (not on the flux scheme) for: `@selectivity`
(`calc_selectivity()`), `@mu_b` (`setExtMort()`) and `@ext_diffusion`
(`setExtDiffusion()`). Anything here that rebuilds one of those arrays must use
the same switch. `src/objective_function.cpp` receives both as separate
`DATA_INTEGER`s — `second_order` for the steady-state flux scheme, `bin_average`
for quadrature. They were once conflated, which made the optimised selectivity
disagree with the model's own by a factor of 200 at the tail of the curve.

Test any new integral with `bin_average` **on** as well as off; the default path
alone proves nothing.

## `bin_average_weight()` is exported and takes two arguments

`mizer::bin_average_weight(K, params)` (mizer >= 3.2.1.9000) does the gating
itself — use it directly, do not wrap it and do not reach for
`mizer:::bin_average_weight(K)`, which is the same object and errors on the
missing `params`. This package briefly carried a private `binAverageWeight()`
duplicate; it is gone.

`mizer:::power_law_bin_average()`, `mizer:::flux_limiter_scheme()` and
`mizer:::flux_limiter_chi()` really are unexported. Prefer calling the public
setter that wraps them (`setExtDiffusion()` reproduces
`D_ext * power_law_bin_average(...)` exactly, in both schemes) over
reimplementing the formula.

## Recording a species parameter as "given"

A MizerParams object holds each species parameter twice: in
`given_species_params` (what was supplied) and `species_params` (what is used).
On any recalculation the second is rebuilt from the first, so a value that only
ever reached `species_params` is reset to its default. Neither public route does
what this package usually needs:

- `given_species_params(p)$x <- v` always calls `setParams()`, **and** rebuilds
  `species_params` from the given table alone — on a model carrying Ecopath
  columns (`consumption_observed`, `production_observed`, `gonad_proportion`, …)
  that drops all of them. Safe only on a params object fresh from
  `newMultispeciesParams()`, as in `newAllometricParams()`.
- `species_params(p, recalculate = FALSE) <- sp_new` does not recalculate, but
  only records parameters whose value it *sees change*. It cannot pin a value
  that is already correct but unrecorded — which is exactly the common case.

When you need to record without recalculating, write
`p@given_species_params$x <- v` and `p@species_params$x <- v` together, with a
comment saying why. `makeNoninteracting()` in `R/helpers.R` is the worked
example.

## Models must be recalculation-safe

Check any MizerParams this package produces or ships with:

```r
q <- p; species_params(q) <- species_params(q)
# every column of species_params should be unchanged
```

`celtic_params` failed this for six parameters, `interaction_resource` and
`pred_kernel_type` among them: one species-parameter change turned resource
feeding back on and swapped the `power_law` kernel for the lognormal default,
moving consumption by up to 66%. Repaired by
`data-raw/repair_celtic_params.R`, which records the calculated values and
verifies the model is otherwise byte-identical.

**`ns_3_spp_model_initial`, `ns_3_spp_model_trial` and `ns_3_spp_model_final`
still have this defect** (and also warn that they predate the installed mizer).
No test uses them, so they were left alone.

## `h = Inf` makes the search volume infinite

`newAllometricParams()` and `newVonBertalanffyParams()` switch off satiation with
`sp$h <- Inf`. mizer derives `gamma` from `h` and `f0`, so that gives an infinite
`gamma` and an infinite `search_vol`, which `check_finite()` now rejects — the
whole `newAllometricParams()` test file was erroring on it. Pin the calculated
`gamma` as a given parameter *before* setting `h`; the search volume has nothing
to do with satiation, and the model still needs a usable one for
`makeInteracting()` later.

## Known remaining issues

- `R/update_params.R:102-107` duplicates `setExtDiffusion()` via
  `mizer:::power_law_bin_average` instead of calling it.
- `R/update_params.R:95` writes `z_ext * w^d` straight into `ext_mort`,
  skipping the bin averaging `setExtMort()` would apply (1.4%).
- The `z_ext` / `D_ext` read-backs (`matchCatch.R:157`, `update_params.R:26`,
  `diffusion_params.R:26`, `ecopathOtherControl.R:26`) invert a possibly
  bin-averaged array with a point value.
- `makeInteracting()` records `interaction_resource` with the diff-based
  `species_params(recalculate = FALSE) <-`, which as noted above cannot pin an
  already-correct value.
- `plotPreyAvailability.R:64` decomposes the encounter integral using
  `getPredKernel()`; `mizer::encounter_kernel()` is the bin-integrated object
  the convolution actually consumes.
