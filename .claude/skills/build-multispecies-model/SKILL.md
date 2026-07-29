---
name: build-multispecies-model
description: >-
  Build a calibrated multi-species mizer model from a species-parameter data
  frame. Use whenever the user wants to create a MizerParams object with
  newMultispeciesParams(), set up an interaction matrix or fishing gears, bring
  the model to steady state with steady(), or calibrate/match it to observed
  biomasses, yields, or growth (calibrateBiomass, matchBiomasses, matchGrowth,
  calibrateYield, setBevertonHolt). Follow this ordered workflow rather than
  guessing at parameters or writing the dynamics by hand.
---

# Building a calibrated multi-species mizer model

Build a model in this order. Each `set...`/`match...`/`calibrate...`/`steady`
function **returns a new `MizerParams` object** — always reassign
(`params <- f(params, ...)`); never modify slots in place. Change species
parameters through `given_species_params(params) <- ...`, which triggers the
recalculation of dependent quantities.

## Step 1 — Assemble the species parameters

`species_params` is a data frame with **one row per species**. Only two things
are truly required:

- **`species`** — the species name.
- **`w_inf`** — the von Bertalanffy asymptotic weight (g). This is the required
  maximum-size parameter; mizer derives `w_max` (the computational grid boundary,
  default `1.5 * w_inf`), `w_repro_max`, and a default `w_mat` from it.

Everything else has a sensible default or is calculated. Commonly supplied:

| Column | Meaning |
|--------|---------|
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default ~100) |
| `sigma` | Width of the lognormal predation kernel (default ~1.3) |
| `k_vb` | von Bertalanffy K — used to derive `h` (and then `gamma`) if `h`/`gamma` absent |
| `h`, `gamma` | Max intake coefficient and search-volume coefficient (alternative to `k_vb`) |
| `alpha` | Assimilation efficiency (default 0.6) |
| `R_max` | Beverton–Holt maximum recruitment (density dependence) |
| `biomass_observed` | Observed biomass, for calibration (Step 4) |
| `yield_observed` | Observed yield, for calibration (Step 4) |

Units: weights in **grams**, lengths in **cm**, time in **years**. A CSV read
with `read.csv()` is a fine source; the package ships an example at
`system.file("extdata", "NS_species_params.csv", package = "mizer")`.

## Step 2 — Create the MizerParams object

```r
params <- newMultispeciesParams(species_params)
```

Optional arguments to `newMultispeciesParams()`:

- **`interaction`** — a species × species matrix of dimensionless overlaps in
  `[0, 1]` (1 = full interaction, the default for every pair). Scales encounter
  and predation mortality.
- **`gear_params`** — a data frame defining fishing gears (columns `gear`,
  `species`, `sel_func`, `catchability`, selectivity params like
  `knife_edge_size`). Omit it and mizer builds a default knife-edge gear that
  catches every species. Change gears later with `gear_params(params) <- ...` or
  `setFishing()`.
- **`no_w`, `min_w`, `max_w`** — size-grid resolution and range (default
  `no_w = 100`).

Inspect the result with `summary(params)`, `species_params(params)`,
`getInteraction(params)`, and `gear_params(params)`.

## Step 3 — Find the steady state

`newMultispeciesParams()` gives only a rough initial spectrum. Bring the model
to a steady state, which also sets the initial values used by later steps and by
`project()`:

```r
params <- steady(params)
```

`steady()` runs the dynamics to convergence **with births held fixed** (which
makes it settle reliably), then re-tunes the reproduction parameters to that
steady state — use its `preserve` argument to pick whether `reproduction_level`
(default), `R_max`, or `erepro` is held fixed. For a steady state that is
dynamically unstable, the experimental `steadyNewton()` solves the steady-state
equation directly.

## Step 4 — Calibrate to observations

Do this only if you have observed biomasses and/or yields (supplied as the
`biomass_observed` / `yield_observed` columns, optionally with
`biomass_cutoff` / `yield_cutoff` size thresholds). Typical order:

```r
params <- calibrateBiomass(params)   # scale kappa so total modelled ≈ total observed biomass
params <- matchBiomasses(params)     # adjust each species' abundance to its observed biomass
params <- matchGrowth(params)        # rescale h, gamma, ks, k so growth hits w_mat/w_inf on time
params <- steady(params)             # re-converge after the adjustments
```

For yield instead of biomass, use `calibrateYield()` (this needs
`yield_observed`). After any `match...`, re-run `steady()` and check
with `plotBiomassObservedVsModel(params)` / `plotYieldObservedVsModel(params)`.

## Step 5 — Tune density-dependent reproduction

Set how strongly reproduction is density-limited. `reproduction_level` is the
fraction of the maximum recruitment realised at steady state (0 = no density
dependence, closer to 1 = strong):

```r
params <- setBevertonHolt(params, reproduction_level = 0.25)
```

You can instead pass `R_max`, `erepro`, or a per-species named vector.

## Step 6 — Verify before projecting

```r
plotSpectra(params)                     # sensible, overlapping size spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)      # points near the 1:1 line?
```

When the model looks right, project it forward and analyse the results — see the
`run-simulation` and `analyse-and-plot` skills:

```r
sim <- project(params, t_max = 20, effort = 1)
```

## Saving and loading a model

A calibrated model is worth persisting so you don't rebuild it every session.
Use mizer's own save/restore functions rather than bare `saveRDS()` — they store
the model in a version-stable form:

```r
saveParams(params, "cod_model.rds")   # write a MizerParams to disk
params <- readParams("cod_model.rds")  # read it back
```

`saveSim()` / `readSim()` do the same for a `MizerSim` object.

## Common pitfalls

- Forgetting to reassign the return value (`steady(params)` without `params <-`)
  silently discards the work.
- Skipping `steady()` after a `match...`/`calibrate...` step leaves the model
  off its steady state.
- Editing `species_params(params) <-` instead of `given_species_params(params) <-`
  when you want dependent quantities (e.g. `h` from `k_vb`) recalculated.

---

## Project notes

If a `NOTES.md` file sits beside this one, read it too before you start: it
records what earlier work in *this* project found, and wins over the guidance
above wherever the two disagree.

Write what you learn about this project into that `NOTES.md`, creating it if
it is not there. Do not edit `SKILL.md`: it is installed by
`mizerAgents::setup_mizer_agent()` and your edits would be reported as a
conflict on the next run rather than kept.

A lesson that is true of mizer in general, rather than of this project,
belongs upstream where every project gets it. Tell the user, and offer to
report it at <https://github.com/sizespectrum/mizerAgents/issues>.
