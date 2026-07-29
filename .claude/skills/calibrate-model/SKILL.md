---
name: calibrate-model
description: >-
  Bring a mizer model to steady state and calibrate it to observed data. Use
  whenever the user wants to find the steady state (steady, projectToSteady,
  steadySingleSpecies), match modelled biomass, yield, or growth to observations
  (calibrateBiomass, matchBiomasses, calibrateYield, matchGrowth),
  set the level of density-dependent reproduction (setBevertonHolt), or diagnose
  why a model will not settle or reproduce observed values.
---

# Calibrating a mizer model to data

This skill covers the tune-to-data loop for an existing `MizerParams` object.
To build the model from scratch first, see the `build-multispecies-model` skill.

Every function here **returns a new `MizerParams`** — always reassign
(`params <- f(params, ...)`). Change species parameters through
`given_species_params(params) <- ...` so dependent quantities recalculate.

## The calibration loop

Observed data lives in the species-parameter columns `biomass_observed` and/or
`yield_observed` (optionally with `biomass_cutoff` / `yield_cutoff` size
thresholds below which observations are not counted). The usual order:

```r
params <- steady(params)             # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale resource level kappa so total
                                     #    modelled biomass matches total observed
params <- matchBiomasses(params)     # 3. adjust each species so its biomass
                                     #    matches its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- steady(params)             # 5. re-converge after the changes
```

Re-run `steady()` after **any** `match...`/`calibrate...` step — those functions
move parameters off the current steady state.

**Yield instead of biomass:** use `calibrateYield()` (scales overall abundance to
total observed yield), which reads `yield_observed`. Yield calibration depends on
the fishing setup, so make sure gears and effort are right first (see the
`set-up-fishing` skill).

**What each does, briefly:**

| Function | Adjusts | To match |
|---|---|---|
| `calibrateBiomass()` | `kappa` (resource level) | total community biomass |
| `matchBiomasses()` | per-species abundance | each `biomass_observed` |
| `calibrateYield()` | overall abundance scale | total community yield |
| `matchGrowth()` | `h`, `gamma`, `ks`, `k` | von Bertalanffy growth to `w_mat`/`w_inf` |

## Finding the steady state

- **`steady(params)`** — runs the dynamics to convergence **with births held
  fixed**, then stores the result as the initial state. Holding births constant
  lets the dynamics settle reliably onto *a* steady state, so this is the default
  first choice during setup and calibration.
- **`steadySingleSpecies(params)`** — sets each species' spectrum to its
  single-species steady form (also with births held fixed) given the current
  rates, without changing the resource; a fast way to get a sensible starting
  spectrum before `steady()`.
- **`projectToSteady(params)`** — the lower-level routine `steady()` builds on,
  but with **births responding dynamically** rather than held fixed; exposes
  `t_max`, `tol`, and `return_sim` if you need to watch convergence.
- **`steadyNewton(params)`** *(experimental)* — solves the steady-state equation
  directly, converging even when the steady state is dynamically unstable.

After converging, `steady()` and `steadySingleSpecies()` re-tune the reproduction
parameters so that density-dependent reproduction reproduces exactly that birth
rate at the new steady state. Use their `preserve` argument to choose whether
`reproduction_level` (default), `R_max`, or `erepro` is held fixed during that
re-tuning.

## Density-dependent reproduction

Set how strongly reproduction is density-limited. `reproduction_level` is the
fraction of maximum recruitment realised at steady state (0 = density
independent, closer to 1 = strongly limited):

```r
params <- setBevertonHolt(params, reproduction_level = 0.25)
```

Alternatively pass `R_max`, `erepro`, or a per-species named vector. This does
not change the steady state itself — it sets how the model responds to
perturbations away from it.

Read the current values back with `getReproductionLevel(params)`, which returns
the realised reproduction level per species — useful to check what a model was
tuned to before changing it.

## Diagnosing calibration problems

- **Model won't settle in `steady()`** — the initial spectrum is likely far off.
  Try `params <- steadySingleSpecies(params)` first, or reduce the step of a
  parameter change and re-run. Persistent instability is a case for
  `steadyNewton()`.
- **A species collapses to near-zero** — its mortality exceeds the growth it can
  fund; check `beta`/`sigma` (predation kernel), the interaction matrix row, and
  whether fishing mortality is too high.
- **Biomass matches but growth is wrong (or vice versa)** — alternate
  `matchGrowth()` and `matchBiomasses()`, re-running `steady()` between them;
  they pull on different parameters and usually converge in a few passes.
- **Check the fit** with `plotBiomassObservedVsModel(params)`,
  `plotYieldObservedVsModel(params)`, and `plotGrowthCurves(params)`.

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides
`tuneParams()`, a Shiny gadget that exposes sliders for the parameters above and
re-runs `steady()` live. It is **not** part of core mizer — install/load
`mizerExperimental` to use it.

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
