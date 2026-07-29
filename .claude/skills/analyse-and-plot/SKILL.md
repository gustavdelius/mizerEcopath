---
name: analyse-and-plot
description: >-
  Analyse and visualise the results of a mizer simulation or the state of a
  MizerParams object. Use whenever the user wants to extract, summarise, or plot
  size spectra, biomass, yield, SSB, abundance, feeding level, mortality, diet,
  trophic level, community indicators, growth curves, or the plankton resource —
  including comparing two simulations or animating spectra through time. Prefer
  the built-in mizer functions described here over writing custom extraction or
  ggplot code.
---

# Analysing and plotting mizer results

mizer ships a large family of extraction, summary, and plotting functions.
**Always prefer these over hand-written array wrangling or custom ggplot code** —
they handle size-range integration, species colours/linetypes, and units for you.

Most functions accept **either** a `MizerSim` object (returning a **time series**)
**or** a `MizerParams` object (returning a **single value** from the initial
state). So `getBiomass(sim)` gives biomass over time, `getBiomass(params)` gives
biomass now.

To get the single value **at one time step** of a simulation, extract a
`MizerParams` snapshot with `finalParams(sim)` (last step), `initialParams(sim)`
(first step), or `getParams(sim, time_range = ...)` (averaged over a range) and
pass that in. Prefer this over indexing the time series with `idxFinalT(sim)`:

```r
getMeanMaxWeight(finalParams(sim))          # equilibrium value at the last step
getMeanMaxWeight(sim)[idxFinalT(sim), ]     # older equivalent
```

If you need a function you don't see here, `grep` for `"plot"` or the specific
name in the bundled API index (path at the end of `MIZER-AGENTS.md`) before
writing custom code — don't read the whole file. The index gives you the name;
read the help page for the arguments.

## 1. Accessing raw arrays

| Function | Returns | Dimensions |
|---|---|---|
| `N(sim)` | species abundance density | time × species × size |
| `NResource(sim)` | resource abundance density | time × size |
| `finalN(sim)` | abundance at last time | species × size |
| `finalNResource(sim)` | resource abundance at last time | size |
| `getEffort(sim)` | fishing effort | time × gear |
| `getTimes(sim)` | saved time steps | vector |

```r
N(sim)[, "Cod", ]        # time × size for Cod
N(sim)["2010", "Cod", ]  # size vector for Cod in 2010
finalN(sim)["Cod", ]     # size vector for Cod at final step
```

## 2. Summary functions

All accept `MizerSim` or `MizerParams`. See `?summary_functions`.

| Function | Returns |
|---|---|
| `getBiomass(sim, min_w, max_w)` | total biomass (time × species) |
| `getSSB(sim)` | spawning stock biomass |
| `getN(sim, min_w, max_w)` | total abundance |
| `getYield(sim)` | yield summed over gears |
| `getYieldGear(sim)` | yield by gear |
| `getFeedingLevel(sim)` | feeding level at size |
| `getPredMort(sim)` | predation mortality at size |
| `getFMort(sim)` | fishing mortality at size |
| `getDiet(params)` | diet resolved by prey (predator × size × prey) |
| `getTrophicLevel(params)` | trophic level at size |

`getBiomass()` and `getN()` take `min_w`/`max_w` (or `min_l`/`max_l`) to restrict
to a size range, e.g. `getBiomass(sim, min_w = 10, max_w = 1e4)`.

## 3. Community indicator functions

All accept `MizerSim` (time series) or `MizerParams` (single value). See
`?indicator_functions`.

`getProportionOfLargeFish(sim, threshold_w = 100)`, `getMeanWeight(sim)`,
`getMeanMaxWeight(sim, measure = "both")`, `getCommunitySlope(sim, min_w, max_w)`.

## 4. Named plotting functions

Every named plot returns a **ggplot2 object** you can extend with `+`, and has a
`plotly...()` counterpart (e.g. `plotlyBiomass()`) for interactive use. Restrict
species with `species = c(...)`. See `?plotting_functions`.

**Against time:** `plotBiomass(sim)`, `plotYield(sim)`, `plotYieldGear(sim)`.
**Against size:** `plotSpectra(sim)`, `plotFeedingLevel(sim)`, `plotPredMort(sim)`,
`plotFMort(sim)`, `plotGrowthCurves(sim)`, `plotDiet(params, species = "Cod")`.
**Calibration:** `plotBiomassObservedVsModel(params)`, `plotYieldObservedVsModel(params)`.
**Overview:** `plot(sim)` combines several panels; `plot(params)` shows the same
panels for a model's steady state.

**Common arguments** — most plot and summary functions (named or array) share:
`species` (subset), `tlim = c(min, max)` (restrict the time axis of a
time-series plot), `time_range` (average over a period, for plots against size),
`wlim`/`ylim` (restrict the visible size/value window), `highlight`, `total`, and
`log_x`/`log_y`. Note `tlim` replaces the deprecated `start_time`/`end_time`, and
`log_x`/`log_y` replace the older single `log`. `wlim`/`ylim` only set the
**visible window** — to change the numbers (e.g. the size range a biomass is
summed over) pass `min_w`/`max_w` (or `min_l`/`max_l`) to the `get...()` function
instead, e.g. `plotBiomass(sim, min_w = 10)`.

## 5. Plotting arrays directly (the compositional API)

The `get...()` functions return classed arrays that have their own `plot()`
methods, so you can visualise **any** quantity without a named plot function or
custom ggplot code:

- `ArrayTimeBySpecies` (e.g. `getSSB(sim)`, `getBiomass(sim)`) → time series per
  species.
- `ArraySpeciesBySize` (e.g. `getFeedingLevel(params)`, `getPredMort(params)`,
  `getFMort(params)`, `getEncounter(params)`) → curve of value against body size.
- `ArrayTimeBySpeciesBySize` (e.g. `getFMort(sim)`) → `plot()` shows one time
  slice (control with `time =`).

```r
plot(getSSB(sim))              # time series per species
plot(getFeedingLevel(params))  # feeding level vs size
```

Compose and enrich these plots:

- **`addPlot(p, arr)`** — overlay another compatible array on an existing ggplot.
- **`plotHover(getBiomass(sim))`** — hover-enabled plotly version of any array
  plot (the compositional API has no `plotly...()` wrappers; use `plotHover()`).
- **`animate(sim)`** or `animate(getFMort(sim))` — animate spectra/rate arrays
  through time (`animateSpectra()` is a retained alias).

## 6. Comparing two objects

To compare two compatible arrays, simulations, or params objects:

- **`plot2(x, y, name1, name2)`** — two compatible arrays on one plot; colour =
  species, linetype = which object. Works for any array class, e.g.
  `plot2(getFMort(params1), getFMort(params2), "Before", "After")`.
- **`plotRelative(x, y)`** — plots the relative difference
  `2 (y - x) / (x + y)` between two compatible arrays. Good for "how much did
  this rate change" questions; supports `wlim`, `species`, `log_x`.
- **`plotSpectra2(object1, object2, name1, name2)`** — compare abundance spectra
  from two sims/params.
- **`plotSpectraRelative(object1, object2)`** — relative difference of two spectra.
- **`plotCDF2(object1, object2, name1, name2)`** — compare cumulative
  distributions from two sims/params.

## 7. Cumulative distributions

**`plotCDF(object, species, power, normalise)`** plots cumulative abundance or
biomass over size — steadier than a density spectrum for eyeballing where biomass
sits. `power = 1` (default) gives biomass, `power = 0` gives numbers; `normalise
= FALSE` plots the cumulative total rather than the proportion.

```r
plotCDF(NS_params, species = c("Cod", "Herring"))
plotCDF(NS_sim, power = 0, normalise = FALSE)
```

## 8. The plankton resource

Resource-related quantities come back as an **`ArrayResourceBySize`** — a numeric
vector over the size grid carrying a `value_name`, `units`, and its `params`, with
`print()`, `summary()`, `as.data.frame()`, and `plot()` methods. Producers include
`NResource(params)` / `finalNResource(sim)`, `getResourceMort(params)`,
`resource_rate(params)` (intrinsic birth rate), `resource_capacity(params)`
(carrying capacity), and `resource_level(params)`.

```r
plot(getResourceMort(params))   # resource mortality vs size
summary(NResource(params))
```

Time-resolved resource data (`NResource(sim)`) is an
`ArrayTimeByResourceBySize`, which `animate()` can play through time.

To include the resource in a species spectrum plot, pass `resource = TRUE`
(supported by `plotSpectra()`, `plotCDF()`, etc.).

## Tips

- Every static plot is a ggplot — add titles/themes with `+ ggtitle(...)` etc.
- Species line colours/types come from the `linecolour`/`linetype` slots of the
  `MizerParams`; override there for consistent styling across plots.
- For interactive exploration prefer the `plotly...()` twin of a named function,
  or `plotHover()` for the compositional array plots.

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
