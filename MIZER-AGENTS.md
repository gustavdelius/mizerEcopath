# Mizer

Mizer is an R package for dynamic multi-species size-spectrum modelling of
fish communities. It tracks the full size distribution of each species and
the plankton resource, computing growth, predation, and mortality from
individual-level physiology.

## Do not write mizer code from memory

Mizer's API has moved on, and most mizer code in your training data predates
the version installed here — recollection that feels solid is often a version
or two stale. Before calling any function you have not looked up in this
session, read its help page from the installed mizer and check the real
signature. This failure is quiet: outdated calls often still run and return
plausible numbers.

Nothing in this repository is a substitute for that. The bundled API index
(path at the end of this file) tells you which functions exist, not how to call
them, and this card is a summary rather than a reference. Argument lists come
from the installed package or they come from a guess.

Memory is most often stale on:

- **maximum size** — `w_inf`, `w_repro_max` and `w_max` are three distinct
  parameters (see below), not one
- **setting species parameters** — use `species_params(params) <- value`, which
  records the change and recalculates dependent quantities
- **reproduction** — `setBevertonHolt()` takes `erepro`, `R_max` *or*
  `reproduction_level`; check which the task calls for

If the installed mizer disagrees with this file, the installed mizer wins.
Check with `?name` and say so rather than quietly working around it.

## Core workflow

```r
library(mizer)

# 1. Create model parameters from a species data frame
params <- newMultispeciesParams(species_params, interaction)

# 2. Find the steady state (sets initial values)
params <- steady(params)

# 3. Calibrate to observed biomasses / yields
params <- calibrateBiomass(params)  # adjusts kappa
params <- matchBiomasses(params)    # adjusts R_max per species
params <- matchGrowth(params)       # adjusts h per species

# 4. Tune density-dependent reproduction
params <- setBevertonHolt(params, reproduction_level = 0.25)

# 5. Project forward in time
sim <- project(params, t_max = 20, effort = 1)

# 6. Analyse results
plot(sim)
getBiomass(sim)
getYield(sim)
plotSpectra(sim)
```

## Key objects

**`MizerParams`** — holds all model parameters. Never modify slots directly.
All setter functions return a new copy: `params <- setFishing(params, ...)`.

**`MizerSim`** — simulation output from `project()`. Arrays: `N(sim)` (time ×
species × size), `NResource(sim)`.

## Species parameters

The `species_params` data frame must have `species` (name) and the
von Bertalanffy asymptotic weight `w_inf`. Everything else has defaults.
Change species parameters with `species_params(params) <- value`, which records
the change and triggers recalculation of dependent quantities. See the
`change-parameters` skill.

| Column | Meaning |
|--------|---------|
| `w_inf` | Von Bertalanffy asymptotic weight (g); accepted maximum-size input, sets `w_repro_max` |
| `w_max` | Computational upper size-grid boundary (g) — purely numerical; defaults to `1.5 * w_inf` |
| `w_repro_max` | Weight beyond which no growth/reproduction |
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default ~100) |
| `sigma` | S.d. of lognormal predation kernel (default ~1.3) |
| `h` | Max intake rate coefficient |
| `alpha` | Assimilation efficiency (default 0.6) |
| `erepro` | Reproductive efficiency |
| `R_max` | Beverton-Holt max recruitment |
| `biomass_observed` | Observed biomass for `calibrateBiomass()` |

## Units

Weights in grams, lengths in cm, time in years.

## Numerical scheme for dynamics

The default `project()` flux scheme (first-order upwind) carries substantial
*numerical* diffusion that silently smears cohorts and travelling waves and can
completely damp real oscillations / limit cycles — a correctness issue, not just
accuracy. For any study of dynamics (oscillations, cohort waves, diffusion),
build the model with `second_order_w = TRUE` (van Leer flux) and project with
`method = "tr_bdf2"` (second order in time). See the `run-simulation` skill.

## Gotchas

- `w_max` defaults to `1.5 * w_inf`. Passing `max_w = w_inf` to
  `newMultispeciesParams()` then errors — set a `w_max` column equal to `w_inf`
  as well.
- The steady-state feeding level is set by the `f0` species parameter (from which
  the default `gamma` is derived), **not** by `h`; `h = Inf` makes `gamma`
  non-finite. See the `change-parameters` skill.
- With growth diffusion on (`D_ext > 0`), set `w_max` well above the sizes you
  analyse so abundance at the boundary stays negligible; the default
  `1.5 * w_inf` is usually enough, raise it if `D_ext` is large.

## Plotting

The return values of most `get...()` functions also have `plot()` methods,
so you can visualise any quantity directly.
Always prefer this over writing custom plotting code.

```r
plot(getSSB(sim))           # ArrayTimeBySpecies  → time series per species
plot(getTrophicLevel(params)) # ArraySpeciesBySize → curve per species
```

In addition, Mizer provides many custom plotting functions. 

```r
plot(sim)              # overview of simulation
plotSpectra(sim)       # size spectra
plotBiomass(sim)       # biomass over time
plotYield(sim)         # yield over time
plotGrowthCurves(sim)  # growth curves
plotFMort(sim)         # fishing mortality
```

Grep for "plot" in the bundled API index (path at the end of this file) to
discover the full list of available plots before writing any custom code, then
read the help page of the one you pick for its arguments.

## Extending mizer

To replace a rate function: `params <- setRateFunction(params, "Encounter", myFun)`.
To add a new ecosystem component: `params <- setComponent(params, "detritus", ...)`.
See https://sizespectrum.org/mizer/articles/extending-mizer.html

## Task skills (read on demand)

Step-by-step guides for common mizer tasks are installed under `.claude/skills/<name>/SKILL.md`. Claude Code loads them automatically; other agents should **read the matching file before starting** such a task rather than working from memory. Triggers:

- **`analyse-and-plot`**: Analyse and visualise the results of a mizer simulation or the state of a MizerParams object. Use whenever the user wants to extract, summarise, or plot size spectra, biomass, yield, SSB, abundance, feeding level, mortality, diet, trophic level, community indicators, growth curves, or the plankton resource — including comparing two simulations or animating spectra through time. Prefer the built-in mizer functions described here over writing custom extraction or ggplot code.
- **`build-multispecies-model`**: Build a calibrated multi-species mizer model from a species-parameter data frame. Use whenever the user wants to create a MizerParams object with newMultispeciesParams(), set up an interaction matrix or fishing gears, bring the model to steady state with steady(), or calibrate/match it to observed biomasses, yields, or growth (calibrateBiomass, matchBiomasses, matchGrowth, calibrateYield, setBevertonHolt). Follow this ordered workflow rather than guessing at parameters or writing the dynamics by hand.
- **`calibrate-model`**: Bring a mizer model to steady state and calibrate it to observed data. Use whenever the user wants to find the steady state (steady, projectToSteady, steadySingleSpecies), match modelled biomass, yield, or growth to observations (calibrateBiomass, matchBiomasses, calibrateYield, matchGrowth), set the level of density-dependent reproduction (setBevertonHolt), or diagnose why a model will not settle or reproduce observed values.
- **`change-parameters`**: Change parameters of an existing mizer model correctly. Use whenever the user wants to modify species parameters, size-dependent rates, fishing, the resource, or interactions — and especially when unsure which accessor to use: given_species_params() vs species_params(), changing a species parameter vs setting a rate array directly (setSearchVolume, setPredKernel, setParams…), or gear_params() vs the resource setters. Follow these rules to avoid changes that silently fail to propagate or get overwritten.
- **`extend-mizer`**: Extend or customise mizer's dynamics — add external food/mortality, replace a built-in rate calculation, or add a new ecosystem component. Use whenever the user wants a custom encounter/growth/mortality/reproduction formulation (setRateFunction), a background food or predation source (setExtEncounter, setExtMort), a new dynamical pool like detritus or carrion (setComponent), or asks how to make mizer do something its standard setters do not cover.
- **`run-simulation`**: Project a mizer model forward in time and set up fishing-effort scenarios. Use whenever the user wants to run a simulation with project(), specify constant or time-varying fishing effort, choose a projection method or time step, run to a new steady state after a change, continue an existing MizerSim, or set up scenario comparisons. For extracting and plotting the results, see the analyse-and-plot skill.
- **`set-up-fishing`**: Set up or change fishing in a mizer model — gears, selectivity curves, catchability, and effort. Use whenever the user wants to define fishing gears, choose or configure a selectivity function (knife_edge, sigmoid_length, double_sigmoid_length, sigmoid_weight), set which gear catches which species, change catchability, or set the baseline fishing effort with setFishing() and the gear_params data frame.

A skill's directory may also hold a `NOTES.md` recording what earlier work in this project found. Read it whenever you read the `SKILL.md`, and treat it as taking precedence. Write new project-specific findings there, creating the file if needed.

Do not edit `SKILL.md` or this card: both are installed by `mizerAgents::setup_mizer_agent()`. Project notes that belong to no single skill go in `AGENTS.md` / `CLAUDE.md`, outside the `<!-- mizerAgents: ... -->` markers. A lesson that is true of mizer in general rather than of this project belongs upstream, where every project gets it: tell the user, and offer to report it at <https://github.com/sizespectrum/mizerAgents/issues>.


## The user's live R session

This project is configured with an MCP server named `r-mizer` (provided by
the [btw](https://posit-dev.github.io/btw/) package) that connects you to the
R session the user is working in. Use it:

- **Look up mizer functions before calling them.** The docs tools
  (`btw_tool_docs_help_page`, `btw_tool_docs_available_vignettes`,
  `btw_tool_docs_vignette`, `btw_tool_docs_package_news`) read the *installed*
  mizer. That is the authority on argument names and defaults, above this
  card, and far above your own recollection, which is very likely stale.
- **Inspect what the user already has.** `btw_tool_env_describe_environment`
  lists the objects in their global environment; do not rebuild a
  `MizerParams` or re-run a simulation that is already sitting there.
- **Read what they are looking at.** `btw_tool_ide_read_current_editor`
  returns the document open in RStudio, which is usually the thing a vague
  request refers to.
- **Run mizer code in their session, not in a scratch script.**
  `btw_tool_run_r` evaluates in their global environment, so results
  persist and the user can carry on with them. Plots come back to you as
  images: after calibrating or projecting, plot the result and *look at it*
  before reporting success.
- **That session holds their work.** Assigning over an existing object
  destroys it, and there is no undo. Assign to a new name, or say what you
  are about to overwrite first. Long projections block their console, so
  keep `t_max` modest unless asked otherwise.
- **This project is an R package**, most likely a mizer extension. Use the
  package tools (`btw_tool_pkg_load_all`, `btw_tool_pkg_document`,
  `btw_tool_pkg_test`, `btw_tool_pkg_check`, `btw_tool_pkg_coverage`) rather
  than shelling out to `R CMD` or `devtools`: they run in the user's session,
  so after `load_all` the new code is live and you can exercise it directly.
  Read the `extend-mizer` skill before changing how a mizer rate is
  calculated.

If these tools are missing, the user has not run `btw::btw_mcp_session()` in
their RStudio console; ask them to, rather than working around it.


## Finding the right mizer function

Two steps, and they use different sources:

1. **Which function do I need?** Grep the bundled API index: every
   exported function with a one-line description, grouped by workflow
   stage (creating a model, tuning the steady state, projecting,
   plotting). Grep it for a keyword; do not read the whole file:
   /home/gustav/R/x86_64-pc-linux-gnu-library/4.6/mizerAgents/llms.txt
2. **How do I call it?** Read the help page for the *installed* mizer with
   `btw_tool_docs_help_page`. The index above deliberately carries no
   argument lists, and this card is not a reference either.

Never supply arguments from memory for a function you have not looked up
in this session. Rendered documentation for the current release is at
<https://sizespectrum.org/mizer/reference/>, but the installed version is
what your code will run against, so prefer the local help page.

