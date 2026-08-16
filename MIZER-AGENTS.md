# Mizer

Mizer is an R package for dynamic multi-species size-spectrum modelling of
fish communities. It tracks the full size distribution of each species and
the plankton resource, computing growth, predation, and mortality from
individual-level physiology.

This card is a routing table, not a reference. It says where the answer lives —
in a task skill, or in the installed mizer's help pages — and deliberately does
not restate the answer, because a summary you can act on without opening the
skill is a summary you will act on while it is out of date.

## Do not write mizer code from memory

Mizer's API has evolved, and most mizer code in your training data predates
the version installed here. Recollection that feels solid is often a version
or two stale, and outdated calls frequently run and return plausible numbers
while doing the wrong thing.

Before calling any function you have not looked up in this session, verify its
signature in the installed mizer's documentation using the lookup tools below.
The bundled API index names available functions and the task skills describe
workflows, but neither provides argument lists.

The one correction worth making before reading anything else: **`w_inf`**,
**`w_repro_max`** and **`w_max`** are three distinct parameters, not three names
for the maximum size. The `build-model` skill has the difference.

If the installed mizer disagrees with this card or any skill, the installed
package wins. Report any discrepancy to the user rather than quietly working
around it.

## Key objects

**`MizerParams`** holds all model parameters; **`MizerSim`** is what `project()`
returns. Never modify slots directly. Every `new…`/`set…`/`match…`/`calibrate…`
function returns a *new* object, so always reassign — `params <- setFishing(params, ...)`.

## Task skills (read on demand)

Step-by-step guides for common mizer tasks are installed under `.claude/skills/<name>/SKILL.md`. Claude Code loads them automatically; other agents should **read the matching file before starting** such a task rather than working from memory. They are the reference this card is not: workflows, argument tables and the failure modes of each step. Triggers:

- **`analyse-and-plot`**: Analyse and visualise the results of a mizer simulation or the state of a MizerParams object. Use whenever the user wants to extract, summarise or plot size spectra, biomass, numbers, yield, SSB, feeding level, mortality, diet, trophic level, community indicators, growth curves or the resource — including comparing two models and animating spectra through time. Also covers choosing what a density plot shows (biomass, per_log_size, size_axis), the ArraySpeciesBySize/ArrayTimeBySpecies wrappers the getters return, and writing a custom indicator with sizeIntegral(). Prefer these functions over hand-written array wrangling or custom ggplot code. If a plotting argument that used to work now errors (power=), see the upgrade-mizer-code skill.
- **`analyse-stability`**: Analyse the dynamic stability of a mizer steady state and characterise limit cycles (experimental). Use whenever the user asks whether a steady state is stable or unstable, wants the spectral radius or leading eigenvalue, the period of an emergent oscillation, a Hopf bifurcation, a limit cycle to build or plot, or a bifurcation diagram over fishing effort — via getStability(), steadyNewton(stability = TRUE), getLimitCycleSim() and plotBifurcation(). This skill and calibrate-model share steadyNewton() and getSteadyResidual(): use calibrate-model to find a steady state, this skill to ask whether the state you found is stable. Assumes the standard semichemostat resource dynamics.
- **`build-model`**: Build a new mizer model from a species-parameter data frame. Use whenever the user wants to create a MizerParams object with newMultispeciesParams() (or newTraitParams, newCommunityParams, newSingleSpeciesParams), decide which species-parameter columns to supply and which to leave to mizer's allometric defaults, set up an interaction matrix, choose the size grid (no_w, min_w, max_w, min_w_pp), or save and reload a finished model. Follow this ordered workflow rather than guessing at parameters or writing the dynamics by hand. To change a model that already exists see the change-parameters skill; fishing is covered by the set-up-fishing skill and steady state and calibration by the calibrate-model skill.
- **`calibrate-model`**: Bring a mizer model to steady state and calibrate it to observed data. Use whenever the user wants to find the steady state (steady, projectToSteady, steadySingleSpecies, steadyNewton), match modelled biomass, numbers, yield or growth to observations (calibrateBiomass, matchBiomasses, calibrateNumber, matchNumbers, matchGrowth), supply those observations (the biomass_observed/biomass_cutoff or number_observed species-parameter columns, the yield_observed gear-parameter column), set the level of density-dependent reproduction (reproduction_level<-), check convergence (getSteadyResidual), or diagnose a model that collapses, explodes or will not settle. To ask whether the steady state you found is dynamically stable, see the analyse-stability skill.
- **`change-parameters`**: Change parameters of an existing mizer model correctly, so that the change propagates downwards and is not silently overwritten. Use whenever the user wants to modify species parameters, size-dependent rates, the resource or the interaction matrix — and especially when unsure which level to work at: given_species_params() vs species_params(), or changing a species parameter vs setting a rate array directly (setSearchVolume, setPredKernel, setParams…). Covers which values get recalculated and which stay put, the freeze trap when a rate array is set by hand, length-vs-weight precedence (l_mat vs w_mat), resource balancing (balance =), and warnings that a change could not take effect. Fishing gears are covered by the set-up-fishing skill, custom rate functions by the extend-mizer skill.
- **`extend-mizer`**: Extend or customise mizer's dynamics — add external food or mortality, replace a built-in rate calculation, or add a new ecosystem component. Use whenever the user wants a custom encounter/growth/mortality/reproduction formulation (setRateFunction() wrapping or replacing mizerEncounter, mizerPredRate, mizerMort, mizerEReproAndGrowth), a background food, diffusion or predation source needing no new state variable (setExtEncounter, setExtMort, setExtDiffusion), a new dynamical pool such as detritus or carrion (setComponent), or asks how to make mizer do something its standard setters do not cover — including how a custom rate must respect the second_order_w quadrature scheme. Pick the lightest mechanism that works: to change an existing rate's parameters see the change-parameters skill.
- **`run-simulation`**: Project a mizer model forward in time and set up fishing-effort scenarios. Use whenever the user wants to run a simulation with project(), choose the time stepping (t_max, dt, t_save, t_start), give constant, time-varying or per-gear fishing effort, continue an existing MizerSim (append = TRUE), carry a simulation's end state into a new run (setInitialValues, finalParams), run to a new steady state after a change, or set up scenario comparisons — including what to do about numerical diffusion and the second_order_w scheme when growth looks smeared. For extracting and plotting the results see the analyse-and-plot skill; for reaching steady state first see the calibrate-model skill.
- **`set-up-fishing`**: Set up or change fishing in a mizer model — gears, selectivity curves, catchability and effort. Use whenever the user wants to define fishing gears in the gear_params data frame (sel_func, catchability, yield_observed), choose and parameterise a selectivity function (knife_edge with knife_edge_size, knife_edge_length, sigmoid_length and double_sigmoid_length with l50/l25/l50_right/l25_right, sigmoid_weight), set which gear catches which species, apply the result with setFishing(), or set the fishing effort (initial_effort<-). Note that catchability sets the units in which effort is measured. Matching modelled yields to observed ones is covered by the calibrate-model skill.
- **`understand-size-spectrum-dynamics`**: Understand how mizer models behave: which quantities you impose and which the model produces for itself, the food and predation feedback loops that couple species, what sets the slope of the steady-state spectrum and the timescale of its dynamics, and the two distinct kinds of density dependence. Use whenever reasoning about why a species collapses, explodes, or oscillates, why growth or mortality is not what you asked for, why changing one species moved another, why a model is insensitive to fishing, or what the Sheldon spectrum, feeding level, reproduction level (R_max) and trophic cascades mean in a mizer model.
- **`upgrade-mizer-code`**: Diagnose and fix user code or models that broke, changed results, or started warning after a mizer upgrade — every documented change from 2.5.4 through 3.3, release by release, with the fix. Use whenever a script that "used to work" now errors, a deprecation warning appears, plots or numbers differ from a previous run, an argument is suddenly unused or rejected (power=, sim=, time_range=, setParams() rejecting an argument it does not use), a function has gone (matchYields, calibrateYield), a parameter change warns that it cannot take effect, an identical() comparison against a saved rate array fails, or `$` on a parameter table stopped matching partially. Starts from a symptom index, so search it by the message the user actually saw.

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
  returns the document open in the editor, which is usually the thing a
  vague request refers to. It works only when the user is in RStudio or
  Positron; in any other editor it errors, which tells you nothing about
  the rest of the session. Ask them which file they mean instead.
- **Run mizer code in their session, not in a scratch script.**
  `btw_tool_run_r` evaluates in their global environment, so results
  persist and the user can carry on with them. Plots come back to you as
  images: after calibrating or projecting, plot the result and *look at it*
  before reporting success.
- **That session holds their work.** You are in *their* global
  environment: assigning over an existing object destroys it, and `load()`
  and `rm()` reach everything they have open. There is no undo and nothing
  warns you. Assign to a new name, or say what you are about to overwrite
  first, and do exploratory work in a throwaway environment or a separate
  `Rscript` -- keep the session itself for the steps whose result the user
  wants to keep. Say so whenever you do change a global. Long projections
  block their console, so keep `t_max` modest unless asked otherwise.

If these tools are missing, the user has not connected their session; ask
them to run `mizerAgents::connect_mizer_agent()` in their R console, rather
than working around it.

## Finding the right mizer function

Two steps, and they use different sources:

1. **Which function do I need?** Grep the bundled API index: every
   exported function with a one-line description, grouped by workflow
   stage (creating a model, tuning the steady state, projecting,
   plotting). Grep it for a keyword; do not read the whole file:
   /home/gustav/R/x86_64-pc-linux-gnu-library/4.6/mizer/llms.txt
2. **How do I call it?** Read the help page for the *installed* mizer with
   `btw_tool_docs_help_page`. The index above deliberately carries no
   argument lists, and this card is not a reference either.

Never supply arguments from memory for a function you have not looked up
in this session. Rendered documentation for the current release is at
<https://sizespectrum.org/mizer/reference/>, but the installed version is
what your code will run against, so prefer the local help page.
