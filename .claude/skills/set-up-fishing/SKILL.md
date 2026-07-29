---
name: set-up-fishing
description: >-
  Set up or change fishing in a mizer model — gears, selectivity curves,
  catchability, and effort. Use whenever the user wants to define fishing gears,
  choose or configure a selectivity function (knife_edge, sigmoid_length,
  double_sigmoid_length, sigmoid_weight), set which gear catches which species,
  change catchability, or set the baseline fishing effort with setFishing() and
  the gear_params data frame.
---

# Setting up fishing

Fishing mortality in mizer is `F[g,i](w) = S[g,i](w) * Q[g,i] * E[g]` — a gear's
**selectivity** at size, times its **catchability** for a species, times its
**effort**. You configure the first two through the `gear_params` data frame and
`setFishing()`; effort at run time is handled in the `run-simulation` skill.

Setter functions return a new `MizerParams` — always reassign.

## The `gear_params` data frame

One row per **gear–species combination** (a gear that catches three species
contributes three rows). Required columns:

| Column | Meaning |
|---|---|
| `gear` | gear name |
| `species` | species this row applies to |
| `sel_func` | name of the selectivity function (e.g. `"sigmoid_length"`) |
| `catchability` | scales `F` for this gear–species pair (default 1) |

Plus **one column per parameter of the chosen `sel_func`**, named exactly like
the function's argument (see below). Read/replace the table with
`gear_params(params)` / `gear_params(params) <- ...`.

```r
gp <- gear_params(params)
gp["Cod, Otter", "catchability"] <- 0.8   # row names are "species, gear"
gear_params(params) <- gp
```

**Setting one up from scratch.** Assign a fresh data frame with one row per
gear–species combination. Only `species` is strictly required: `gear` defaults to
the species name, `sel_func` to `knife_edge`, `catchability` to 1, and the
`knife_edge` cut-off to `w_mat`. You must, however, supply the parameter columns
of whatever `sel_func` you choose. mizer generates the `"species, gear"` row names
for you:

```r
gear_params(params) <- data.frame(
    gear         = c("Otter", "Beam"),
    species      = c("Cod",   "Cod"),
    sel_func     = "sigmoid_length",   # recycled to both rows
    l50          = c(25, 20),          # 50% selected at this length (cm)
    l25          = c(20, 15),          # 25% selected at this length (cm)
    catchability = 1
)
```

If each species is caught by only one gear, you may instead put the gear columns
directly in `species_params` when building the model; mizer copies them into
`gear_params`. (Later edits to those `species_params` columns will **not**
propagate — edit `gear_params` after construction.)

## Selectivity functions

Every selectivity function takes `w` as its first argument and returns a value in
`[0, 1]` at each size. The other arguments must appear as columns in
`gear_params`:

| `sel_func` | Parameter column(s) | Shape |
|---|---|---|
| `knife_edge` (default) | `knife_edge_size` | 0 below the size, 1 above (default size = `w_mat`) |
| `sigmoid_length` | `l50`, `l25` | smooth; lengths (cm) at 50% and 25% selection |
| `double_sigmoid_length` | `l50`, `l25`, `l50_right`, `l25_right` | dome-shaped (selects a length band) |
| `sigmoid_weight` | `sigmoidal_weight`, `sigmoidal_sigma` | smooth transition in weight |

`sigmoid_length` is the most commonly used. You can also write your own function
(first argument `w`, returns selectivity at size) and name it in `sel_func`.

```r
gp <- gear_params(params)
gp$sel_func <- "sigmoid_length"
gp$l50 <- 25            # 50% selected at 25 cm
gp$l25 <- 20            # 25% selected at 20 cm
gear_params(params) <- gp
```

## The selectivity and catchability arrays

Behind the scenes mizer turns the `gear_params` table into two numeric arrays,
the ones that enter the fishing-mortality formula directly. You can read them,
and — when a `sel_func` cannot express the shape you need — set them by hand:

| Function | Returns | Dimensions |
|---|---|---|
| `catchability(params)` / `getCatchability(params)` | `Q[g,i]` | gear × species |
| `selectivity(params)` / `getSelectivity(params)` | `S[g,i](w)` | gear × species × size |

The bare and `get`-prefixed names are equivalent. Each has a matching setter that
pushes an array straight into the model (this routes through `setFishing()`, so
validation still runs):

```r
selectivity(params)["Otter", "Cod", ]      # the S curve for one gear–species pair

sel <- getSelectivity(params)
sel["Otter", "Cod", ] <- my_curve          # length = number of size bins, in [0, 1]
selectivity(params) <- sel                 # triggers recalculation via setFishing()
```

Setting an array by hand **freezes** it: mizer marks it manual and stops
recalculating it from `gear_params`, so later edits to the gear table leave it
untouched (you'll see a message saying so). To discard the hand-set array and
rebuild from `gear_params`, call `setFishing(params, reset = TRUE)`.

## `setFishing()`

`setFishing(params, selectivity = NULL, catchability = NULL, reset = FALSE,
initial_effort = NULL, ...)` recomputes the fishing setup after a change.
Assigning `gear_params(params) <- ...` already triggers recalculation, so you
usually only call `setFishing()` directly when supplying a `selectivity` or
`catchability` **array**, setting a baseline effort, or rebuilding from scratch:

```r
params <- setFishing(params, initial_effort = c(Otter = 1, Beam = 0.5))
params <- setFishing(params, reset = TRUE)   # rebuild arrays from gear_params
```

## Baseline effort

The model stores a baseline effort per gear, read with `initial_effort(params)`
(a named vector) and set via `initial_effort(params) <- ...` or the
`initial_effort` argument of `setFishing()`. This is the effort used when
`project()` is called without an explicit `effort` argument.

```r
initial_effort(params) <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
```

To vary effort **through time** in a run, pass a time × gear array to
`project()` — see the `run-simulation` skill.

## Catchability sets the units of effort

Because `F = S * Q * E`, only the product `Q * E` is pinned down by the fishing
mortality — so **catchability defines what one unit of effort means**. This gives
you two common conventions:

- **Effort = fishing mortality rate.** Set `catchability = 1`; then an effort of
  `E` produces `F = E` on fully-selected sizes. Simplest when you want to drive
  the model directly with F values (e.g. from a stock assessment).
- **Effort in real-world units** (vessel-days, kW-days, …). Fold the conversion
  into catchability: `Q` is the fraction of the selected stock taken per unit of
  whatever effort you supply. A useful trick is to set each species' catchability
  to its `F` in a chosen reference year, so an effort of 1 reproduces that year's
  fishing mortality and other years' efforts are relative to it.

Either way, if you rescale catchability you must rescale effort inversely to keep
the same `F`. This is why yield calibration (`calibrateYield`) depends on the
fishing setup being fixed first.

## Checking the setup

```r
gear_params(params)                 # the gear table
catchability(params)                # Q array (gear × species)
selectivity(params)                 # S array (gear × species × size)
initial_effort(params)              # baseline effort per gear
plotFMort(params)                   # realised fishing mortality at size
getFMortGear(params)                # F by gear × species × size
```

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
