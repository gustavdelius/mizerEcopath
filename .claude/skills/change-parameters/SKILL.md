---
name: change-parameters
description: >-
  Change parameters of an existing mizer model correctly. Use whenever the user
  wants to modify species parameters, size-dependent rates, fishing, the
  resource, or interactions — and especially when unsure which accessor to use:
  given_species_params() vs species_params(), changing a species parameter vs
  setting a rate array directly (setSearchVolume, setPredKernel, setParams…), or
  gear_params() vs the resource setters. Follow these rules to avoid changes
  that silently fail to propagate or get overwritten.
---

# Changing model parameters

Guiding principle: **change the model at the highest level that expresses the
intent, and let mizer propagate the change downwards.** Drop to a lower level
only to deliberately override mizer's calculation. Every setter returns a **new**
`MizerParams` — always reassign (`params <- setResource(params, ...)`).

A mizer model is built in layers, and almost every change is a choice of which
layer to reach into:

1. **Size-independent parameters** — the high-level inputs: per-species scalars
   (`w_inf`, `beta`, `gamma`, `h`, `erepro`, …), fishing gears (`gear_params`),
   resource scalars (`resource_params`), and the interaction matrix. Most are used
   to *calculate* the size-dependent rate arrays below. **This is where almost all
   everyday work happens.**
2. **Size-dependent rates** — arrays over size (search volume, metabolic rate,
   predation kernel, selectivity, resource capacity, …) that mizer builds from the
   level-1 parameters. Reach in here only when you need a size-dependence mizer
   does not produce by default.
3. **Rate functions** — the functions mizer calls to compute rates during a
   simulation. Replace one with `setRateFunction()` (e.g. to make a rate
   time-dependent), or add a whole new dynamical component with `setComponent()`.
   See the `extend-mizer` skill.

## Species parameters: which accessor

| Accessor | Returns |
|---|---|
| `given_species_params(params)` | only the parameters the user supplied |
| `calculated_species_params(params)` | parameters mizer derived or defaulted |
| `species_params(params)` | everything (given, with calculated filling gaps) |

**Rule (mizer ≥ 3.2): change species parameters with `species_params(params) <-`.**
As of mizer 3.2 this is the recommended setter for scripts: it detects what you
changed, records it as *given* (so defaults can no longer overwrite it), and
triggers recalculation of the derived scalars **and** the size-dependent rate
arrays that depend on it.

```r
species_params(params)$beta <- 150   # recorded as given; also rebuilds the predation kernel
```

`given_species_params(params) <-` does the same recording and recalculation and
is preferable in **interactive** sessions, because it additionally *warns* when
you change a parameter whose effect is overridden by another parameter you have
already given.

> **Version note.** Older guidance said to avoid `species_params(params) <-`
> because it bypassed the `given_species_params` protection and skipped
> recalculation. That was fixed in mizer 3.2. On mizer **< 3.2**, still prefer
> `given_species_params(params) <-` for edits.

Columns come back as named vectors: `species_params(params)$w_mat` is named by
species; `given_species_params(params)$gamma` is `NA` where the user never set it.

## How a species-parameter change propagates

Many species parameters exist only to set up a rate array; changing one re-runs
the relevant setter automatically:

| Species parameter(s) | Sets up | via |
|---|---|---|
| `gamma`, `q` | search volume | `setSearchVolume()` |
| `h`, `n` | maximum intake rate | `setMaxIntakeRate()` |
| `k`, `ks`, `p` | metabolic rate | `setMetabolicRate()` |
| `z0`, `z_ext`, `d` | external mortality | `setExtMort()` |
| `beta`, `sigma`, `pred_kernel_type` | predation kernel | `setPredKernel()` |
| `w_mat`, `w_mat25`, `w_repro_max`, `m` | reproduction allocation | `setReproduction()` |

Other species parameters are used **directly** and build no array (changing them
just changes the model): `alpha`, `w_min` (egg size), `erepro`, `R_max`,
`interaction_resource`, and the length–weight parameters `a`, `b`.

## The feeding level is set by `f0`, not `h`

When `gamma` (the search-volume coefficient) is not supplied, mizer **derives it
from the target feeding level `f0`** — roughly `gamma ∝ h · f0 / (1 - f0)` — so
that the modelled steady-state feeding level equals `f0`. A consequence often
missed: raising `h` does **not** lower the feeding level, because the derived
`gamma` compensates and the feeding level stays pinned at `f0`. To make growth
more or less resource-dependent (a stronger or weaker density dependence / the
"phantom-jam" feedback), change **`f0`**, not `h`: low `f0` makes juvenile growth
strongly resource-limited, `f0` near 1 makes it nearly saturated and
resource-insensitive.

Corollary: `h = Inf` (a deliberately "no-satiation" model) makes the derived
`gamma` non-finite and throws `search_vol must not contain non-finite values` —
supply `gamma` explicitly in that case.

## Setting a rate array directly — and the freeze trap

Each size-dependent rate has its own setter. Called with only `params`, it
**recomputes** the array from the current species parameters (unless the array is
frozen — see below):

```r
params <- setMetabolicRate(params)         # recompute from k, ks, p
params <- setParams(params)                # rebuild ALL rate arrays at once
```

Each rate array also has a **direct setter/getter** — `metab(params) <-`,
`search_vol(params) <-`, etc. (and `metab(params)`, `search_vol(params)` to read
it) — that writes an array straight in. This is equivalent to passing the array
to the `set…()` function, except the direct setter modifies `params` **in place**
while `set…()` returns a **new** object:

```r
metab(params) <- my_array                             # direct setter, in place
params <- setMetabolicRate(params, metab = my_array)  # same, returns new object
```

**Trap:** either way, supplying an array **freezes** it. mizer marks it "set
manually" and will no longer update it from species parameters. After this,
changing the feeding species parameter has *no effect* on that rate:

```r
params <- setSearchVolume(params, search_vol = my_array)  # frozen
given_species_params(params)$gamma <- 2 * gamma           # search volume UNCHANGED now
```

To hand control back to mizer you must pass **`reset = TRUE`**. A bare
`set…(params)` will **not** recompute a frozen array — it leaves the manual value
in place and warns:

```r
params <- setSearchVolume(params, reset = TRUE)   # drop the override, recompute
```

If a user reports "I changed `gamma`/`h`/`beta` but nothing happened," suspect a
manually-set (frozen) rate array — recompute it with `set…(params, reset = TRUE)`.
The `reset` argument works the same way for every rate setter, including
`setFishing()` and `setResource()`.

## Fishing

Gears, selectivity, and catchability live in `gear_params(params)`, one row per
gear–species pair (row names `"species, gear"`). Assigning to it **does**
recompute the fishing arrays.

```r
gp <- gear_params(params)
gp["Cod, Otter", "catchability"] <- 0.8
gear_params(params) <- gp                  # rebuilds fishing mortality
```

Use `setFishing()` for supplying selectivity/catchability arrays directly or
setting baseline effort. See the `set-up-fishing` skill.

## Resource — same model as species parameters

`resource_params(params)` returns scalars (`kappa`, `lambda`, `r_pp`, `n`,
`w_pp_cutoff`) that set up the resource size-spectrum arrays, and — like
`given_species_params<-` — assigning to it **rebuilds those arrays**:

```r
resource_params(params)$kappa  <- 1e11     # rebuilds the carrying capacity (cc_pp)
resource_params(params)$lambda <- 2.05     # rebuilds cc_pp (slope)
resource_params(params)$r_pp   <- 10       # rebuilds the replenishment rate (rr_pp)
```

`kappa`, `lambda`, `w_pp_cutoff` drive the capacity `cc_pp`; `r_pp`, `n` drive
the rate `rr_pp`.

Setting the size-resolved arrays directly **freezes** them ("set manually"), the
same as the species rate arrays — a later scalar change then leaves them alone.
Recompute from the scalars with `setResource(params, reset = TRUE)`:

```r
resource_capacity(params) <- my_capacity   # array over size; now frozen
resource_params(params)$kappa <- 1e11      # ignored while cc_pp is frozen
params <- setResource(params, reset = TRUE)  # unfreeze: recompute from scalars
resource_rate(params)  <- my_rate          # freeze rr_pp instead
resource_level(params) <- my_level
```

## Interaction matrix

```r
inter <- getInteraction(params)
inter["Cod", "Herring"] <- 0.5
params <- setInteraction(params, inter)
```

The resource interaction is the `interaction_resource` species-parameter column
instead.

## Which to use

| To change… | Use |
|---|---|
| a per-species value | `species_params(params) <- …` (mizer ≥ 3.2; `given_species_params(params) <-` interactively or on older mizer) |
| a rate, keeping it tied to the parameters | change the underlying species parameter |
| a rate to a bespoke array (freezing it) | `metab(params) <- …` (direct) or the matching `set…(params, array)` |
| a frozen rate back to its default form | `set…(params, reset = TRUE)` |
| everything after several edits | `setParams(params)` |
| fishing gears / selectivity / catchability | `gear_params(params) <- …` |
| baseline effort or selectivity/catchability arrays | `setFishing(params, …)` |
| the resource (kappa, lambda, r_pp, …) | `resource_params(params) <- …` |
| the resource capacity/rate as a bespoke array (freezing it) | `resource_capacity(params) <- …` / `resource_rate(params) <- …` / `setResource(params, …)` |
| species interactions | `setInteraction(params, …)` |
| how a rate is *computed* (e.g. to make it time-dependent) | `setRateFunction(params, …)` (see `extend-mizer`) |

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
