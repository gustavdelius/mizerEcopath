# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```r
devtools::load_all()        # Load package for interactive use
devtools::test()            # Run all tests
devtools::check()           # Full R CMD check
roxygen2::roxygenise()      # Regenerate documentation from roxygen comments
lintr::lint_package()       # Lint (4-space indentation enforced via .lintr)
```

Run a single test file:
```r
testthat::test_file("tests/testthat/test-matchCatch.R")
```

Run tests matching a pattern:
```r
devtools::test(filter = "matchCatch")
```

## Architecture

**mizerEcopath** bridges [Ecopath](https://ecopath.org/) (biomass-pool ecosystem models) and [mizer](https://sizespectrum.org/mizer/) (size-spectrum models). The goal is to initialise and calibrate mizer models using Ecopath parameters (B, Q, P, C, DC) by matching the size-integrated mizer quantities to their Ecopath counterparts.

### Core workflow

1. **`newAllometricParams()`** — create an allometric (power-law) MizerParams object as starting point
2. **`addEcopathParams()`** — attach Ecopath-derived parameters to species_params
3. **`matchEcopath()`** — master calibration: calls `matchConsumption()` then `matchProduction()` iteratively
4. **`matchCatch()`** — adjust gear selectivity and catchability to match observed catch size distributions (uses TMB for optimisation, see `src/objective_function.cpp`)
5. **`matchDiet()`** — set predation preferences to match diet matrix
6. **`alignResource()`** — scale resource spectrum to match fish abundance
7. **`tuneEcopath()`** — interactive Shiny gadget for manual fine-tuning

### Key design patterns

- Most functions take a `MizerParams` object and return a modified `MizerParams` (functional style).
- "Ecopath quantities" (B, Q/B, P/B, catch) are computed by integrating over the mizer size spectrum — see `getConsumption()`, `getSomaticProduction()`, `getGonadicProduction()`, `getProduction()`, `getM0B()`, `getM2B()`, `getZB()`, `getDietMatrix()`.
- `matchExtMortOnce()`, `matchGonadicProportionOnce()`, `matchRespirationOnce()` are single-step helpers called iteratively by the higher-level matching functions.
- `bindParams()` combines multiple single-species MizerParams into a multi-species model.

### Diffusion / age-structured module

A separate mathematical framework (PDE-based) models cohort growth diffusion for fitting age-at-size data:
- `project_diffusion()` → `solve_diffusion_pde()` → `solve_double_sweep()` (numerical PDE solver)
- `simulateCohort()` for single-cohort dynamics
- `plotAgeLikelihood()`, `getLogLik()` for age-data likelihood fitting
- Supporting biology: `age_mat_vB()`, `von_mises_pdf()`, `spawning_density()`, `length_rebinning_matrix()`

### Shiny tuning app (`tuneEcopath()`)

The interactive gadget is built from composable pieces in `R/`:
- `*Tab.R` files define UI panels (growthTab, dietTab, deathTab, spectraTab, surveyTab, ageTab, catchTab, reproTab)
- `*Control.R` files define control widgets (growthControl, matchControl, exponentControl, fishingControl, diffusionControl, otherControl, reproductionControl)

### C++ optimisation

`src/objective_function.cpp` implements the TMB objective function used by `matchCatch()` to optimise gear selectivity curves against observed catch size distributions.

### Test data

Example datasets in `data/` used in tests:
- `celtic_params.rda`, `celtic_catch.rda` — Celtic Sea 12-species model
- `ns_3_spp_model_*.rda` — North Sea 3-species variants
- `species_params_example.rda`, `ecopath_diet_example.rda`, `diet_matrix_example.rda`

### Documentation

Vignettes are in `vignettes/` as Quarto (`.qmd`) files. The most important is `ecopath_to_mizer.qmd`, which contains the full mathematical derivation of how Ecopath quantities map to mizer integrals. The pkgdown site is pre-built in `docs/`.
