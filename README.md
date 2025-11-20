# mizerEcopath: Extending Ecopath to Mizer

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The `mizerEcopath` package provides a framework for setting up and calibrating multi-species size spectrum models using the [mizer](https://sizespectrum.org/mizer/) package. The size spectra are obtained by balancing the effects of growth and mortality. In that sense it is similar to the [Ecopath](https://ecopath.org/) model, but where Ecopath is only concerned with balancing at the species level, mizerEcopath balances at each size class of each species.

mizerEcopath offers tools to assist in translating empirical observations into steady-state size spectrum models that are ecologically coherent. If one already has an Ecopath model, it is possible to start from parameters of that model.

## Goal

The goal is to make the process of building, tuning, and validating mizer models more accessible, robust, and transparent. These methods support applications ranging from academic research to exploratory ecosystem analysis and policy-relevant fisheries modelling.

## Installation

You can install the development version of mizerEcopath from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("gustavdelius/mizerEcopath")
```

## Workflow

The typical workflow for creating a mizer model from Ecopath data involves the following steps:

1.  **Export Ecopath Parameters**: Open your Ecopath model in "Ecopath with Ecosim" and export the parameters to CSV files.
2.  **Prepare Species Table**: Create a table in R with one row for each species you want to include in the mizer model. This table must contain:
    *   `species`: The name of the species.
    *   `w_max`: The largest observed size of the species.
    *   `w_mat`: The weight at maturity.
    *   `ecopath_groups`: The corresponding Ecopath group(s).
3.  **Initialize Model**: Use `newAllometricParams()` to create a mizer model based on these parameters. This sets up a model with allometric encounter and mortality rates, initially with no interactions.
    ```r
    params <- newAllometricParams(species_params)
    ```
4.  **Match Observations**: Use various functions provided by the package to `matchBiomass()` to match the total biomass for each species, `matchGrowth()` to match the growth rates, `matchCatch()` to match the total yield and the size distribution of the catch, and `matchConsumption()` to match the total consumption for each species.

5.  **Tune the Model**: Use the interactive Shiny gadget `tuneEcopath()` to fine-tune the model to reproduce the observed size structure and other Ecopath parameters.
    ```r
    params <- tuneEcopath(params)
    ```
    This gadget allows you to interactively adjust parameters and see the effects on growth, yield, and consumption in real-time.

6.  **Set predation preferences**: Use `matchDiet()` to set predation preferences in the model to match the observed diet compositions.
    ```r
    params <- matchDiet(params)
    ```

## Background

Ecopath models are based on biomass pools, while mizer models are based on size spectra. `mizerEcopath` bridges this gap by using the Ecopath parameters (Biomass, Consumption, Production, Catch) to constrain and inform the size-dependent parameters required by mizer.

For a detailed explanation of the theoretical background and the translation of parameters, see the vignette `vignettes/ecopath_to_mizer.qmd`.
