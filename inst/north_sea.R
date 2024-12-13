library(mizerEcopath)
library(rfishbase)
library(dplyr)
library(readr)

## Prepare catch data ----
# We do this first because we will use it to set maximum sizes for species
catch <- readRDS("inst/extdata/ns_catch.rds")
# Select "total" gear only
catch <- catch[catch$gear == "total", ]

## Prepare species parameters ----

sp <- NS_species_params |>
    select(species, w_mat, w_max)

# We will need better weight-length relationship parameters
species_to_latin <- list(
    "Sprat"   = "Sprattus sprattus",
    "Sandeel" = "Ammodytes spp",
    "N.pout"  = "Trisopterus esmarkii",
    "Herring" = "Clupea harengus",
    "Dab"     = "Limanda limanda",
    "Whiting" = "Gadus merlangus",
    "Sole"    = "Solea solea",
    "Gurnard" = c("Triglidae", "Eutrigula gurnardus"),
    "Plaice"  = "Pleuronectes platessa",
    "Haddock" = "Gadus aeglefinus",
    "Cod"     = "Gadus morhua",
    "Saithe"  = "Pollachius virens"
)

sp$a <- 0.01
sp$b <- 3

# Set maximum weights to at least 1.2 times the maximum observed weight
sp <- catch |>
    group_by(species) |>
    summarise(max_observed_length = max(length + dl)) |>
    right_join(sp, by = c("species" = "species")) |>
    mutate(max_observed_weight = a * max_observed_length^b) |>
    mutate(w_max = pmax(w_max, 1.2 * max_observed_weight, na.rm = TRUE))

# Add gonadic proportion
sp$gonad_proportion <- 0.2

## Add Ecopath species parameters ----

basic_estimates <-
    read_csv("inst/extdata/North Sea-Basic estimates.csv")

# Dictionary between species and ecopath groups
species_to_groups <- list(
    "Sprat"   = "Sprat",
    "Sandeel" = "Sandeels",
    "N.pout"  = "Norway pout",
    "Herring" = c("Herring (juvenile 0, 1)", "Herring (adult)"),
    "Dab"     = "Dab",
    "Whiting" = c("Juvenile Whiting (0-1, 0-20cm)", "Whiting (adult)"),
    "Sole"    = "Sole",
    "Gurnard" = "Gurnards",
    "Plaice"  = "Plaice",
    "Haddock" = c("Juvenile Haddock (0-1, 0-20cm)", "Haddock (adult)"),
    "Cod"     = c("Juvenile Cod(0-2, 0-40cm)", "Cod (adult)"),
    "Saithe"  = c("Juvenile Saithe (0-3, 0-40cm)", "Saithe (adult)")
)

sp <- addEcopathParams(sp, basic_estimates, species_to_groups)

## Create allometric model ----
p <- newAllometricParams(sp)

sp <- p@species_params
no_sp <- nrow(sp) # Number of species

## Set up gear params ----
ecopath_catch <- read_csv("inst/extdata/North Sea-Catch.csv")

p <- addEcopathCatchTotal(p, ecopath_catch)

# Get new steady state with desired biomass and growth
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

## Match catch ----
p <- matchCatch(p, catch = catch)
# Fix the species with unrealistic reproductive efficiency
p <- tuneEcopath(p, catch = catch)

## Match consumption ----
p <- matchConsumption(p)

## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv("inst/extdata/North Sea-Diet composition.csv")
dm <- reduceEcopathDiet(sp, ecopath_diet)
