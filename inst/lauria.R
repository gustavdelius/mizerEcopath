# The task is to create a mizer model with a steady state that agrees with
# the Ecopath model for the Celtic Sea by Hernvann et.al.

# We start by creating the species with power-law consumption and mortality,
# then setting the predation kernel using stomach data and
# and the interaction matrix so as to reproduce the ecopath diet matrix.

library(mizerEcopath)
library(mizerExperimental)
library(dplyr)

## Set up species parameters ----

sp <- read.csv("../CelticSea/ecopath/Lauria/species_params.csv")
# keep only selected parameters
sp <- sp |> select(species, w_max, w_mat, a, b, LWRSource, age_mat)

no_sp <- nrow(sp) # Number of species
sp["n"] <- 0.7 # Exponent of consumption
sp["p"] <- 0.7 # Exponent of respiration (metabolism)
sp["d"] <- -0.3 # Exponent of mortality (in Mizer n - 1)
sp["alpha"] <- 0.8  # Ecopath default (conversion efficiency)
sp$gonad_proportion <- 0.2 # Proportion of total production which goes into gonads

# Add predation kernel params from stomach data
fits <- readRDS("../CelticSea/data/stomach_data_fit.rds")
sp$pred_kernel_type <- "power_law"
lambda <- 2
sp$kernel_exp <- fits$alpha + 4/3 - lambda
sp$kernel_l_l <- fits$ll
sp$kernel_u_l <- fits$ul
sp$kernel_l_r <- fits$lr
sp$kernel_u_r <- fits$ur

# Dictionary between species and ecopath groups
species_to_groups <- list(
    "Herring" = "Herring",
    "Cod" = c("Cod", "Juvenile cod"),
    "Haddock" = c("Haddock", "Juvenile haddock"),
    "Whiting" = c("Whiting", "Juvenile whiting"),
    "Blue whiting" = c("Blue whiting", "Juvenile blue whiting"),
    "Hake" = c("Hake", "Juvenile hake"),
    "Monkfish" = c("Monkfish", "Juvenile monkfish"),
    "Horse mackerel" = "Horse mackerel",
    "Mackerel" = "Mackerel",
    "Plaice" = c("Plaice", "Juvenile plaice"),
    "Megrim" = c("Megrim", "Juvenile megrim"),
    "Sole" = "Sole"
)

## Add Ecopath species parameters ----
ecopath_params <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Basic estimates.csv")
ecopath_params <- validEcopathParams(ecopath_params, species_to_groups)
sp <- addEcopathParams(sp, ecopath_params, species_to_groups)

## Create MizerParams object ----
p <- newAllometricParams(sp, no_w = 200, lambda = lambda)
sp <- p@species_params

## Set up gear params ----
# For the lauria model we do not yet have the catch data
# We'll calculate it from the fishing mortalities and biomases
catch <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Fishing mortality rates.csv")
catch$totalF <- rowSums(catch[, -(1:2)], na.rm = TRUE)
catch <- left_join(ecopath_params, catch, by = c("Group.name" = "Fleet.group")) |>
    mutate(TotalCatch..t.km..year. = totalF * Biomass..t.km..) |>
    select(Group.name, TotalCatch..t.km..year.)

p <- addEcopathCatchTotal(p, catch)

# Get new steady state
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

p_backup <- p

## Match consumption and yield ----
p <- p_backup |>
    matchConsumption() |> matchYield(keep = "biomass") |>
    matchConsumption() |> matchYield(keep = "biomass") |>
    matchConsumption() |> matchYield(keep = "biomass")

## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Diet composition.csv")
dm <- reduceEcopathDiet(sp, ecopath_diet)

# Prepare catch data
catch <- readRDS("../CelticSea/data/ecopath_catch_v2.rds")
catch$dl <- 1
catch$gear <- as.character(catch$gear)
catch$gear[catch$gear == "commercial"] <- "total"

# Change case for Horse Mackerel:
catch <- catch %>%
    mutate(species = ifelse(species == "Horse Mackerel", "Horse mackerel", species))

## Tune ----
p <- tuneEcopath(p, catch = catch, diet = dm,
                 tabs = c("Spectra", "Catch", "Growth", "Repro", "Diet", "Death"))

## Switch on interactions
p <- matchDiet(p, dm)

## Switch on satiation
p <- setFeedingLevel(p, 0.6)

# Check that steady state has not changed
ps <- p |> steadySingleSpecies()
waldo::compare(initialN(p), initialN(ps), tolerance = 1e-9)
# Check that ecopath is still matched
isEcopathMatched(p)
# Check that diet matrix is matched
Kn <- getDietMatrix(p)[, 1:no_sp]
# Convert ecopath diet matrix from proportions to absolute consumption
Q <- sp$ecopath_consumption
dm <- dm * Q / rowSums(dm)
# Drop the "other" column
D <- dm[, -ncol(dm)]
all.equal(Kn, D, tolerance = 1e-2)


# Tune dynamics ----

# We now have an interacting model with the correct steady state size spectra.
# From now on we only want to tune the parameters that affect the dynamics
# away from steady state. These are
# - reproduction levels
# - feeding levels
# - resource level
