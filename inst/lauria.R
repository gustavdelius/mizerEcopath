# The task is to create a mizer model with a steady state that agrees with
# the Ecopath model for the Celtic Sea by Lauria et.al.

# We start by creating the species with power-law consumption and mortality,
# then setting the predation kernel using stomach data and
# and the interaction matrix so as to reproduce the ecopath diet matrix.

library(mizerEcopath)
library(dplyr)

## Prepare catch data ----
# We do this first because we will use it to set maximum sizes for species
catch <- readRDS("../CelticSea/data/ecopath_catch_v2.rds")
catch$dl <- 1
catch$gear <- as.character(catch$gear)
catch$gear[catch$gear == "commercial"] <- "total"
# Change case for Horse Mackerel:
catch <- catch %>%
    mutate(species = ifelse(species == "Horse Mackerel", "Horse mackerel", species))
# Select "total" gear only
catch <- catch[catch$gear == "total", ]

## Set up species parameters ----

sp <- read.csv("../CelticSea/ecopath/Lauria/species_params.csv")

# Set maximum weights to at least 1.2 times the maximum observed weight
sp <- catch |>
    group_by(species) |>
    summarise(max_observed_length = max(length + dl)) |>
    right_join(sp, by = c("species" = "species")) |>
    mutate(max_observed_weight = a * max_observed_length^b) |>
    mutate(w_max = pmax(w_max, 1.2 * max_observed_weight))

# keep only selected parameters
sp <- sp |> select(species, w_max, w_mat, a, b, LWRSource, age_mat)
# Add gonadic proportion
sp$gonad_proportion <- 0.2

# Change age_mat of Mackerel
# sp$age_mat[sp$species == "Mackerel"] <- 2.8

## Add Ecopath species parameters ----

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

ecopath_params <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Basic estimates.csv")
sp <- addEcopathParams(sp, ecopath_params, species_to_groups)

## Create allometric model ----
p <- newAllometricParams(sp)

sp <- p@species_params
no_sp <- nrow(sp) # Number of species

## Set up gear params ----
# For the lauria model we do not yet have the catch data
# Because we do not have access to the Ecopath Catch data frame, we calculate it
# from the fishing mortalities and biomases
fmort <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Fishing mortality rates.csv")
fmort$totalF <- rowSums(fmort[, -(1:2)], na.rm = TRUE)
ecopath_catch <- left_join(ecopath_params, fmort, by = c("Group.name" = "Fleet.group")) |>
    mutate(TotalCatch..t.km..year. = totalF * Biomass..t.km.2.) |>
    select(Group.name, TotalCatch..t.km..year.)

p <- addEcopathCatchTotal(p, ecopath_catch)

# Get new steady state with desired biomass and growth
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

## Match catch ----
p <- matchCatch(p, catch = catch)
# Fix the species with unrealistic reproductive efficiency (Monkfish can't be fixed)
p <- tuneEcopath(p, catch = catch)
p_backup <- p

## Match consumption ----
p <- matchConsumption(p)

## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv("../CelticSea/ecopath/Lauria/Ecopath-Diet composition.csv")
dm <- reduceEcopathDiet(sp, ecopath_diet)

## Add predation kernel params from stomach data ----
fits <- readRDS("../CelticSea/data/stomach_data_fit.rds")
sp <- species_params(p)
sp$pred_kernel_type <- "power_law"
lambda <- 2
sp$kernel_exp <- fits$alpha + 4/3 - lambda
sp$kernel_l_l <- fits$ll
sp$kernel_u_l <- fits$ul
sp$kernel_l_r <- fits$lr
sp$kernel_u_r <- fits$ur
species_params(p) <- sp

## Tune ----
p <- tuneEcopath(p, catch = catch, diet = dm,
                 tabs = c("Spectra", "Catch", "Growth", "Repro", "Diet", "Death"))

## Switch on interactions ----
p <- matchDiet(p, dm)

## Switch on satiation ----
p <- setFeedingLevel(p, 0.6)

# Check that steady state has not changed
ps <- p |> steadySingleSpecies()
waldo::compare(initialN(p), initialN(ps), tolerance = 1e-9)
# Check that ecopath is still matched
isEcopathMatched(p)
# Check that diet matrix is matched
Kn <- getDietMatrix(p)[, 1:no_sp]
# Convert ecopath diet matrix from proportions to absolute consumption
Q <- sp$consumption_observed
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
