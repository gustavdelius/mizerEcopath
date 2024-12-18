# The task is to create a mizer model with a steady state that agrees with
# the Ecopath model for the Celtic Sea by Hernvann et.al.

# We start by creating the species with power-law consumption and mortality,
# then setting the predation kernel using stomach data and
# and the interaction matrix so as to reproduce the ecopath diet matrix.

library(mizerEcopath)
library(mizerExperimental)

## Set up species parameters ----

sp <- read.csv(system.file("extdata/species_params_celtic_sea.csv",
                                       package = "mizerEcopath"))

no_sp <- nrow(sp) # Number of species
sp["n"] <- 0.7 # Exponent of encounter
sp["p"] <- 0.7 # Exponent of respiration (metabolism)
sp["d"] <- -0.3 # Exponent of mortality (in Mizer n - 1)
sp["M"] <- 1 # Mortality rate at 1g
sp["alpha"] <- 0.8  # Ecopath default (conversion efficiency)
sp$gonad_proportion <- 0.2 # Proportion of total production which goes into gonads

# Add predation kernel params from stomach data ----
fits <- readRDS("../CelticSea/data/stomach_data_fit.rds")
# This dataframe uses latin names for the species, so we add these to the
# species_params dataframe.
latin_names <- c(
    "Herring" = "Clupea harengus",
    "Sprat" = "Gadus morhua",
    "Cod" = "Melanogrammus aeglefinus",
    "Haddock" = "Merlangius merlangus",
    "Whiting" = "Micromesistius poutassou",
    "Blue whiting" = "Merluccius merluccius",
    "Norway Pout" = "Lophius piscatorius",
    "European Hake" = "Trachurus trachurus",
    "Horse Mackerel" = "Scomber scombrus",
    "Mackerel" = "Pleuronectes platessa",
    "Plaice" = "Lepidorhombus whiffiagonis",
    "Megrim" = "Solea solea",
    "Sole" = "Solea solea"
)
sp$latin_names <- latin_names[sp$species]
row.names(fits) <- fits$species
sp$pred_kernel_type <- "power_law"
lambda <- 2
sp$kernel_exp <- fits[sp$latin_names, "alpha"] + 4/3 - lambda
sp$kernel_l_l <- fits[sp$latin_names, "ll"]
sp$kernel_u_l <- fits[sp$latin_names, "ul"]
sp$kernel_l_r <- fits[sp$latin_names, "lr"]
sp$kernel_u_r <- fits[sp$latin_names, "ur"]

## Add Ecopath species parameters ----
ecopath_params <- read.csv(
    system.file("extdata/celtic_sea_hernvann_et_al/basic_estimates.csv",
                package = "mizerEcopath"))
# Dictionary between species and ecopath groups
species_to_groups <- list(
    "Herring" = "Herring",
    "Sprat" = "Sprat",
    "Cod" = c("Cod large", "Cod small"),
    "Haddock" = "Haddock",
    "Whiting" = "Whiting",
    "Blue whiting" = "Blue whiting",
    "Norway Pout" = "Pouts",
    "European Hake" = c("Hake large", "Hake small"),
    "Horse Mackerel" = "Horse mackerel",
    "Mackerel" = "Mackerel",
    "Plaice" = "Plaice",
    "Megrim" = "Megrim",
    "Sole" = "Sole"
)
sp <- addEcopathParams(sp, ecopath_params, species_to_groups)

## Create MizerParams object ----
p <- newAllometricParams(sp, no_w = 200, lambda = lambda)
sp <- p@species_params

## Set up gear params ----
catch <- read.csv(system.file("extdata/celtic_sea_hernvann_et_al/catch.csv",
                              package = "mizerEcopath"))
p <- addEcopathCatchTotal(p, catch)

# Get new steady state
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

plotlySpectra(p)
p_backup <- p

## Match to Ecopath params ----
p <- p_backup |> matchEcopath()

# Set resource to be in line with fish
total <- colSums(initialN(p))
fish_sel <- p@w_full >= p@w[1]
ratio <- max(total / initialNResource(p)[fish_sel])
p <- scaleDownBackground(p, 1/ratio)
plotSpectra(p)

# saveParams(p, here("ecopath/Hernvann/Celtic_13_ecopath.rds"))

## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv(
    system.file("extdata/celtic_sea_hernvann_et_al/diet_composition.csv",
                package = "mizerEcopath"))
dm <- reduceEcopathDiet(sp, ecopath_diet)

# Prepare catch data
catch <- readRDS("../CelticSea/data/ecopath_catch_v1.rds")
catch$dl <- 1
catch$gear <- as.character(catch$gear)
catch$gear[catch$gear == "commercial"] <- "total"

## Tune ----
p <- tuneEcopath(p, catch = catch, diet = dm,
                 tabs = c("Spectra", "Catch", "Growth", "Repro", "Diet", "Death"))

## Switch on interactions
p <- matchDiet(p, dm)

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
