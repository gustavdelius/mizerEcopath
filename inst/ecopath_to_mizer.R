# The task is to create a mizer model with a steady state that agrees with
# the Ecopath model for the Celtic Sea by Hernvann et.al.

# We start by creating the species with power-law consumption and mortality,
# then setting the predation kernel using stomach data and
# and the interaction matrix so as to reproduce the ecopath diet matrix.

library(mizerEcopath)
library(mizerExperimental)

sp <- read.csv(system.file("extdata/species_params_celtic_sea.csv",
                                       package = "mizerEcopath"))

no_sp <- nrow(sp) # Number of species
sp["n"] <- 0.7 # Exponent of consumption
sp["p"] <- 0.7 # Exponent of respiration (metabolism)
sp["d"] <- -0.3 # Exponent of mortality (in Mizer n - 1)
sp["alpha"] <- 0.8  # Ecopath default (conversion efficiency)
sp$gonad_proportion <- 0.2 # Proportion of total production which goes into gonads

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

## Add Ecopath species parameters ----
ecopath_params <- read.csv(
    system.file("extdata/celtic_sea_hernvann_et_al/basic_estimates.csv",
                package = "mizerEcopath"))
sp <- addEcopathParams(sp, ecopath_params, species_to_groups)

## Set up gear params ----
catch <- read.csv(system.file("extdata/celtic_sea_hernvann_et_al/catch.csv",
                              package = "mizerEcopath"))
gp <- data.frame(
    species = sp$species,
    gear = "total",
    sel_func = "sigmoid_length",
    l50 = w2l(sp$w_mat, sp),
    l25 = w2l(sp$w_mat, sp) * 0.8,
    catchability = 1,
    yield_observed = 0
)
for (i in seq_len(nrow(sp))) {
    for (group in sp$ecopath_groups[[i]]) {
        yield <- catch$TotalCatch..t.km..year.[catch$Group.name == group]
        gp$yield_observed[i] <- gp$yield_observed[i] + yield
    }
}
gp$catchability <- gp$yield_observed / sp$biomass_observed

## Create model ----
p <- newMultispeciesParams(sp, gear_params = gp, initial_effort = 1,
                           no_w = 200, lambda = 2, info_level = 0)
sp <- p@species_params

## Make non-interacting ----
# Extend the resource to maximum size and switch off dynamics
p <- setResource(p, w_pp_cutoff = max(p@w) * 0.999,
                 resource_dynamics = "resource_constant")
# Use power law all the way
lambda <- p@resource_params$lambda
kappa <- initialNResource(p)[1] * p@w_full[1]^lambda
initialNResource(p) <- kappa * p@w_full^(-lambda)

# Power-law mortality
comment(p@mu_b) <- "Using power-law mortality"
species_params(p)$z0 <- 0
for (species in sp$species) {
    spc <- sp[species, ]
    # calculate power-law mortality with value 0.4 at maturity size
    power_mort <- 0.4 * (p@w / spc$w_mat)^(spc$n - 1)
    # add it as background mortality
    ext_mort(p)[species, ] <- power_mort
}

# Switch off all species interactions
p@interaction[] <- 0
ext_encounter(p) <- getEncounter(p)
species_params(p)$interaction_resource[] <- 0

# Get new steady state
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

plotlySpectra(p)
p_backup <- p

## Match to Ecopath params ----
p <- p_backup |> match()

# Turn off satiation
species_params(p)$h <- Inf
ext_encounter(p) <- ext_encounter(p) * 0.4

# Set resource to be in line with fish
total <- colSums(initialN(p))
fish_sel <- p@w_full >= p@w[1]
ratio <- max(total / initialNResource(p)[fish_sel])
p <- scaleDownBackground(p, 1/ratio)
plotSpectra(p)

# saveParams(p, here("ecopath/Hernvann/Celtic_13_ecopath.rds"))

getEcotrophicEfficiency(p)
getM0B(p)

## Calibrate size spectra ----
# This step currently needs to be done with tuneEcopath(),
# using the tune.R script


# saveParams(p, here("ecopath/Hernvann/Celtic_13_ecopath_catch.rds"))

# Interacting model ----

# From here onwards we do not want to change steady-state growth and mortality
# rates any more in order to not disturb the steady-state size spectra.

## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv(
    system.file("extdata/celtic_sea_hernvann_et_al/diet_composition.csv",
                package = "mizerEcopath"))
dm <- reduceEcopathDiet(ecopath_diet, species_to_groups)

p <- matchDiet(p, dm)

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

# tuneParams ----
source("../new_tuning/singleControl.R")
p <- tuneParams(p, controls = c("single", "predation", "fishing",
                                "reproduction", "other"),
                match = c("growth", "biomass"))
