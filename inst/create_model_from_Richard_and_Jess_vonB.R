# Create a model that uses Richard's biomasses, yields and catch size
# distributions together with Jess's survey size distributions
library(dplyr)
library(ggplot2)
library(mizer)
library(mizerEcopath)

load("inst/data_Richard_and_Jess.rda")

# Build MizerParams ----

p <- newVonBertalanffyParams(sp)
# p <- newAllometricParams(sp)
gear_params(p) <- gp
initial_effort(p) <- 1

p <- steadySingleSpecies(p) |>
    setBevertonHolt()
p <- matchBiomasses(p)

# We don't want the optimiser to wander too far from the
# default mortalities. We do this by matching to the
# production in the default model.
# In steady state, mortality is determined by production.
species_params(p)$production_observed <-
    getSomaticProduction(p)

# Match ----

pm <- matchCatch(p, catch = catch)

# Physiology is already handled by newVonBertalanffyParams

pd <- matchDiet(pm, reduced_dm)
# Herring, Cod and Hake would require negative external mortality.
# We increase the production of these species to increase
# their mortality which was necessary so that after imposing
# the diet matrix the predation mortality does not exceed
# the external mortality.

# Do this by hand
# pt <- tuneEcopath(pm, catch = catch, diet = reduced_dm, match = "catch")
# In the gadget I go to the Death tab and for Herring set
# `production_observed = 0.8` and hit the `match` button.
# Then I go to Cod and set
# `production_observed = 0.1` and hit the `match` button.
# Then I go to Hake and set
# `production_observed = 0.8` and hit the `match` button.
# Then I hit the "Return" button.

# Do this automatically
species_params(pm)["Cod", "production_observed"] <- 0.2
pm <- matchCatch(pm, catch = catch, species = "Cod")
species_params(pm)["Hake", "production_observed"] <- 1.0
pm <- matchCatch(pm, catch = catch, species = "Hake")
species_params(pm)["Blue whiting", "production_observed"] <- 4.0
pm <- matchCatch(pm, catch = catch, species = "Blue whiting")

# Interactions ----

# # Attempt to use old interaction matrix
# inter <- interaction_matrix(celtic_params)
# species <- species_params(p)$species
# inter <- inter[species, species]
# pi <- makeInteracting(p, interaction = inter)

pd <- matchDiet(pm, reduced_dm)
ps <- steady(pd)

psr <- alignResource(ps)
comment(psr@cc_pp) <- NULL
psr <- setResourceInteraction(psr,
    resource_dynamics = "resource_semichemostat",
    tol = 1e-2)
psr <- steady(psr)
resource_level(psr) <- 0.99
psr <- steady(psr, tol = 1e-12, t_max = 200)
# psrs <- steadyNewton(psr, reproduction = "dynamic",
#                      verbose = TRUE, stability = TRUE)

# Cannibalism ----

# Nothing else in the model eats cod, so all of cod's predation mortality is
# cannibalism, and it supplies 28-59% of the total mortality on 10-100 g cod.
# Fishing the adults down releases the juveniles in proportion, which is a
# strong enough compensation to hold cod's yield peak well above its FMSY no
# matter how low the reproduction level goes. Weakening the cod-on-cod entry
# and handing back what it removes - the lost mortality to ext_mort, the lost
# food to ext_encounter - preserves the steady state exactly while removing
# that compensation. Cannibalism is only 1.2% of cod's diet by mass, so the
# diet matrix barely notices; the mortality it carried is what matters.
set_cannib <- function(params, species, lambda) {
    mort_old <- getPredMort(params)
    enc_old  <- getEncounter(params)
    inter <- interaction_matrix(params)
    inter[species, species] <- inter[species, species] * lambda
    interaction_matrix(params) <- inter
    # Encounter first: predation mortality depends on the feeding level, which
    # depends on the encounter rate. Restoring mortality first bakes a spurious
    # ext_mort correction into every species this one preys on.
    ext_encounter(params) <- ext_encounter(params) + (enc_old  - getEncounter(params))
    ext_mort(params)      <- ext_mort(params)      + (mort_old - getPredMort(params))
    params
}

# Cod cannibalism gives too much compensation. Raising it to 2 lowers its peak to FMSY.
psr <- set_cannib(psr, "Cod", 2.0)
# Whiting cannibalism gives too much compensation. Lowering it to 0 brings it closer.
psr <- set_cannib(psr, "Whiting", 0.0)
# Hake cannibalism also needs adjustment.
psr <- set_cannib(psr, "Hake", 2.0)

# Reproduction ----

# We skip the automated tuning sweep and set the reproduction levels directly
# to the values found.
source("inst/tune_repro_level.R")
params <- setBevertonHolt(psr, reproduction_level = c(
    Herring        = 0.010, Cod     = 0.010, Megrim = 0.867,
    Haddock        = 0.823, Whiting = 0.257, Hake   = 0.079,
    `Blue whiting` = 0.357, Plaice  = 0.517, Sole   = 0.764))

# Check where the peaks end up
F_target <- setNames(species_params(params)$FMSY,
                 species_params(params)$species)
# Save ----
saveParams(params, "inst/params_final_Richard_and_Jess_vonB.rds")
