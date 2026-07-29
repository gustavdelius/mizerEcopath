# Create a model that uses Richard's biomasses, yields and catch size
# distributions together with Jess's survey size distributions
library(dplyr)
library(ggplot2)
library(mizer)
library(mizerEcopath)

load("inst/data_Richard_and_Jess.rda")

# Build MizerParams ----

#p <- newVonBertalanffyParams(sp)
p <- newAllometricParams(sp)
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

# Feeding levels. Consumption and metabolic loss depend on f and f_c only
# through their ratio, so keeping f_c = f/3 leaves Q/B and respiration (and the
# whole steady state) unchanged whatever we do to the level of f. What the
# level does change is how strongly growth responds to a change in food:
# dlog(E_r)/dlog(E) = (1 - f) / (1 - f_c/f). Whiting is raised to 0.8406 to damp
# that response, which is what brings its yield curve to peak at its FMSY; see
# the "What the feeding level buys" section of vignettes/Richard_and_Jess.qmd.
# NB setFeedingLevels() assigns f positionally - the names below are for our
# benefit only and the order must match species_params(pm).
f <- setNames(rep(0.6, nrow(sp)), sp$species)
f[["Whiting"]] <- 0.8406

pm <- setFeedingLevels(pm, f = f, f_c = f / 3)

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
species_params(pm)["Herring", "production_observed"] <- 0.8
pm <- matchCatch(pm, catch = catch, species = "Herring")
species_params(pm)["Cod", "production_observed"] <- 0.1
pm <- matchCatch(pm, catch = catch, species = "Cod")
species_params(pm)["Hake", "production_observed"] <- 0.8
pm <- matchCatch(pm, catch = catch, species = "Hake")

# Interactions ----

# # Attempt to use old interaction matrix
# inter <- interaction_matrix(celtic_params)
# species <- species_params(p)$species
# inter <- inter[species, species]
# pi <- makeInteracting(p, interaction = inter)

pd <- matchDiet(pm, reduced_dm)
ps <- steady(pd)

psr <- alignResource(ps)
resource_params(psr)$w_pp_cutoff <- 1
initialNResource(psr)[w_full(psr) > 1] <- 0
comment(psr@cc_pp) <- NULL
psr <- setResourceInteraction(psr,
    resource_dynamics = "resource_semichemostat",
    tol = 1e-2)
psr <- steady(psr)
resource_level(psr) <- 0.5
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

psr <- set_cannib(psr, "Cod", 0.1875)

# Reproduction ----

# Tune all species. Changing the reproduction level of one species shifts the
# yield curves of the others, so we sweep twice. The second sweep does not
# simply reproduce the first here - hake moves from 0.36 to 0.14 - so treat two
# passes as a truncated iteration rather than a converged answer.
source("inst/tune_repro_level.R")
params <- setBevertonHolt(psr, reproduction_level = 0.5)

F_target <- setNames(species_params(params)$FMSY,
                 species_params(params)$species)
missing <- is.na(F_target)
gear_sel <- gear_params(params)$gear == "commercial"
F_cur <- gear_params(params)$catchability[gear_sel]
F_target[missing] <- F_cur[missing]

res <- list()
# The loop variable is `s`, not `sp`: `sp` is the species data frame loaded at
# the top and the feeding-level block above needs it to survive.
for (sweep in 1:2) {
    print(paste("Sweep ", sweep))
    for (s in names(F_target)) {
        r <- tune(params, s, F_target[[s]])
        params <- set_rl(params, s, r$rl)
        res[[s]] <- r
        print(paste("Reproduction level for", s, ":", r$rl))
    }
}
do.call(rbind, lapply(res, as.data.frame))

# Check where the peaks end up
t(vapply(names(F_target), function(sp) peak(params, sp, F_target[[sp]]), numeric(2)))

# Save ----
saveParams(params, "inst/params_final_Richard_and_Jess.rds")
