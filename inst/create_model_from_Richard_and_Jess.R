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
# We increase the production of some species to increase
# their mortality which was necessary so that after imposing
# the diet matrix the predation mortality does not exceed
# the external mortality.
species_params(p)["Hake", "production_observed"] <- 0.8
species_params(p)["Herring", "production_observed"] <- 0.8
species_params(p)["Cod", "production_observed"] <- 0.1

# Match ----

pm <- matchCatch(p, catch = catch)

pm <- setFeedingLevels(pm, f = 0.6, f_c = 0.2)
pm <- steady(pm)

# pt <- tuneEcopath(pm, catch = catch, diet = reduced_dm, match = "catch")

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
resource_level(psr) <- 0.1
psr <- steady(psr, tol = 1e-12, t_max = 200)
# psrs <- steadyNewton(psr, reproduction = "dynamic",
#                      verbose = TRUE, stability = TRUE)

# Plots
# params <- psr
# plotYieldVsSpecies(params)
# plotBiomassVsSpecies(params)
# plotlySpectra(params, power = 2, resource = FALSE)
# plotlyFMort(params)
# plotGrowthCurves(params, species_panel = TRUE)
# plotDiet(params, species = "Megrim")

# Reproduction ----

# Tune all species. Changing the reproduction level of one species shifts the
# yield curves of the others a little, so we sweep twice; in practice the
# second sweep returns the same values.
params <- psr

Fcur <- setNames(gp$catchability[gp$gear == "commercial"],
                 species_params(params)$species)
res <- list()
for (sweep in 1:2) {
    print(paste("Sweep ", sweep))
    for (sp in names(Fcur)) {
        r <- tune(params, sp, Fcur[[sp]])
        params <- set_rl(params, sp, r$rl)
        res[[sp]] <- r
        print(paste("Reproduction level for", sp, ":", r$rl))
    }
}
do.call(rbind, lapply(res, as.data.frame))

# Result of the sweeps:
#
#                    rl    side          status
#   Herring       0.010  +0.347  at lower bound
#   Cod           0.990  -0.046  at upper bound
#   Megrim        0.582  -0.009  converged
#   Haddock       0.638  +0.005  converged
#   Whiting       0.010  +0.308  at lower bound
#   Hake          0.990  -0.187  at upper bound
#   Blue whiting  0.010  +0.356  at lower bound
#   Plaice        0.638  +0.014  converged
#   Sole          0.901  -0.007  converged
#
# For Herring, Whiting and Blue whiting the peak stays above F_cur however low
# we push the reproduction level, and the peak barely moves while we do it
# Herring and Blue whiting are
# effectively unfished, so F_MSY is bound to lie far above F_cur. Pushing them
# to the bound would make their reproduction almost density independent for no
# gain, so we leave them at the neutral 0.5.
params <- setBevertonHolt(params, reproduction_level =
                              c(Herring = 0.5, Whiting = 0.5,
                                `Blue whiting` = 0.5))

# Check where the peaks end up
t(vapply(names(Fcur), function(sp) peak(params, sp, Fcur[[sp]]), numeric(2)))

#                 F_peak  F_peak / F_cur
#   Herring       0.0266  >= 8
#   Cod           0.5477     0.90
#   Megrim        0.0873     1.04
#   Haddock       0.1561     0.98
#   Whiting       0.3450     1.72
#   Hake          0.2063     0.56
#   Blue whiting  0.0958  >= 8
#   Plaice        0.1558     1.03
#   Sole          0.1806     0.98
#
# Hake is the one properly fished species we cannot fix: even with the
# reproduction level at 0.99, i.e. with recruitment as insensitive to spawning
# stock as Beverton-Holt allows, its yield peaks at 0.56 * F_cur. The model
# therefore says hake is currently fished well above F_MSY. If that is not
# wanted it has to be addressed elsewhere than in the reproduction level.


# Save ----
saveParams(params, "inst/params_Richard_and_Jess.rds")
