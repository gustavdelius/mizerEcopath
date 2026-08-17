library(dplyr)
library(mizer)
library(mizerEcopath)

load("inst/data_Richard_and_Jess.rda")

p <- newVonBertalanffyParams(sp)
gear_params(p) <- gp
initial_effort(p) <- 1
p <- steadySingleSpecies(p) |> setBevertonHolt()
p <- matchBiomasses(p)
species_params(p)$production_observed <- getSomaticProduction(p)
pm <- matchCatch(p, catch = catch)

species_params(pm)["Cod", "production_observed"] <- 0.2
pm <- matchCatch(pm, catch = catch, species = "Cod")
species_params(pm)["Hake", "production_observed"] <- 1.0
pm <- matchCatch(pm, catch = catch, species = "Hake")
species_params(pm)["Blue whiting", "production_observed"] <- 4.0
pm <- matchCatch(pm, catch = catch, species = "Blue whiting")

pd <- matchDiet(pm, reduced_dm)
ps <- steady(pd)
psr <- alignResource(ps)
comment(psr@cc_pp) <- NULL
psr <- setResourceInteraction(psr, resource_dynamics = "resource_semichemostat", tol = 1e-2)
psr <- steady(psr)
resource_level(psr) <- 0.99
psr <- steady(psr, tol = 1e-12, t_max = 200)

set_cannib <- function(params, species, lambda) {
    mort_old <- getPredMort(params)
    enc_old  <- getEncounter(params)
    inter <- interaction_matrix(params)
    inter[species, species] <- inter[species, species] * lambda
    interaction_matrix(params) <- inter
    ext_encounter(params) <- ext_encounter(params) + (enc_old  - getEncounter(params))
    ext_mort(params)      <- ext_mort(params)      + (mort_old - getPredMort(params))
    params
}

# we'll do cannibalism lambda for Cod. Let's try 0 (no cannibalism) or some other value.
# But wait, with vonB the cannibalism might not be the same. 
# Let's test with no cannibalism for now.
psr <- set_cannib(psr, "Cod", 2.0)
psr <- set_cannib(psr, "Whiting", 0.0)
psr <- set_cannib(psr, "Hake", 2.0)

source("inst/tune_repro_level.R")

params <- setBevertonHolt(psr, reproduction_level = 0.5)

F_target <- setNames(species_params(params)$FMSY, species_params(params)$species)
missing <- is.na(F_target)
gear_sel <- gear_params(params)$gear == "commercial"
F_cur <- gear_params(params)$catchability[gear_sel]
F_target[missing] <- F_cur[missing]

res <- list()
for (sweep in 1:2) {
    for (s in names(F_target)) {
        r <- tune(params, s, F_target[[s]])
        params <- set_rl(params, s, r$rl)
        res[[s]] <- r
    }
}
res_df <- do.call(rbind, lapply(res, as.data.frame))
saveRDS(list(params=params, res=res_df), "inst/tune_results_vonB.rds")
print(res_df)
