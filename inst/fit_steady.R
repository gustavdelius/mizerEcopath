# In this script we will use `optim()` to fit the parameters of a mizer model to
# given ecopath parameters and given size-distribution data.

library(mizerEcopath)
library(here)

# Extract a single-species model from an existing model
p <- readRDS(here("../CelticSea/ecopath/Lauria/lauria_HOM_2.rds"))
species <- species_params(p)$species
remove <- species[species != "Cod"]
p <- removeSpecies(p, remove)
sp <- species_params(p)
gp <- gear_params(p)

# We will attempt to estimate the parameters `ks`, `mu`, `catchability`, `L50`,
# `L25`, `w_repro_max`, and `m` for each species in the model.
#
# In this first iteration we will take `w_mat` as known and fixed. In future
# we might estimate it as well while penalising deviations from the specified
# value.

# Because we can later rescale area to get the specified total Biomass and we
# can rescale time to get the specified Consumption, we can keep two parameters
# fixed during optimisation. We choose to fix the encounter rate and the
# number density at w_min

# The objective function will determine the steady state with given parameters
# and calculate the distance from the following data:
# P/Q, C/Q, Catch-size-distribution
# It will also penalise deviations from a reproductive efficiency of 0.01.
objective_fn <- function(pars) {
    sp$ks <- pars[1]
    sp$mu <- pars[2]
    gp$catchability <- pars[3]
    gp$L50 <- pars[4]
    gp$L25 <- pars[5] * L50
    sp$w_repro_max <- pars[6]
    # sp$m <- pars[7]
    species_params(p) <- sp
    gear_params(p) <- gp
    p <- steadySingleSpecies(p)

    Q <- getConsumption(p)
    P <- getProduction(p)
    C <- getYield(p)
    getY
}
