# In this script we will use `optim()` to fit the parameters of a mizer model to
# given ecopath parameters and given size-distribution data.

library(mizerExperimental)

# We will attempt to estimate the parameters `ks`, `mu`, `catchability`, `L50`,
# `L25`, `w_repro_max`, and `m` for each species in the model.
#
# In this first iteration we will take `w_mat` as known and fixed. In future
# we might estimate it as well while penalising deviations from the specified
# value.
#
# The objective function will determine the steady state with given parameters
# and calculate the distance from the following data:
# Production, Catch, Catch-size-distribution
# It will also penalise deviations from a given `age_mat` and from a
# reproductive efficiency of 0.01.
objective_fn <- function(pars) {
    sp <- p@species_params[i, ]
    gp <- p@gear_params[i, ]
    sp$ks <- pars[1]
    sp$mu <- pars[2]
    gp$catchability <- pars[3]
    gp$L50 <- pars[4]
    gp$L25 <- pars[5] * L50
    sp$w_repro_max <- pars[6]
    sp$m <- pars[7]
    species_params(p)[i, ] <- sp
    gear_params(p)[i, ] <- gp

    p <- steadySingleSpecies(p, species = i)
    p <- matchBiomasses(p, species = i)
    p <- matchConsumption(p, species = i)

    P <- getProduction(p)
    C <- getYield(p)
}
