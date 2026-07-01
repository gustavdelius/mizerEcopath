library(mizerEcopath)
library(dplyr)
library(mizer)
newVonBertalanffyParams <- function(species_params, no_w = 200, max_w = NULL) {
    sp <- validGivenSpeciesParams(species_params)

    # Impose relation between exponents
    sp <- set_species_param_default(sp, "n", 1 - 1 / sp$b)
    sp$p <- sp$n
    sp <- set_species_param_default(sp, "d", sp$n - 1)

    # Set default assimilation efficiency
    sp <- set_species_param_default(sp, "alpha", 0.8)

    # Switch off metabolic respiration
    sp$ks <- 0
    # Switch off constant mortality
    sp$z0 <- 0

    # Generate a default mizer model with the desired species We extend the
    # resource spectrum over the entire size range to ensure that all species
    # have sufficient prey throughout their life.

    if (is.null(max_w)) {max_w <- max(sp$w_max)}

    if (max_w < max(sp$w_max)) {
        warning("The maximum weight provided (max_w) is lower than the
                maximum size of the fish. The model has been generated with
                the latter")
    }

    max_w <- max(max_w, sp$w_max)

    p <- newMultispeciesParams(sp, no_w = no_w, info_level = 0,
                               # extend resource over entire size range
                               max_w = max_w,
                               w_pp_cutoff = max_w * (1 + 1e-9),
                               resource_dynamics = "resource_constant",
                               n = 0.7)
    sp <- p@species_params

    # Switch off all interactions
    interaction_matrix(p)[] <- 0
    sp$interaction_resource <- 0

    # Switch off satiation
    sp$h <- Inf
    intake_max(p)[] <- Inf

    # Set constant metabolic rate
    sp$k <- sp$b * sp$k_vb

    # Set power-law encounter rate (the coefficient will be adjusted below)
    sp$E_ext <- sp$b * sp$k_vb * sp$w_inf^(1/sp$b) / sp$alpha

    species_params(p) <- sp

    # Subtract reproduction investment from metabolic loss rate
    metab(p) <- metab(p) - getERepro(p)

    # Set power-law mortality
    # Choose a positive coefficient so that the juvenile biomass density
    # has a slightly negative slope (of -0.2).
    sp$z_ext <- sp$alpha * sp$E_ext * (1 + 0.2 - sp$n)

    species_params(p) <- sp
    # Match Biomasses
    p <- matchBiomasses(p)
    # Set to steady state
    p <- steadySingleSpecies(p, keep = "biomass")
    p <- setBevertonHolt(p, reproduction_level = 0)

    return(p)
}

Jess_sp <- celtic_params@species_params |>
    select(species, a, b, age_mat, w_mat, w_max,
           biomass_cutoff, biomass_observed,
           pred_kernel_type, beta, sigma,
           kernel_exp, kernel_l_l, kernel_u_l, kernel_l_r, kernel_u_r)
James_sp <- readRDS("~/Git/mizerEcopath/inst/James_sp.rds")
sp <- James_sp |>
    select(species, w_inf, k_vb, t0,
           production_observed, consumption_observed) |>
    inner_join(Jess_sp, by = "species")

sp$D_ext <- 1
p <- newVonBertalanffyParams(sp)

gp <- gear_params(celtic_params) |>
    filter(species %in% sp$species)
gear_params(p) <- gp
initial_effort(p) <- 1
p <- steadySingleSpecies(p) |>
    setBevertonHolt()

p <- matchCatch(p, catch = celtic_catch, production_lambda = 0,
                yield_lambda = 100)

p <- tuneEcopath(p, catch = celtic_catch, diet = reduced_dm)

plotSpectra(p, power = 2, resource = FALSE)
plotBiomassVsSpecies(p)
plotYieldVsSpecies(p)
plot_catch(p, catch = celtic_catch, species = "Cod")
plotProductionVsSpecies(p)
p@species_params$z_ext
p@species_params$D_ext
sim <- project(p, t_max = 10)
plotBiomass(sim)
