#' Set up parameters for a model with allometric encounter and mortality rates
#'
#' This function sets up a model with allometric encounter and mortality rates
#' for each species and no interaction between the species or with the resource.
#' The initial abundances are set to the steady-state abundances with the total
#' biomass for each species matched to observations.
#' The metabolic respiration rate, the feeding level and the reproduction level
#' are set to zero.
#'
#' If the exponent `n` of the power-law encounter rate is not provided as a
#' species parameter, it is set to 0.7.  If the exponent `d` of the power-law
#' mortality rate is  not provided as a species parameter, it is set to `n – 1`.
#'
#' If the species parameter `alpha`, which gives the proportion of the
#' consumption that is assimilated, is not provided, it is set to 0.8, the
#' default value used by Ecopath.
#'
#' The encounter rate coefficient for each species is chosen so that the
#' species reaches its maturity weight `w_mat` by the age at maturity
#' `age_mat`.  The mortality-rate coefficient is then derived from that growth
#' rate and set so that the juvenile biomass spectrum has a slope of –0.2.
#'
#' The function uses `matchBiomasses()` to match the biomass to the observations
#' and `mizerEcopath::steadySingleSpecies()` to bring each species to steady state. It then
#' calls `setBevertonHolt()` to set the reproduction level to zero
#'
#' Because the model does not make use of the resource spectrum, the resource
#' dynamics is switched off.
#'
#' @param species_params A data frame with species parameters
#' @param no_w The number of weight bins to use in the model
#' @return A MizerParams object
#' @export
newAllometricParams <- function(species_params, no_w = 200) {
    sp <- validGivenSpeciesParams(species_params)

    # Impose relation between exponents
    sp <- set_species_param_default(sp, "n", 0.7)
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
    max_w <- max(sp$w_max)
    p <- newMultispeciesParams(sp, no_w = no_w, info_level = 0,
                               # extend resource over entire size range
                               max_w = max_w,
                               w_pp_cutoff = max_w * (1 + 1e-9),
                               resource_dynamics = "resource_constant")
    sp <- p@species_params

    # Switch off all interactions
    interaction_matrix(p)[] <- 0
    sp$interaction_resource <- 0

    # Switch off satiation
    sp$h <- Inf
    intake_max(p)[] <- Inf

    species_params(p) <- sp

    # Set power-law encounter rate (the coefficient will be adjusted below)
    ext_encounter(p) <- t(outer(p@w, sp$n, "^"))
    # Set encounter rate coefficient to produce desired growth rate
    sp <- set_species_param_default(sp, "age_mat", age_mat_vB(sp))
    factor <- age_mat(p) / sp$age_mat
    ext_encounter(p) <- sweep(ext_encounter(p), 1, factor, "*")
    # Determine power-law coefficient for encounter rate
    e0 <- getEGrowth(p)[, 1] / w(p)[1] ^ sp$n[1]

    # Set power-law mortality
    # Choose a positive coefficient so that the juvenile biomass density
    # has a slightly negative slope (of -0.2).
    mu0 <- e0 * (1 + 0.2 - sp$n)
    ext_mort(p) <- sweep(t(outer(p@w, sp$d, "^")), 1, mu0, "*")

    # Match Biomasses
    p <- matchBiomasses(p)
    # Set to steady state
    p <- mizerEcopath::steadySingleSpecies(p, keep = "biomass")
    p <- setBevertonHolt(p, reproduction_level = 0)

    return(p)
}

#' Test whether a model has allometric encounter and mortality rates
#'
#' This function returns TRUE when the model has allometric encounter and
#' mortality rates and FALSE otherwise.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance.
#' @return TRUE if the model has allometric encounter and mortality rates
#' @export
isAllometric <- function(params, tol = 1e-5) {
    sp <- params@species_params
    sp <- set_species_param_default(sp, "d", sp$n - 1)
    m <- getMort(params) - getFMort(params)
    e <- getEncounter(params)
    for (species in sp$species) {
        spc <- sp[species, ]
        mc <- m[species, ] / params@w ^ spc$d
        if (anyNA(mc) || any(is.nan(mc)) || any(mc < 0) || all(mc == 0)) {
            stop("The model has invalid mortality rates.")
        }
        if (abs((max(mc) - min(mc)) / min(mc)) > tol) {
            return(FALSE)
        }
        ec <- e[species, ] / params@w ^ spc$n
        if (anyNA(ec) || any(is.nan(ec)) || any(ec < 0) || all(ec == 0)) {
            stop("The model has invalid encounter rates.")
        }
        if (abs((max(ec) - min(ec)) / min(ec)) > tol) {
            return(FALSE)
        }
    }
    return(TRUE)
}
