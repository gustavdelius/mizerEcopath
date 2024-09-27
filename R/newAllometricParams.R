#' Set up parameters for a model with allometric encounter and mortality rates
#'
#' This function sets up a multispecies model with allometric encounter and
#' mortality rates.
#'
#' @param species_params A data frame with species parameters
#' @param no_w The number of weight bins to use in the model
#' @param lambda The exponent of the Sheldon spectrum
#' @return A MizerParams object
#' @export
newAllometricParams <- function(species_params, no_w = 200, lambda = 2) {
    p <- newMultispeciesParams(
        sp, no_w = no_w, lambda = lambda, info_level = 0)
    sp <- p@species_params

    # Extend the resource to maximum size
    p <- setResource(p, w_pp_cutoff = max(p@w) * 0.999,
                     resource_dynamics = "resource_constant")
    # Use power law all the way
    lambda <- p@resource_params$lambda
    kappa <- initialNResource(p)[1] * p@w_full[1]^lambda
    initialNResource(p) <- kappa * p@w_full^(-lambda)

    # Switch off all interactions
    p@interaction[] <- 0
    ext_encounter(p) <- getEncounter(p)
    species_params(p)$interaction_resource[] <- 0

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

    return(p)
}
