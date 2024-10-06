#' Set up parameters for a model with allometric encounter and mortality rates
#'
#' This function sets up a model with allometric encounter and
#' mortality rates for each species and no interaction between the species
#' or with the resource. The model has power-law mortality with a value of 0.4
#' at maturity size and no satiation.
#'
#' @param species_params A data frame with species parameters
#' @param no_w The number of weight bins to use in the model
#' @param lambda The exponent of the Sheldon spectrum
#' @return A MizerParams object
#' @export
newAllometricParams <- function(species_params, no_w = 200, lambda = 2) {

    sp <- species_params
    sp <- set_species_param_default(sp, "n", 0.7)
    sp <- set_species_param_default(sp, "p", sp$n)
    sp <- set_species_param_default(sp, "d", 1 - sp$n)
    max_w <- max(sp$w_max)
    p <- newMultispeciesParams(sp, no_w = no_w, info_level = 0,
                               # extend resource over entire size range
                               max_w = max_w,
                               w_pp_cutoff = max_w * (1 + 1e-9),
                               lambda = lambda,
                               resource_dynamics = "resource_constant")
    sp <- p@species_params

    # Switch off all interactions
    p <- makeNoninteracting(p)

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

    # Turn off satiation
    p <- setFeedingLevel(p, 0)

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
isAllometric <- function(params, tol = 1e-6) {
    sp <- params@species_params
    sp <- set_species_param_default(sp, "d", sp$n - 1)
    m <- getMort(params)
    e <- getEncounter(params)
    for (species in sp$species) {
        spc <- sp[species, ]
        mc <- m[species, ] / params@w ^ spc$d
        if (anyNA(mc) || any(is.nan(mc)) || any(mc < 0) || all(mc == 0)) {
            stop("The model has invalid mortality rates.")
        }
        if (abs((max(mc) - min(mc)) / max(mc)) > tol) {
            return(FALSE)
        }
        ec <- e[species, ] / params@w ^ spc$n
        if (anyNA(ec) || any(is.nan(ec)) || any(ec < 0) || all(ec == 0)) {
            stop("The model has invalid encounter rates.")
        }
        if (abs((max(ec) - min(ec)) / max(ec)) > tol) {
            return(FALSE)
        }
    }
    return(TRUE)
}
