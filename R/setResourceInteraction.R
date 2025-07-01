#' Set resource interaction to absorb external encounter rate
#'
#' This function sets the resource interaction for each species to absorb as
#' much of the external encounter rate as possible. It sets the resource
#' dynamics and adjusts the resource capacity so that the resource can sustain
#' the extra encounter without changing the current resource abundance or the
#' resource level.
#'
#' If the initial abundances were at steady state, then they will remain at
#' steady state. This is because the total encounter rate is not modified by
#' this function and because the extra mortality on the resource is exactly
#' balanced by the changed resource capacity.
#'
#' @param params A mizer params object
#' @param resource_dynamics Optional. Name of the function that determines the
#'   resource dynamics by calculating the resource spectrum at the next time
#'   step from the current state. If not provided, the function will
#'   use `resource_semichemostat()`.
#' @return The modified mizer params object with increased `interaction_resource`
#'   species parameters and correspondingly decreased external encounter rate.
#' @export
setResourceInteraction <- function(params,
                                   resource_dynamics = "resource_semichemostat") {
    # Save the resource level so we can restore it later
    resource_level <- resource_level(params)
    # If the resource level is exactly 1 then decrease it slightly
    # to make it valid
    resource_level[resource_level == 1] <- 0.999999999

    # Calculate the encounter rate achieved with interaction_resource = 1
    temp_params <- params
    species_params(temp_params)$interaction_resource <- 1
    encounter <- getResourceEncounterRate(temp_params)

    # Then calculate the ratio of the external encounter rate to the
    # encounter rate achieved with interaction_resource = 1
    ratio <- params@ext_encounter / encounter
    ratio[is.nan(ratio)] <- 0
    # The minimum ratio is the maximum that can be absorbed from the external
    # encounter rate
    ratio <- apply(ratio, 1, min)
    absorbed_encounter <- encounter * ratio

    # Increase the resource interaction
    params@species_params$interaction_resource <-
        params@species_params$interaction_resource + ratio

    # Subtract the absorbed encounter rate from the external encounter rate
    params@ext_encounter <- params@ext_encounter - absorbed_encounter

    # Set the resource capacity so that the the steady-state resource
    # abundance stays the same
    params <- setResource(params,
                          resource_level = resource_level,
                          resource_dynamics = resource_dynamics)

    return(params)
}

#' Get resource encounter rate
#'
#' This function returns the resource encounter rate for predators.
#'
#' @param params A mizer params object
#' @return A matrix (species x species) of resource encounter rates
#' @export
getResourceEncounterRate <- function(params) {
    pred_kernel <- getPredKernel(params)
    n_pp <- initialNResource(params)
    phi_resource <- params@species_params$interaction_resource *
        rowSums(sweep(
            pred_kernel, 3, params@dw_full * params@w_full * n_pp,
            "*",
            check.margin = FALSE
        ), dims = 2)
    encounter <- params@search_vol * phi_resource
    return(encounter)
}
