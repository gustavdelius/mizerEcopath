#' Set resource interaction to absorb external encounter rate
#'
#' This function sets the resource interaction for each species to absorb as
#' much of the external encounter rate as possible. It adjusts the resource
#' capacity so that the resource can sustain the extra encounter without
#' changing the current resource abundance or the resource level.
#'
#' If the initial abundances were at steady state, then they will remain at
#' steady state. This is because the total encounter rate is not modified by
#' this function and because the extra mortality on the resource is exactly
#' balanced by the changed resource capacity.
#'
#' Optionally, you can specify the resource dynamics to use (e.g.,
#' "resource_semichemostat" or "resource_constant") by setting the
#' resource_dynamics argument. This allows you to switch between dynamic and
#' constant resource regimes directly through this function. For example:
#'   params <- setResourceInteraction(params, resource_dynamics = "resource_semichemostat")
#' will set the resource dynamics to semichemostat.
#'
#' Note: The function ensures that the resource level passed to setResource()
#' is strictly between 0 and 1, as required by mizer. If any value of
#' resource_level is greater than or equal to 1, it is set to just below 1
#' (0.99999999). If any value is less than or equal to 0, it is set to a small
#' positive value (1e-8). This prevents errors and ensures biological realism.
#'
#' @param params A mizer params object
#' @param resource_dynamics Optional. Character string specifying the resource dynamics to use. If NULL (default), does not change the current resource dynamics.
#' @return The modified mizer params object with increased `interaction_resource`
#'   species parameters and correspondingly decreased external encounter rate.
#' @export
setResourceInteraction <- function(params, resource_dynamics = NULL) {
    # Optionally set resource dynamics if specified by the user
    if (!is.null(resource_dynamics)) {
        resource_dynamics(params) <- resource_dynamics
    }
    # Save the resource level so we can restore it later
    # The resource_level must be strictly between 0 and 1 for setResource().
    # Values >= 1 are set to just below 1, and values <= 0 are set to a small positive value.
    resource_level <- resource_level(params)
    resource_level[resource_level >= 1] <- 0.99999999
    resource_level[resource_level <= 0] <- 1e-8

    # Calculate the encounter rate achieved with interaction_resource = 1
    temp_params <- params
    species_params(temp_params)$interaction_resource <- 1
    encounter <- getResourceEncounterRate(temp_params)

    # Then calculate the ratio of the external encounter rate to the
    # encounter rate achieved with interaction_resource = 1
    ratio <- params@ext_encounter / encounter
    # This can have zeros if encounter is zero, which can happen at very large
    # sizes that we are not interested in. Set these zeros to Inf so they do not
    # affect the minimum calculation below.
    ratio[ratio == 0] <- Inf
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
                          resource_level = resource_level)

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
