#' Enable predator–prey interactions using an interaction matrix
#'
#' This function should be called with a non-interacting model (where there is
#' only external mortality and external encounter) and an interaction matrix.
#' The external encounter and mortality rates will then be decreased by the
#' amount of the predation encounter and mortality that is created by that
#' interaction matrix.
#'
#' If predation would exceed total mortality or encounter for a species, a warning is issued
#' and negative values are truncated to zero. This can cause the resulting model to
#' diverge from the steady state.
#'
#' Typically used internally by [matchDiet()]
#'
#' @param params A MizerParams object
#' @param interaction An interaction  matrix where entry *(i, j)* gives the strength
#' with which species *i* preys on *j*.
#' @param interaction_resource A numeric vector of length equal to the number of
#'   species giving the interaction strength of each species with the resource.
#'   Defaults to zero (no resource interaction). Pass the original values when
#'   reinstating resource interactions that were zeroed by [makeNoninteracting()].
#' @return The modified MizerParams object
#' @export
makeInteracting <- function(params, interaction,
                            interaction_resource = NULL) {
    if (any(params@interaction != 0) ||
        any(params@species_params$interaction_resource != 0)) {
        stop("This function should be called with a non-interacting model.")
    }
    sp <- params@species_params
    no_sp <- nrow(sp)

    if (is.null(interaction_resource)) {
        interaction_resource <- rep(0, no_sp)
    }

    interaction_matrix(params) <- interaction
    params@species_params$interaction_resource <- interaction_resource

    # Temporarily switch off external encounter to measure predation and
    # resource encounter
    pe <- params
    pe@ext_encounter[] <- 0
    pred_encounter <- getEncounter(pe)

    # Reduce external encounter rate by predation encounter rate
    new_ext_encounter <- params@ext_encounter - pred_encounter
    # Check that predator encounter rate is less than total encounter rate
    # for all species up to at least their largest size
    for (i in 1:no_sp) {
        if (any(new_ext_encounter[i, ] < 0) &&
            (min(params@w[new_ext_encounter[i, ] < 0]) <=
             max(params@w[initialN(params)[i, ] > 1e-10]))) {
            warning("Implied negative external encounter rate for ",
                    params@species_params$species[i],
                    ". Setting it to zero. This may alter the steady state.")
        }
    }
    # Don't allow negative encounter rates
    new_ext_encounter[new_ext_encounter < 0] <- 0
    params@ext_encounter <- new_ext_encounter

    ## Change external mortality rate
    pred_mort <- getPredMort(params)
    new_ext_mort <- params@mu_b - pred_mort
    # Check that predator mortality rate is less than total mortality rate
    # for all species up to at least their largest size
    for (i in 1:no_sp) {
        if (any(new_ext_mort[i, ] < 0) &&
            (min(params@w[new_ext_mort[i, ] < 0]) <=
             max(params@w[initialN(params)[i, ] > 0]))) {
            warning("Implied negative external mortality rate for ",
                    params@species_params$species[i],
                    ". Setting it to zero. This may alter the steady state.")
        }
    }
    # Don't allow negative mortality rates
    new_ext_mort[new_ext_mort < 0] <- 0
    params@mu_b <- new_ext_mort

    return(params)
}
