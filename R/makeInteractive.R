#' Make species interact with given interaction matrix
#'
#' This function should be called with a non-interacting model (where there is
#' only external mortality and external encounter) and an interaction matrix.
#' The external encounter and mortality rates will then be decreased by the
#' amount of the predation encounter and mortality that is created by that
#' interaction matrix.
#'
#' @param params A MizerParams object
#' @param interaction An interaction matrix where the entry in row i and
#'   column j determines how strongly species i predates on species j.
#' @return The modified MizerParams object
#' @export
makeInteractive <- function(params, interaction) {
    if (any(params@interaction != 0) || any(params@interaction_resource != 0)) {
        stop("This function should be called with a non-interacting model.")
    }

    interaction_matrix(params) <- interaction

    # Temporarily switch off external encounter to measure predation encounter
    pe <- params
    ext_encounter(pe)[] <- 0
    pred_encounter <- getEncounter(pe)

    # Reduce external encounter rate by predation encounter rate
    new_ext_encounter <- ext_encounter(params) - pred_encounter
    # Check that predator encounter rate is less than total encounter rate
    # for all species up to at least their largest size
    for (i in 1:no_sp) {
        if (any(new_ext_encounter[i, ] < 0) &&
            (min(params@w[new_ext_encounter[i, ] < 0]) <=
             max(params@w[initialN(params)[i, ] > 1e-10]))) {
            warning("Negative external encounter rate required for ",
                    params@species_params$species[i])
        }
    }
    # Don't allow negative encounter rates
    new_ext_encounter[new_ext_encounter < 0] <- 0
    ext_encounter(params) <- new_ext_encounter

    ## Change external mortality rate
    pred_mort <- getPredMort(params)
    new_ext_mort <- ext_mort(params) - pred_mort
    # Check that predator mortality rate is less than total mortality rate
    # for all species up to at least their largest size
    for (i in 1:no_sp) {
        if (any(new_ext_mort[i, ] < 0) &&
            (min(params@w[new_ext_mort[i, ] < 0]) <=
             max(params@w[initialN(params)[i, ] > 0]))) {
            warning("Negative external mortality rate required for ",
                    params@species_params$species[i])
        }
    }
    # Don't allow negative mortality rates
    new_ext_mort[new_ext_mort < 0] <- 0
    ext_mort(params) <- new_ext_mort

    return(params)
}
