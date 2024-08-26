#' Match the diet matrix of the model to the Ecopath diet matrix
#'
#' This function adjusts the interaction matrix of the model so that its diet
#' matrix matches the Ecopath diet matrix. It then adjusts the external
#' encounter and mortality rates so that the steady state is not changed.
#'
#' The function will raise an error if negative external encounter rates are
#' required for any species and issues a warning for any species that requires
#' negative external mortality rates.
#'
#' Note that the function does not currently work properly when the mizer model
#' makes use of additional ecosystem components that add to the encounter rate.
#'
#' @section Centering:
#' Ecopath diet matrices often contain many zeros, which would mean that it is
#' impossible for two species to interact. This is because the Ecopath diet
#' matrix does not account for small prey items that do not contribute
#' significantly to a predator's diet. The `centering` parameter allows one to
#' experiment with adding an estimate of the diet of these small prey to the
#' diet matrix.
#'
#' For this purpose we assume that the Ecopath diet matrix did not account for
#' prey items with a weight below `w_cutoff`. To estimate their contribution to
#' the diet we make the assumption that for such small prey items the predator
#' does not distinguish between species and feeds solely according to its prey
#' size preference. We multipoly the resulting encounter rate of these small
#' prey items by the `centering` argument and add it to the Ecopath diet matrix
#' before matching.
#'
#' @param params A MizerParams object
#' @param diet_matrix The Ecopath diet matrix as produced by
#'   `reduceEcopathDiet()`.
#' @param centering The centering parameter, see "Centering" section below for
#'   details.
#' @param w_prey_cutoff See "Centering" section below for details.
#'
#' @return A MizerParams object with the diet matrix matched.
#' @export
matchDiet <- function(params, diet_matrix, centering = 0,
                             w_prey_cutoff = 1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    no_sp <- nrow(params@species_params)

    # This function is to be used on a model with infinite maximum intake rate
    if (any(intake_max(params) != Inf)) {
        stop("This function is for models with infinite maximum intake rate.")
    }

    if (!is(diet_matrix, "matrix")) {
        stop("diet_matrix must be a matrix.")
    }
    if (nrow(diet_matrix) != no_sp ||
        ncol(diet_matrix) != no_sp + 1) {
        stop("diet_matrix has the wrong dimensions.")
    }
    # Convert diet matrix from proportions to absolute consumption
    Q <- params@species_params$ecopath_consumption
    dm <- diet_matrix * Q / rowSums(diet_matrix)
    # Drop the "other" column
    D <- dm[, -ncol(dm)]

    params <- makeNoninteracting(params)

    # Set resource to be in line with fish
    total <- colSums(initialN(params))
    fish_sel <- params@w_full >= params@w[1]
    ratio <- max(total / initialNResource(params)[fish_sel])
    params <- mizerExperimental::scaleDownBackground(params, 1/ratio)

    # Get diet matrix when interaction matrix is all 1. This is also the
    # Encounter matrix because we have switched off satiation.
    # Set `gamma` so that this interaction matrix produces the desired
    # consumption
    p <- params
    interaction_matrix(p)[] <- 1
    species_params(p)$interaction_resource <- 1
    dmm <- getDietMatrix(p)
    species_params(params)$gamma <- species_params(params)$gamma *
        rowSums(dm) / rowSums(dmm)

    # Diet matrix restricted to large prey
    E <- getDietMatrix(p, min_w_prey = w_prey_cutoff)[, 1:no_sp]
    # Diet matrix restricted to small prey
    S <- getDietMatrix(p, max_w_prey = w_prey_cutoff)[, 1:no_sp]

    # If more of the encounter comes from small prey than from large prey then
    # we are more free to keep that entry of the interaction matrix close to 1
    centering <- centering * S / (E + S)

    # If the encounter matrix has zeros where the Ecopath diet does not, then
    # there is no way to match the diet matrix. This is unlikely to happen.
    if (any(E == 0 & D > 0)) {
        stop("Diet matrix cannot be matched because the mizer predation kernel does lead to zero predator-prey interaction between some species.")
    }

    # We want to set the interaction matrix so that the diet E of non-larval
    # prey in the mizer model is close to the ecopath diet D. In other we want
    # the following ratios to be close to 1:
    ratios <- D / E
    ratios[is.nan(ratios)] <- 1  # Avoid NaNs when E is zero
    # If at the same time we want to keep species preferences approximately
    # equal we would use:
    # theta <- (ratios + centering * rowMeans(ratios)) / (1 + centering)
    # The problem with that is that when the ecopath diet matrix entry is zero
    # then the corresponding entry in the interaction matrix will be zero.
    # We will instead penalise deviations from 1 in theta
    theta <- (ratios + centering) / (1 + centering)
    interaction_matrix(params) <- theta

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
            stop("Negative external encounter rate required for ",
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
