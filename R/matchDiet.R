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
#' size preference. We multiply the resulting encounter rate of these small
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
#' @family match functions
#' @export
matchDiet <- function(params, diet_matrix, centering = 0,
                             w_prey_cutoff = 1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    no_sp <- nrow(sp)

    # This function is to be used on a model with infinite maximum intake rate
    if (any(intake_max(params) != Inf)) {
        stop("This function is for models with infinite maximum intake rate.")
    }

    # Check diet matrix ----
    if (!is(diet_matrix, "matrix")) {
        stop("`diet_matrix` must be a matrix.")
    }
    if (!is.numeric(diet_matrix)) {
        stop("`diet_matrix` must be numeric.")
    }
    if (any(is.nan(diet_matrix))) {
        stop("`diet_matrix` contains NaNs.")
    }
    if (any(is.na(diet_matrix))) {
        stop("`diet_matrix` contains NAs.")
    }
    if (any(rowSums(diet_matrix) == 0)) {
        stop("According to the diet matrix, some species do not eat anything.")
    }

    # Convert diet matrix from proportions to absolute consumption
    Q <- sp$consumption_observed
    dm <- diet_matrix * Q / rowSums(diet_matrix)
    # Keep only the part corresponding to species
    D <- dm[sp$species, sp$species]
    if (nrow(D) != no_sp || ncol(D) != no_sp) {
        stop("diet_matrix does not include all model species.")
    }

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
    interaction_matrix(params)[] <- 1
    dmm <- getDietMatrix(params)
    ratios <- rowSums(dm) / rowSums(dmm)
    species_params(params)$gamma <- species_params(params)$gamma * ratios
    params@search_vol <- params@search_vol * ratios
    E <- getDietMatrix(params)
    interaction_matrix(params)[] <- 0

    # If the encounter matrix has zeros where the Ecopath diet does not, then
    # there is no way to match the diet matrix. This is unlikely to happen.
    if (any(E == 0 & D > 0)) {
        stop("Diet matrix cannot be matched because the mizer predation kernel does lead to zero predator-prey interaction between some species.")
    }

    theta <- D / E
    theta[is.nan(theta)] <- 1  # Avoid NaNs when E is zero

    params <- makeInteractive(params, theta)

    if (centering == 0) {
        return(params)
    }
    stop("Centering not implemented yet.")

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

    params <- makeInteractive(params, theta)

    return(params)
}
