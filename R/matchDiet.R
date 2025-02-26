#' Match the diet matrix of the model to the Ecopath diet matrix
#'
#' This function adjusts the interaction matrix of the model so that its diet
#' matrix matches the Ecopath diet matrix. It then adjusts the external
#' encounter and mortality rates so that the steady state is not changed.
#'
#' If `min_w_pred` or `max_w_pred` are set, then the diet matrix is assumed to
#' represent the flow of biomass into the part of the predator species
#' population with sizes in the specified range.
#'
#' The function will raise a warning for any species that requires a negative
#' external mortality rate or negative external encounter rate.
#'
#' @param params A MizerParams object
#' @param diet_matrix The Ecopath diet matrix as produced by
#'   `reduceEcopathDiet()`.
#' @inheritParams getDietMatrix
#'
#' @return A MizerParams object with the interaction matrix set so that the
#'   desired diet is achieved.
#' @family match functions
#' @export
matchDiet <- function(params, diet_matrix,
                      min_w_pred = 0, max_w_pred = Inf) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    no_sp <- nrow(sp)

    # Convert the diet matrix to absolute consumption
    D <- convertDietMatrix(diet_matrix,
                           min_w_pred = min_w_pred,
                           max_w_pred = max_w_pred,
                           sp, no_sp)

    # Make the interaction matrix non-interacting
    params <- makeNoninteracting(params)

    # Get diet matrix when interaction matrix is all 1. This is also the
    # Encounter matrix because we have switched off satiation.
    interaction_matrix(params)[] <- 1
    E <- getDietMatrix(params, min_w_pred = min_w_pred,
                       max_w_pred = max_w_pred)[1:no_sp, 1:no_sp]
    interaction_matrix(params)[] <- 0

    # If the encounter matrix has zeros where the Ecopath diet does not, then
    # there is no way to match the diet matrix. This is unlikely to happen.
    if (any(E == 0 & D > 0)) {
        stop("Diet matrix cannot be matched because the mizer predation kernel ",
             "leads to zero predator-prey interaction between some species.")
    }

    theta <- D / E
    theta[is.nan(theta)] <- 1  # Avoid NaNs when E is zero

    params <- makeInteracting(params, theta)

    return(params)
}

#' Check if a diet matrix is valid
#'
#' This function performs validation checks on a diet matrix to ensure it meets
#' the requirements for use in mizer.
#'
#' @param diet_matrix The diet matrix to check
#' @return NULL invisibly. Throws an error if any check fails.
#' @keywords internal
checkDietMatrix <- function(diet_matrix) {
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
    invisible(NULL)
}

#' Convert diet matrix to absolute consumption
#'
#' This function converts a diet matrix from proportions to absolute consumption
#' and extracts the part corresponding to species.
#'
#' @param diet_matrix The diet matrix to convert
#' @inheritParams getDietMatrix
#' @param sp Species parameters data frame
#' @param no_sp Number of species
#' @return The converted diet matrix
#' @keywords internal
convertDietMatrix <- function(diet_matrix, params, min_w_pred, max_w_pred) {
    sp <- params@species_params
    no_sp <- nrow(sp)
    checkDietMatrix(diet_matrix)
    # Convert diet matrix from proportions to absolute consumption
    Q <- getConsumption(params,
                        min_w_pred = min_w_pred,
                        max_w_pred = max_w_pred)
    dm <- diet_matrix * Q / rowSums(diet_matrix)
    # Keep only the part corresponding to species
    D <- dm[sp$species, sp$species]
    if (nrow(D) != no_sp || ncol(D) != no_sp) {
        stop("diet_matrix does not include all model species.")
    }
    return(D)
}
