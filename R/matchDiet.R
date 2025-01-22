#' Match the diet matrix of the model to the Ecopath diet matrix
#'
#' This function adjusts the interaction matrix of the model so that its diet
#' matrix matches the Ecopath diet matrix. It then adjusts the external
#' encounter and mortality rates so that the steady state is not changed.
#'
#' The function will raise a warning for any species that requires a negative
#' external mortality rate or negative external encounter rate.
#'
#' @param params A MizerParams object
#' @param diet_matrix The Ecopath diet matrix as produced by
#'   `reduceEcopathDiet()`.
#'
#' @return A MizerParams object with the interaction matrix set so that the
#'   desired diet is achieved.
#' @family match functions
#' @export
matchDiet <- function(params, diet_matrix, centering = 0,
                             w_prey_cutoff = 1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    no_sp <- nrow(sp)

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

    # Get diet matrix when interaction matrix is all 1. This is also the
    # Encounter matrix because we have switched off satiation.
    interaction_matrix(params)[] <- 1
    E <- getDietMatrix(params)[1:no_sp, 1:no_sp]
    interaction_matrix(params)[] <- 0

    # If the encounter matrix has zeros where the Ecopath diet does not, then
    # there is no way to match the diet matrix. This is unlikely to happen.
    if (any(E == 0 & D > 0)) {
        stop("Diet matrix cannot be matched because the mizer predation kernel leads to zero predator-prey interaction between some species.")
    }

    theta <- D / E
    theta[is.nan(theta)] <- 1  # Avoid NaNs when E is zero

    params <- makeInteractive(params, theta)

    return(params)
}
