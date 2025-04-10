#' Match the diet matrix of the model to an observed diet matrix
#'
#' This function adjusts the interaction matrix of the model so that its diet
#' matrix matches the given diet matrix (typically derived using`reduceEcopathDiet()`).
#'
#' This function first converts the input diet matrix from proportions to
#' absolute consumption using species-specific consumption rates. It then
#' calculates the interaction strengths required to reproduce this diet and
#' applies them using `makeInteracting()`. External mortality and encounter
#' terms are adjusted to preserve the model’s steady state.
#'
#' If `min_w_pred` or `max_w_pred` are provided, the diet matrix is assumed to
#' represent consumption by individuals within that size range
#'
#' If matching the diet requires interaction strengths that would imply
#' negative external mortality or encounter rates for any species, warnings
#' will be issued by `makeInteracting()`. These indicate that the model
#' assumptions may be biologically unrealistic or that the input diet matrix
#' is incompatible with other model parameters.
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
    D <- convertDietMatrix(diet_matrix = diet_matrix,
                           params = params,
                           min_w_pred = min_w_pred,
                           max_w_pred = max_w_pred)

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
#' Checks that a diet matrix is valid for use in mizer.
#' Ensures there are no missing values (NA/NaN), and that each predator species
#' has a non-zero total diet.
#'
#' @param diet_matrix A numeric matrix, typically the output of
#'   `reduceEcopathDiet()`, with predators as rows and prey as columns.
#' @return Invisibly returns `NULL`. Throws an error if the matrix is invalid.
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

#' Convert a proportional diet matrix to absolute consumption
#'
#' This internal helper function converts a diet matrix expressed in proportions
#' (i.e. each row summing to 1) into a matrix of absolute consumption
#' (i.e. biomass consumed per year), based on the species' total consumption
#' rates. It also restricts the matrix to include only the modelled species
#' (i.e. drops the "other" column and any unmatched rows). The conversion is done
#' by multiplying each row of the input diet matrix by the total consumption of the
#' corresponding predator species. This gives the absolute biomass consumed by each
#' predator from each prey.
#'
#' This functions is intended for internal use and is called by `matchDiet()`.
#'
#' @inheritParams getDietMatrix
#' @param diet_matrix A numeric matrix where rows correspond to predator species
#'   and columns to prey species (plus an optional "other" column). Rows must sum
#'   to 1. Should include all model species as both row and column names.
#' @return A numeric matrix with dimensions [n_species x n_species], giving
#'   absolute consumption rates (g/year) for each predator–prey pair.
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
