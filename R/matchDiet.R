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
#' @param diet_matrix The diet matrix to match. Rows correspond to predator
#'   species and columns to prey species. Values give the diet composition of
#'   each predator. Any prey columns whose names do not match species in
#'   `params` are summed into an "other" category. Rows are normalised to sum
#'   to 1, so the entries need not be pre-normalised proportions.
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

    # Make the model non-interacting
    params <- makeNoninteracting(params)

    # Get diet matrix when interaction matrix is all 1, keeping the total
    # encounter rate unchanged. Setting interaction_matrix[] <- 1 would add
    # species-species encounter on top of the already-full ext_encounter, so we
    # subtract that excess to hold the total constant. This makes the diet
    # exactly linear in theta, so that theta = D / E recovers the original
    # interaction strengths when D is the current diet matrix.
    interaction_matrix(params)[] <- 1
    params_E <- params
    ext_encounter(params_E) <- 2 * ext_encounter(params_E) - getEncounter(params_E)
    E <- getDietMatrix(params_E, min_w_pred = min_w_pred,
                       max_w_pred = max_w_pred)[1:no_sp, 1:no_sp, drop = FALSE]
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
    if (any(diet_matrix < 0)) {
        stop("`diet_matrix` contains negative values.")
    }
    if (any(rowSums(diet_matrix) == 0)) {
        stop("According to the diet matrix, some species do not eat anything.")
    }
    invisible(NULL)
}

#' Convert a proportional diet matrix to absolute consumption
#'
#' This internal helper function converts a diet matrix into a matrix of
#' absolute consumption (i.e. biomass consumed per year), based on the species'
#' total consumption rates. Any prey columns not matching model species are
#' summed into an "other" category, rows are normalised to sum to 1, and then
#' each row is scaled by the total consumption of the corresponding predator.
#' The result is restricted to the modelled species on both axes.
#'
#' This functions is intended for internal use and is called by `matchDiet()`.
#'
#' @inheritParams getDietMatrix
#' @param diet_matrix A numeric matrix where rows correspond to predator species
#'   and columns to prey species. Any prey columns whose names do not match
#'   species in `params` are summed into a single "other" column. Rows are
#'   normalised to sum to 1. Should include all model species as row names.
#' @return A numeric matrix with dimensions (n_species x n_species), giving
#'   absolute consumption rates (g/year) for each predator–prey pair.
#' @keywords internal
convertDietMatrix <- function(diet_matrix, params, min_w_pred, max_w_pred) {
    sp <- params@species_params
    no_sp <- nrow(sp)
    checkDietMatrix(diet_matrix)

    # Sum all prey columns not matching model species into "other"
    species_prey_cols <- intersect(colnames(diet_matrix), sp$species)
    other_prey_cols <- setdiff(colnames(diet_matrix), sp$species)
    other_sum <- if (length(other_prey_cols) > 0) {
        rowSums(diet_matrix[, other_prey_cols, drop = FALSE])
    } else {
        rep(0, nrow(diet_matrix))
    }
    diet_matrix <- cbind(
        diet_matrix[, species_prey_cols, drop = FALSE],
        other = other_sum
    )

    # Normalise rows to sum to 1
    diet_matrix <- diet_matrix / rowSums(diet_matrix)

    # Convert diet matrix from proportions to absolute consumption
    Q <- getConsumption(params,
                        min_w_pred = min_w_pred,
                        max_w_pred = max_w_pred)
    if (!all(sp$species %in% rownames(diet_matrix))) {
        stop("diet_matrix does not include all model species as rows.")
    }
    dm <- diet_matrix[sp$species, , drop = FALSE] * Q
    # Keep only the species-on-species part
    D <- dm[, sp$species, drop = FALSE]
    return(D)
}
