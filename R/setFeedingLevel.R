#' Set constant feeding level for each species in a non-interacting model
#'
#' This function adjusts the maximum intake rate in order to achieve the desired
#' feeding level for each species. It also adjust the external encounter rate
#' so that the consumption rate does not change in spite of the change of
#' feeding level.
#'
#' This function only works for models where all encounter is external
#' encounter. Try calling `makeNoninteracting()` first.
#'
#' @param params A MizerParams object
#' @param feeding_level The feeding level to set for each species. If missing,
#'   the feeding level is taken from the `f0` column in the species_params
#'   object. If `f0` is not present, an error is thrown. If given, it must be
#'   a vector with one value for each species or a single value that is used
#'   for all species.
#' @return A MizerParams object with the feeding level set
#' @export
setFeedingLevel <- function(params, feeding_level) {
    sp <- params@species_params
    # Check that params describes a non-interacting model
    if (!isTRUE(all.equal(getEncounter(params), getExtEncounter(params)))) {
        stop("This function only works for models where all encounter is external encounter. Try calling `makeNoninteracting()` first.")
    }
    if (missing(feeding_level)) {
        if (hasName(sp, "f0")) {
            feeding_level <- sp$f0
        } else {
            stop("You need to supply the desired feeding_level.")
        }
    }
    assert_that(is.numeric(feeding_level))
    # If feeding_level is a single value, make it a vector
    if (length(feeding_level) == 1) {
        feeding_level <- rep(feeding_level, nrow(sp))
    }
    # Check that feeding_level is a vector with one value for each species
    if (length(feeding_level) != nrow(sp)) {
        stop("The feeding_level vector has the wrong length")
    }
    # Check that the feeding level is in the correct range
    if (any(feeding_level < 0 | feeding_level >= 1)) {
        stop("Feeding level must be positive and strictly less than 1.")
    }

    # Save feeding level in species params
    # Don't use `given_species_params<-` because we do not want to trigger
    # a recalculation of anything
    params@given_species_params$f0 <- feeding_level
    params@species_params$f0 <- feeding_level

    # Adjust maximum intake rate to achieve the desired feeding level
    E_old <- getExtEncounter(params)
    f_old <- getFeedingLevel(params)
    intake_max(params) <- (1 - f_old) / feeding_level * E_old
    # Division by zero may have produced NaNs, replace them with Inf
    intake_max(params)[is.nan(intake_max(params))] <- Inf

    # Adjust the encounter rate so that consumption does not change
    ext_encounter(params) <- (1 - f_old) / (1 - feeding_level) * E_old

    params@time_modified <- lubridate::now()
    return(params)
}
