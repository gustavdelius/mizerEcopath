#' Match the gonadic proportion of a species to a given value.
#'
#' This function adjusts the `w_repro_max` parameters to match the gonadic
#' proportion of the species.
#'
#' @param params A MizerParams object
#' @param steady Whether to return the model to a steady state after adjusting
#'  the gonadic proportion. Default is TRUE.
#' @return A MizerParams object with the gonadic proportion matched
#' @family match functions
#' @export
matchGonadicProportionOnce <- function(params, steady = TRUE) {
    assert_that(is.flag(steady))
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }

    # Adjust gonadic production
    current <- getGonadicProduction(params) / getProduction(params)
    ratio <- sp$gonad_proportion / current
    sp <- set_species_param_default(sp, "w_repro_max", sp$w_max)
    sp <- set_species_param_default(sp, "m", 1)
    w_maxratio <- ratio ^ (1 / (sp$n - sp$m))
    sp$w_repro_max <- sp$w_repro_max * w_maxratio
    if (any(sp$w_repro_max < sp$w_mat)) {
        stop("The gonadic proportion leads to a `w_repro_max` smaller than `w_mat`.")
    }
    species_params(params)$w_repro_max <- sp$w_repro_max

    if (steady) {
        # Determine new steady state
        params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    return(params)
}
