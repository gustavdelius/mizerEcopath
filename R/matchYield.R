#' Match the observed yield by adjusting catchability
#'
#' This function adjusts the catchability parameters of the model to match the
#' observed yield.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' @param params A MizerParams object
#' @param steady Whether to return the model to a steady state after adjusting
#'   the catchability. Default is TRUE.
#' @return A MizerParams object with the catchability adjusted to match the
#'   observed yield
#' @family match functions
#' @export
matchYieldOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    # Adjust fishing mortality
    # TODO: I am just lazy here at the moment. Code for making yield match
    # for each gear and each species will probably require looping
    if (nrow(gp) != nrow(sp) ||
        !all(gp$species == sp$species)) {
        stop("This code currently requires a single gear for each species.")
    }
    Cratio <- gp$yield_observed / getYield(params)
    sel <- !is.na(Cratio)
    gp$catchability[sel] <- gp$catchability[sel] * Cratio
    gear_params(params)$catchability <- gp$catchability

    if (steady) {
        # Determine new steady state
        params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    return(params)
}
