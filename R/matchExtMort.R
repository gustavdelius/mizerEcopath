#' Match the external mortality to that required by the ecopath production rate.
#'
#' Adjust the external mortality so that the biomass loss due to mortality
#' (including the fisheries yield) matches the somatic production.
#'
#' @param params A MizerParams object
#' @param steady Whether to return the model to a steady state after adjusting
#'  the external mortality. Default is TRUE.
#' @return A MizerParams object with the external mortality matched
#' @family match functions
#' @export
matchExtMortOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    gp <- validGearParams(params@gear_params, params@species_params)
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    # Adjust external mortality so that the loss due to mortality matches
    # the somatic production: M0B + C = Ps
    current <- getM0B(params)
    desired <- getSomaticProduction(params) - gp$yield_observed
    if (any(desired < 0)) {
        stop("Somatic production does not cover observed yield.")
    }
    Zratio <- desired / current
    sel <- !is.na(Zratio)
    ext_mort(params)[sel, ] <- ext_mort(params)[sel, ] * Zratio

    if (steady) {
        # Determine new steady state
        params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    params@time_modified <- lubridate::now()
    return(params)
}
