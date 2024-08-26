matchExtMortOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    gp <- validGearParams(params@gear_params)
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
    if (!is.null(comment(ext_mort(params)))) {
        warning("This function has rescaled the external mortality that was set manually.")
    }
    ext_mort(params)[sel, ] <- ext_mort(params)[sel, ] * Zratio

    if (steady) {
        # Determine new steady state
        params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    return(params)
}
