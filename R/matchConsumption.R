#' Match the consumption of the model to the Ecopath consumption
#'
#' This function adjusts the consumption parameters of the model to match the
#' Ecopath consumption It does this by adjusting the `gamma`, `h`, `ks` and
#' `ext_mort` parameters.
#'
#' The `matchConsumption` function calls `matchConsumptionOnce()`
#' repeatedly until the production matches within the tolerance specified by
#' the `tol` parameter.
#' The function will return a warning if the maximum number of iterations is
#' reached without converging.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance for the match. Default is 0.1.
#' @param max_iter The maximum number of iterations. Default is 10.
#'
#' @return A MizerParams object with the consumption matched
#' @export
matchConsumption <- function(params, tol = 0.1, max_iter = 10) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    if (!hasName(sp, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }

    params |> steadySingleSpecies() |>
        matchBiomasses() |> steadySingleSpecies()

    for (i in seq_len(max_iter)) {
        if (isConsumptionMatched(params, tol)) {
            break
        }
        params <- matchConsumptionOnce(params)
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Consumption converged in ", i - 1, " iterations.")
    }
    params
}

#' @rdname matchConsumption
matchConsumptionOnce <- function(params) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    if (!hasName(sp, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }

    Q <- getConsumption(params)
    Qratio <- sp$ecopath_consumption / Q
    sel <- !is.na(Qratio)

    # Adjust consumption
    species_params(params)$gamma[sel] <- species_params(params)$gamma[sel] * Qratio
    species_params(params)$h[sel] <- species_params(params)$h[sel] * Qratio
    species_params(params)$ks[sel] <- species_params(params)$ks[sel] * Qratio
    ext_mort(params)[sel, ] <- ext_mort(params)[sel, ] * Qratio
    ext_encounter(params)[sel, ] <- ext_encounter(params)[sel, ] * Qratio

    # Return new steady state
    params |> steadySingleSpecies() |>
        matchBiomasses() |> steadySingleSpecies()
}

#' @rdname matchConsumption
isConsumptionMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    if (!hasName(sp, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }

    Q <- getConsumption(params)
    Qratio <- sp$ecopath_consumption / Q

    return(max(abs(Qratio - 1)) < tol)
}
