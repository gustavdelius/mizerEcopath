#' Match the production of the model to the Ecopath production
#'
#' This function adjusts the parameters of the model to match the
#' Ecopath production. Specifically it adjusts:
#' - `ks` to get the desired respiration,
#' - `w_repro_max` to get the desired gonadic proportion,
#' - `catchability` to get the desired yield, and
#' - `ext_mort` to get the desired mortality.
#'
#' Note that adjusting the `catchability` parameter to match the yield may not
#' lead to a convergent iteration because it currently assumes that increasing
#' the catchability will increase the yield. This may not be the case if the
#' increased mortality on large individuals truncates the size spectrum too
#' much.
#'
#' The `matchProduction` function calls `matchProductionOnce()`
#' repeatedly until the production matches within the tolerance specified by
#' the `tol` parameter.
#' The function will return a warning if the maximum number of iterations is
#' reached without converging.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance for the match. Default is 0.1.
#' @param max_iter The maximum number of iterations. Default is 10.
#'
#' @return A MizerParams object with the production matched
#' @family match functions
#' @export
matchProduction <- function(params, tol = 0.1, max_iter = 10) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    for (i in seq_len(max_iter)) {
        # Break if tolerance achieved
        if (isProductionMatched(params, tol)) {
            break
        }
        params <- matchProductionOnce(params)
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Production converged in ", i - 1, " iterations.")
    }
    params
}

#' @rdname matchProduction
#' @param steady Whether to return the model to a steady state after adjusting
#'   the production. Default is TRUE.
#' @export
matchProductionOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    params <- matchRespirationOnce(params, steady = FALSE)

    params <- matchGonadicProportionOnce(params, steady = FALSE)

    params <- matchYieldOnce(params, steady = FALSE)

    params <- matchExtMortOnce(params, steady = FALSE)

    if (steady) {
        # Determine new steady state
        params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    return(params)
}

#' @rdname matchProduction
isProductionMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    # Calculate discrepancy in production
    Pratio <- sp$ecopath_production / getSomaticProduction(params)
    # Calculate discrepancy in gonadic production
    Gratio <- getGonadicProduction(params) / getProduction(params) /
        sp$gonad_proportion
    # Calculate discrepancy in yields
    Cratio <- gp$yield_observed / getYield(params)
    # Calculate discrepancy in mortality
    current <- getM0B(params) + getM2B(params)
    desired <- getSomaticProduction(params) - gp$yield_observed
    Zratio <- desired / current

    return(max(abs(Cratio - 1)) < tol &&
               max(abs(Zratio - 1)) < tol &&
               max(abs(Pratio - 1)) < tol &&
               max(abs(Gratio - 1)) < tol)
}
