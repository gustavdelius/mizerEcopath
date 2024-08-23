getEcopathRatios <- function(params) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(sp, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    # Calculate discrepancy in consumption
    Qratio <- sp$ecopath_consumption / getConsumption(params)
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

    df <- data.frame(species = sp$species, Qratio = Qratio,
                     Pratio = Pratio, Gratio = Gratio,
                     Cratio = Cratio, Zratio = Zratio)
    return(df)
}

#' Match the model to the Ecopath data
#'
#' This function adjusts the model parameters to match the Ecopath data. It does
#' this by adjusting the consumption, production, and yield parameters.
#'
#' The `match` function calls `matchConsumption` and
#' `matchProduction` repeatedly until both consumption and production
#' matches within the tolerance specified by the `tol` parameter. The function
#' will return a warning if the maximum number of iterations is reached without
#' converging.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance for the match. Default is 0.1.
#' @param max_iter The maximum number of iterations. Default is 100.
#'
#' @return A MizerParams object matched to the Ecopath data
#' @export
match <- function(params, tol = 0.1, max_iter = 100) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }

    params <- matchConsumption(params, tol)

    for (i in seq_len(max_iter)) {
        # Break if tolerance achieved
        if (isMatched(params, tol)) {
            break
        }
        params <- params |>
            matchConsumptionOnce() |>
            matchProductionOnce()
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Converged in ", i - 1, " iterations.")
    }
    params
}

#' @rdname match
isMatched <- function(params, tol = 0.1) {
    isConsumptionMatched(params, tol) &&
        isProductionMatched(params, tol)
}
