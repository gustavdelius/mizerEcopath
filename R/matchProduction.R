#' Match the production of the model to the Ecopath production
#'
#' This function adjusts the parameters of the model to match the
#' Ecopath production. Specifically it adjusts:
#' - `ks` to get the desired respiration,
#' - `w_repro_max` to get the desired gonadic proportion,
#' - `catchability` to get the desired yield, and
#' - `ext_mort` to get the desired mortality.
#'
#' Note that adjusting the `catchability` parameter to match the yield may no
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
#' @export
#' @md
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

matchProductionOnce <- function(params) {
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

    # Adjust respiration
    R <- getRespiration(params)
    Q <- getConsumption(params)
    # ecopath_production is only the somatic production. To get total production
    # we need to divide by 1 - gonad_proportion
    P_desired <- sp$ecopath_production / (1 - sp$gonad_proportion)
    R_desired <- sp$alpha * Q - P_desired
    if (any(R_desired < 0)) {
        stop("Negative respiration required.")
    }
    Rratio <- R_desired / R
    sp$ks <- sp$ks * Rratio
    species_params(params)$ks <- sp$ks

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

    # Return new steady state
    params |> steadySingleSpecies() |>
        matchBiomasses() |> steadySingleSpecies()
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
