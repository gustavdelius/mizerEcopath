#' Match the respiration rate to that required by the ecopath production rate.
#'
#' This function sets the respiration rate (called the metabolic rate in mizer)
#' according to the first master equation of Ecopath, which is:
#' \deqn{Q = P + R + U}
#' where \eqn{Q} is the consumption rate, \eqn{P} is the production rate,
#' \eqn{R} is the respiration rate, and \eqn{U} is the unassimilated part of the
#' consumption rate, i.e., \eqn{U = (1 - \alpha) Q}, where \eqn{\alpha} is the
#' assimilation efficiency. Solving this for \eqn{R} we get:
#' \deqn{R = Q \alpha - P}
#'
#' The production rate of a species is the sum of the somatic production rate and
#' the gonad production rate. The somatic production rate is the production rate
#'
#' The Ecopath production rate is the somatic production rate, i.e., the production
#' rate of the somatic tissue. To get the total production rate we need to divide
#' the somatic production rate by \eqn{1 - f}, where \eqn{f} is the gonad
#' proportion. This is because the gonad production is not included in the Ecopath
#' production rate.
#'
#' This function calculates \eqn{Q} using `getConsumption()`
#'
#' @param params A MizerParams object
#' @param steady Whether to return the model to a steady state after adjusting
#'  the respiration. Default is TRUE.
#' @return A MizerParams object with the respiration matched
#' @family match functions
#' @export
matchRespirationOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }

    R <- getRespiration(params)
    Q <- getConsumption(params)
    # ecopath_production is only the somatic production.
    P_desired <- sp$ecopath_production + getGonadicProduction(params)
    R_desired <- sp$alpha * Q - P_desired
    if (any(R_desired < 0)) {
        warning("Negative respiration required.")
    }
    Rfactor <- R_desired / R
    # If the metabolic rate was set manually, we need to rescale it, otherwise we
    # rescale the species parameters used to calculate the metabolic rate
    if (!is.null(comment(metab(params)))) {
        warning("This function has rescaled the metabolic rate that was set manually.")
        params@metab[] <- params@metab[] * Rfactor
    } else {
        species_params(params)$ks <- sp$ks * Rfactor
        if (hasName(sp, "k")) {
            species_params(params)$k <- sp$k * Rfactor
        }
    }

    if (steady) {
        # Determine new steady state
        params <- params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    params@time_modified <- lubridate::now()
    return(params)
}
