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
    return(params)
}
