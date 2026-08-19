#' Get somatic production for each species
#'
#' For each species, returns the rate at which somatic biomass is produced by all
#' individuals of that species.
#'
#' The somatic production rate \eqn{P_{s.i}} is calculated as
#' \deqn{P_{s.i} = \int_{w_0}^{w_{max}} g_i(w) N_i(w) dw}
#' where \eqn{g(w)} is the somatic growth rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getEGrowth}})
#' and \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' @section Discretisation:
#' The integral is evaluated on mizer's finite-volume size grid as
#' \eqn{\sum_j N_{ij}\bar K_{ij}\Delta w_j}, where \eqn{K_i(w)} is everything in
#' the integrand except the number density \eqn{N_i(w)}, including any
#' size-window mask coming from the `min_w`/`max_w` arguments. By default
#' \eqn{\bar K_{ij}} is the value of \eqn{K_i} at the left edge of the bin,
#' which is what mizer's own summary functions use. If the model has bin
#' averaging switched on with [mizer::second_order_w()], then \eqn{\bar K_{ij}}
#' is instead the trapezoidal average of \eqn{K_i} over the bin, making the
#' result second order in the bin width and consistent with the quadratures
#' mizer uses in that case.
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of somatic production for each species
#' @export
#' @family rate functions
#' @examples
#' getSomaticProduction(NS_params)
getSomaticProduction <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(getEGrowth(params) * sel, params)
    Ps <- as.vector((N * K) %*% dw(params))
    names(Ps) <- params@species_params$species
    return(Ps)
}

#' Get gonadic production for each species
#'
#' For each species, returns the rate at which gonad biomass is produced by all
#' individuals of that species.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of gonadic production for each species
#' @export
#' @family rate functions
#' @examples
#' getGonadicProduction(NS_params)
getGonadicProduction <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(getERepro(params) * sel, params)
    Pg <- as.vector((N * K) %*% dw(params))
    names(Pg) <- params@species_params$species
    return(Pg)
}


#' Get total production rate for each species
#'
#' For each species, returns the rate at which biomass is produced by all
#' individuals of that species. This includes both somatic and gonadic
#' production, which can be obtained individually with [getSomaticProduction()]
#' and [getGonadicProduction()]. This differs from the Ecopath definition of
#' production, see [getProduction()].
#'
#' Used by [matchConsumption()] to calculate how much assimilated energy is left
#' over after production, in order to estimate the required metabolic loss.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of gonadic production for each species
#' @export
#' @family rate functions
#' @examples
#' getTotalProduction(NS_params)
getTotalProduction <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(getEReproAndGrowth(params) * sel, params)
    Pg <- as.vector((N * K) %*% dw(params))
    names(Pg) <- params@species_params$species
    return(Pg)
}

#' Get production rate as defined by Ecopath
#'
#' Unlike [getTotalProduction()], this function follows the Ecopath definition
#' whereby the production rate in a
#' mass-balanced model is equal to the rate at which biomass is removed by
#' mortality. Thus the production rate includes the somatic production rate
#' that is due to growth of individuals and it includes the rate at which
#' offspring biomass is produced during reproduction. It does not include the
#' full rate of gonad production because most of that is lost during the
#' inefficient reproduction process.
#' The Ecopath production rate is the sum of the somatic production rate
#' obtained with \code{\link{getSomaticProduction}} and the rate of production
#' of offspring biomass obtained with \code{\link{getOffspringProduction}}.
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of production for each species
#' @export
#' @family rate functions
#' @examples
#' getProduction(NS_params)
getProduction <- function(params, ...) {
    getSomaticProduction(params, ...) + getOffspringProduction(params)
}

#' Get consumption rate for each species
#'
#' For each species, this function returns the rate at which food is consumed
#' by all individuals of that species. This corresponds to the consumption rate
#' \eqn{Q} used by Ecopath models.
#'
#' The consumption is calculated as
#' \deqn{Q_i = \int_{w_0}^{w_{max}} E_{e.i}(w) (1 - f_i(w)) N_i(w) dw}
#' where
#' \eqn{E_{e.i}(w)} is the encounter rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getEncounter}}), \eqn{f_i(w)}
#' is the feeding level of an individual of species \eqn{i} and weight \eqn{w}
#' (calculated with \code{\link{getFeedingLevel}}) and \eqn{N_i(w)} is the
#' number density of species \eqn{i} at weight \eqn{w}.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritParams getDietMatrix
#' @return A named vector of consumption rate for each species
#' @export
#' @family rate functions
#' @examples
#' getConsumption(NS_params)
getConsumption <- function(params, min_w_pred = 0, max_w_pred = Inf) {
    N <- initialN(params)
    sel <- params@w >= min_w_pred & params@w <= max_w_pred
    K <- sweep(getEncounter(params) * (1 - getFeedingLevel(params)),
               2, sel, "*")
    K <- bin_average_weight(K, params)
    Q <- drop((N * K) %*% dw(params))
    return(Q)
}

#' Get metabolic respiration for each species
#'
#' For each species, returns the rate at which energy is lost to metabolic respiration
#' by all individuals of that species. Note that this is different from the
#' respiration rate \eqn{R} used by Ecopath models, see [getRespiration()].
#'
#' The respiration rate \eqn{K_i} is calculated as
#' \deqn{K_i = \int_{w_0}^{w_{max}} k_i(w) N_i(w) dw}
#' where \eqn{k_i(w)} is the metabolic rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getMetabolicRate}}) and
#' \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' Used internally by [matchConsumption()] to help scale metabolic loss to match observed consumption.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of respiration for each species
#' @export
#' @family rate functions
#' @examples
#' getMetabolicRespiration(NS_params)
getMetabolicRespiration <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(getMetabolicRate(params) * sel, params)
    R <- as.vector((N * K) %*% dw(params))
    names(R) <- params@species_params$species
    return(R)
}

#' Get rate of biomass loss due to reproduction
#'
#' For each species, returns the difference between the rate of gonad production
#' (calculated with [getGonadicProduction]) and the rate of offspring biomass
#' production (calculated with [getOffspringProduction]).
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of reproduction biomass loss for each species
#' @export
#' @family rate functions
#' @examples
#' getReproductiveLoss(NS_params)
getReproductiveLoss <- function(params, ...) {
    getGonadicProduction(params, ...) - getOffspringProduction(params)
}

#' Get respiration rate as defined by Ecopath
#'
#' Unlike [getMetabolicRespiration()], this function follows the Ecopath
#' definition of respiration as the rate at which assimilated food is NOT used
#' for production. This is expressed by the first master equation of Ecopath:
#' \deqn{Q = P + R + U} where \eqn{Q} is the consumption rate, \eqn{P} is the
#' production rate, \eqn{R} is the respiration rate, and \eqn{U} is the
#' unassimilated part of the consumption rate, i.e., \eqn{U = (1 - \alpha) Q},
#' where \eqn{\alpha} is the assimilation efficiency. Solving this for \eqn{R}
#' we get: \deqn{R = \alpha Q - P}
#'
#' In mizer, this respiration rate has two components: the metabolic respiration
#' and the loss due to gonad production that does not result in offspring
#' biomass. Thus the respiration rate is the sum of the metabolic respiration
#' rate (calculated with [getMetabolicRespiration()]) and the rate of
#' biomass loss due to reproduction (calculated with
#' [getReproductiveLoss()]).
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of respiration for each species
#' @export
#' @family rate functions
#' @examples
#' getRespiration(NS_params)
getRespiration <- function(params, ...) {
    getMetabolicRespiration(params, ...) + getReproductiveLoss(params, ...)
}

#' Get unassimilated food
#'
#' For each species, returns the rate at which food is not assimilated.
#' This is calculated as
#' \deqn{U_i = (1 - \alpha_i) Q_i}
#' where \eqn{\alpha_i} is the assimilation efficiency of species \eqn{i}
#' and \eqn{Q_i} is the consumption rate of species \eqn{i} (calculated with
#' \code{\link{getConsumption}}).
#'
#' @param params A MizerParams object
#' @return A named vector of unassimilated food for each species
#' @export
#' @family rate functions
#' @examples
#' getUnassimilated(NS_params)
getUnassimilated <- function(params) {
    sp <- species_params(params)
    Q <- getConsumption(params)
    U <- (1 - sp$alpha) * Q
    return(U)
}

#' Get rate at which biomass is lost due to mortality
#'
#' For each species, returns the rate at which biomass is lost due to mortality.
#'
#' This is calculated as
#' \deqn{ZB_i = \int_{w_0}^{w_{max}} \mu_i(w) N_i(w) w dw}
#' where \eqn{\mu_i(w)} is the mortality rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getMort}}) and \eqn{N_i(w)} is the
#' number density of species \eqn{i} at weight \eqn{w}.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of biomass loss rate due to mortality for each species
#' @export
#' @family rate functions
#' @examples
#' getZB(NS_params)
getZB <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(sweep(getMort(params) * sel, 2, w(params), "*"),
                          params)
    ZB <- as.vector((N * K) %*% dw(params))
    names(ZB) <- params@species_params$species
    return(ZB)
}

#' Get rate at which biomass is lost due to external mortality
#'
#' For each species, returns the rate at which biomass is lost due to external
#' mortality.
#'
#' This is calculated as
#' \deqn{M0B_i = \int_{w_0}^{w_{max}} \mu_{ext.i}(w) N_i(w) w dw}
#' where \eqn{\mu_{ext.i}(w)} is the external mortality rate of an individual
#' of species \eqn{i} and weight \eqn{w} (obtained with \code{\link{getExtMort}})
#' and \eqn{N_i(w)} is the number density of species i at weight w.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of biomass loss rate due to external mortality for
#'   each species
#' @export
#' @family rate functions
#' @examples
#' getM0B(NS_params)
getM0B <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(sweep(getExtMort(params) * sel, 2, w(params), "*"),
                          params)
    M0B <- as.vector((N * K) %*% dw(params))
    names(M0B) <- params@species_params$species
    return(M0B)
}

#' Get rate at which biomass is lost due to predation mortality
#'
#' For each species, returns the rate at which biomass is lost due to predation
#' mortality.
#'
#' This is calculated as
#' \deqn{M2B_i = \int_{w_0}^{w_{max}} \mu_{p.i}(w) N_i(w) w dw}
#' where \eqn{\mu_{p.i}(w)} is the predation mortality rate of an individual
#' of species \eqn{i} and weight \eqn{w} (obtained with \code{\link{getPredMort}})
#' and \eqn{N_i(w)} is the number density of species i at weight w.
#'
#' @inheritSection getSomaticProduction Discretisation
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of biomass loss rate due to predation mortality for
#'   each species
#' @export
#' @family rate functions
#' @examples
#' getM2B(NS_params)
getM2B <- function(params, ...) {
    N <- initialN(params)
    sel <- get_size_range_array(params, ...)
    K <- bin_average_weight(sweep(getPredMort(params) * sel, 2, w(params), "*"),
                          params)
    M2B <- as.vector((N * K) %*% dw(params))
    names(M2B) <- params@species_params$species
    return(M2B)
}

#' Get Ecotrophic Efficiency for each species
#'
#' For each species, returns the ecotrophic efficiency, the proportion of
#' production that is not lost to mortality external to the model:
#' \deqn{EE_i = 1-\frac{M0B_i}{P_i}}
#' where \eqn{M0B_i} is the rate of biomass loss due to external mortality
#' and \eqn{P_i} is the rate of production.
#'
#' @param params A MizerParams object
#' @inheritDotParams mizer::get_size_range_array -params
#' @return A named vector of ecotrophic efficiency for each species
#' @export
#' @family rate functions
#' @examples
#' getEcotrophicEfficiency(NS_params)
getEcotrophicEfficiency <- function(params, ...) {
    P <- getProduction(params, ...)
    M0B <- getM0B(params, ...)
    EE <- 1 - M0B / P
    return(EE)
}

#' Get rate at which offspring biomass is produced
#'
#' For each species, returns the rate at which offspring biomass is produced is
#' the product of the rate of offspring production (obtained with [mizer::getRDD()])
#' and the offspring biomass \eqn{w_{0.i}}.
#'
#' @param params A MizerParams object
#' @return A named vector of offspring biomass production for each species
#' @export
#' @family rate functions
#' @examples
#' getOffspringProduction(NS_params)
getOffspringProduction <- function(params) {
    sp <- species_params(params)
    offspring_biomass <- getRDD(params) * sp$w_min
    return(offspring_biomass)
}

#' Get the efficiency of reproduction
#'
#' For each species, returns the fraction of the energy invested into
#' reproduction that arrives in the size spectrum as offspring biomass. For a
#' `MizerParams` object this is calculated from the initial state, for a
#' `MizerSim` object it is returned for every saved time.
#'
#' Mature individuals invest energy into reproduction at the total rate
#' \eqn{E_{R.i}} (grams per year), obtained with [getGonadicProduction()].
#' Mizer assumes that half of that investment is made by females, so
#' \eqn{E_{R.i} / 2} is what is available for making eggs. Offspring biomass is
#' produced at the rate \eqn{R_{dd.i} w_{0.i}}, obtained with
#' [getOffspringProduction()]. The reproductive efficiency is the ratio of the
#' two:
#' \deqn{\frac{R_{dd.i}\, w_{0.i}}{E_{R.i} / 2} =
#'       \epsilon_i \, \frac{R_{dd.i}}{R_{di.i}},}
#' where \eqn{\epsilon_i} is the `erepro` species parameter and
#' \eqn{R_{di.i}} is the density-independent rate of egg production, see
#' [mizer::getRDI()]. The second factor \eqn{R_{dd.i}/R_{di.i}} is the fraction
#' of eggs that survives density dependence and equals
#' `1 - reproduction_level`, see [mizer::setBevertonHolt()].
#'
#' The function evaluates the right-hand side of the above equation, because
#' that is the form mizer itself uses when it projects reproduction. The two
#' sides agree exactly, whatever [mizer::second_order_w()] is set to, because
#' [getGonadicProduction()] discretises the reproduction integral with the same
#' quadrature that [mizer::getRDI()] uses.
#'
#' A value above 1 is energetically impossible: the model would be producing
#' more offspring biomass than the adults invested. Use
#' [plotReproductiveEfficiency()] to see how far each species is from that
#' limit and how much of the loss is due to density dependence. Species that
#' do not reproduce at all (\eqn{R_{di.i} = 0}) get `NA`.
#'
#' @param object A MizerParams or MizerSim object
#' @return If called with a MizerParams object, a named vector with the
#'   reproductive efficiency of each species in the model. If called with a
#'   MizerSim object, an `ArrayTimeBySpecies` object (time x species) with the
#'   reproductive efficiency at each saved time. The values are dimensionless
#'   fractions.
#' @export
#' @family rate functions
#' @seealso [plotReproductiveEfficiency()], [getOffspringProduction()],
#'   [getGonadicProduction()]
#' @examples
#' getReproductiveEfficiency(celtic_params)
getReproductiveEfficiency <- function(object) {
    UseMethod("getReproductiveEfficiency")
}

#' @export
getReproductiveEfficiency.MizerParams <- function(object) {
    params <- validParams(object)
    rdi <- getRDI(params)
    rdd <- getRDD(params)
    efficiency <- params@species_params$erepro * ifelse(rdi > 0, rdd / rdi,
                                                        NA_real_)
    names(efficiency) <- params@species_params$species
    efficiency
}

#' @export
getReproductiveEfficiency.MizerSim <- function(object) {
    sim <- object
    params <- sim@params
    rdi <- getRDI(sim)
    rdd <- getRDD(sim)
    efficiency <- matrix(ifelse(rdi > 0, rdd / rdi, NA_real_),
                         nrow = nrow(rdi), dimnames = dimnames(rdi))
    efficiency <- sweep(efficiency, 2, params@species_params$erepro, "*")
    ArrayTimeBySpecies(efficiency, value_name = "Reproductive efficiency",
                       params = params)
}
