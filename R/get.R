#' Get somatic production
#'
#' For each species returns the rate at which somatic biomass is produced. This
#' is calculated as
#' \deqn{P_{s.i} = \int_{w_0}^{w_{max}} g_i(w) N_i(w) dw}
#' where \eqn{g(w)} is the somatic growth rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getEGrowth}})
#' and \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @return A named vector of somatic production for each species
#' @export
#' @family rate functions
getSomaticProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    Ps <- as.vector((getEGrowth(params) * N) %*% dw)
    names(Ps) <- params@species_params$species
    return(Ps)
}

#' Get gonadic production
#'
#' For each species returns the rate at which gonad biomass is produced.
#'
#' @param params A MizerParams object
#' @return A named vector of gonadic production for each species
#' @export
#' @family rate functions
getGonadicProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    Pg <- as.vector((getERepro(params) * N) %*% dw)
    names(Pg) <- params@species_params$species
    return(Pg)
}


#' Get total production
#'
#' For each species returns the rate at which biomass is produced. This includes
#' both somatic and gonadic production.
#'
#' @param params A MizerParams object
#' @return A named vector of gonadic production for each species
#' @export
#' @family rate functions
getTotalProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    Pg <- as.vector((getEReproAndGrowth(params) * N) %*% dw)
    names(Pg) <- params@species_params$species
    return(Pg)
}

#' Get production
#'
#' For each species returns the rate at which biomass is produced.
#'
#' We follow the Ecopath definition whereby the production rate in a
#' mass-balanced model is equal to the rate at which biomass is removed by
#' mortality. Thus the production rate includes the somatic production rate
#' that is due to growth of individuals and it includes the rate at which
#' offspring biomass is produced during reproduction. It does not include the
#' full rate of gonad production because most of that is lost during the
#' inefficient reproduction process.
#'
#' The production rate is the sum of the somatic production rate
#' obtained with \code{\link{getSomaticProduction}} and the rate of production
#' of offspring biomass obtained with \code{\link{getOffspringProduction}}.
#'
#' @param params A MizerParams object
#' @return A named vector of production for each species
#' @export
#' @family rate functions
getProduction <- function(params) {
    getSomaticProduction(params) + getOffspringProduction(params)
}

#' Get consumption
#'
#' For each species, this function returns the rate at which food is consumed.
#' This is calculated as
#' \deqn{Q_i = \int_{w_0}^{w_{max}} E_{e.i}(w) (1 - f_i(w)) N_i(w) dw}
#' where
#' \eqn{E_{e.i}(w)} is the encounter rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getEncounter}}), \eqn{f_i(w)}
#' is the feeding level of an individual of species \eqn{i} and weight \eqn{w}
#' (calculated with \code{\link{getFeedingLevel}}) and \eqn{N_i(w)} is the
#' number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @param w_min The minimum weight of prey species to include in the consumption
#'  rate calculation
#' @param w_max The maximum weight of prey species to include in the consumption
#'  rate calculation
#' @return A named vector of consumption rate for each species
#' @export
#' @family rate functions
getConsumption <- function(params, w_min = 0, w_max = Inf) {
    N <- initialN(params)
    q <- sweep(getEncounter(params) * (1 - getFeedingLevel(params)) * N,
               2, dw(params), "*")
    sel <- params@w >= w_min & params@w <= w_max
    Q <- rowSums(q[, sel, drop = FALSE])
    return(Q)
}

#' Get metabolic respiration
#'
#' For each species returns the rate at which energy is used for respiration.
#' This is calculated as
#' \deqn{K_i = \int_{w_0}^{w_{max}} k_i(w) N_i(w) dw}
#' where \eqn{k_i(w)} is the metabolic rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getMetabolicRate}}) and
#' \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @return A named vector of respiration for each species
#' @export
#' @family rate functions
getMetabolicRespiration <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    R <- as.vector((getMetabolicRate(params) * N) %*% dw)
    names(R) <- params@species_params$species
    return(R)
}

#' Get rate of biomass loss due to reproduction
#'
#' Returns the difference between the rate of gonad production (calculated with
#' `getGonadicProduction`) and the rate of offspring biomass production (calculated with
#' `getOffspringProduction`).
#'
#' @param params A MizerParams object
#' @return A named vector of reproduction biomass loss for each species
#' @export
#' @family rate functions
getReproductiveLoss <- function(params) {
    getGonadicProduction(params) - getOffspringProduction(params)
}

#' Get rate of metabolic respiration
#'
#' We use the Ecopath definition of respiration as the rate at which assimilated
#' food is NOT used for production. This is expressed by the first master
#' equation of Ecopath: \deqn{Q = P + R + U} where \eqn{Q} is the consumption
#' rate, \eqn{P} is the production rate, \eqn{R} is the respiration rate,
#' and \eqn{U} is the unassimilated part of the consumption rate, i.e.,
#' \eqn{U = (1 - \alpha) Q}, where \eqn{\alpha} is the assimilation efficiency.
#' Solving this for \eqn{R} we get: \deqn{R = \alpha Q - P}
#'
#' In mizer, this respiration rate has two components: the metabolic respiration
#' and the loss due to gonad production that does not result in offspring
#' biomass. Thus the respiration rate is the sum of the metabolic respiration
#' rate (calculated with `getMetabolicRespiration()`) and the rate of
#' biomass loss due to reproduction (calculated with
#' `getReproductiveLoss()`).
#'
#' @param params A MizerParams object
#' @return A named vector of respiration for each species
#' @export
#' @family rate functions
getRespiration <- function(params) {
    getMetabolicRespiration(params) + getReproductiveLoss(params)
}

#' Get unassimilated food
#'
#' For each species returns the rate at which food is not assimilated.
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
getUnassimilated <- function(params) {
    sp <- species_params(params)
    Q <- getConsumption(params)
    U = (1 - sp$alpha) * Q
    return(U)
}

#' Get rate at which biomass is lost due to mortality
#'
#' For each species returns the rate at which biomass is lost due to mortality
#' is calculated as
#' \deqn{ZB_i = \int_{w_0}^{w_{max}} \mu_i(w) N_i(w) w dw}
#' where \eqn{\mu_i(w)} is the mortality rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getMort}}) and \eqn{N_i(w)} is the
#' number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @return A named vector of biomass loss rate due to mortality for each species
#' @export
#' @family rate functions
getZB <- function(params) {
    N <- initialN(params)
    w <- w(params)
    dw <- dw(params)
    ZB <- as.vector((getMort(params) * N) %*% (w * dw))
    names(ZB) <- params@species_params$species
    return(ZB)
}

#' Get rate at which biomass is lost due to external mortality
#'
#' For each species returns the rate at which biomass is lost due to external
#' mortality is calculated as
#' \deqn{M0B_i = \int_{w_0}^{w_{max}} \mu_{ext.i}(w) N_i(w) w dw}
#' where \eqn{\mu_{ext.i}(w)} is the external mortality rate of an individual
#' of species \eqn{i} and weight \eqn{w} (obtained with \code{\link{getExtMort}})
#' and \eqn{N_i(w)} is the number density of species i at weight w.
#'
#' @param params A MizerParams object
#' @return A named vector of biomass loss rate due to external mortality for
#'   each species
#' @export
#' @family rate functions
getM0B <- function(params) {
    N <- initialN(params)
    w <- w(params)
    dw <- dw(params)
    M0B <- as.vector((getExtMort(params) * N) %*% (w * dw))
    names(M0B) <- params@species_params$species
    return(M0B)
}

#' Get rate at which biomass is lost due to predation mortality
#'
#' For each species returns the rate at which biomass is lost due to predation
#' mortality is calculated as
#' \deqn{M2B_i = \int_{w_0}^{w_{max}} \mu_{p.i}(w) N_i(w) w dw}
#' where \eqn{\mu_{p.i}(w)} is the predation mortality rate of an individual
#' of species \eqn{i} and weight \eqn{w} (obtained with \code{\link{getPredMort}})
#' and \eqn{N_i(w)} is the number density of species i at weight w.
#'
#' @param params A MizerParams object
#' @return A named vector of biomass loss rate due to predation mortality for
#'   each species
#' @export
#' @family rate functions
getM2B <- function(params) {
    N <- initialN(params)
    w <- w(params)
    dw <- dw(params)
    M2B <- as.vector((getPredMort(params) * N) %*% (w * dw))
    names(M2B) <- params@species_params$species
    return(M2B)
}

#' Get Ecotrophic Efficiency
#'
#' For each species returns the ecotrophic efficiency, the proportion of
#' production that is not lost to mortality external to the model:
#' \deqn{EE_i = 1-\frac{M0B_i}{P_i}}
#' where \eqn{M0B_i} is the rate of biomass loss due to external mortality
#' and \eqn{P_i} is the rate of production.
#'
#' @param params A MizerParams object
#' @return A named vector of ecotrophic efficiency for each species
#' @export
#' @family rate functions
getEcotrophicEfficiency <- function(params) {
    P <- getProduction(params)
    M0B <- getM0B(params)
    EE <- 1 - M0B / P
    return(EE)
}

#' Get rate at which offspring biomass is produced
#'
#' For each species the rate at which offspring biomass is produced is the
#' product of the rate of offspring production and the offspring biomass
#' \eqn{w_{0.i}}.
#'
#' @param params A MizerParams object
#' @return A named vector of egg biomass production for each species
#' @export
#' @family rate functions
getOffspringProduction <- function(params) {
    sp <- species_params(params)
    Offspring_biomass <- getRDD(params) * sp$w_min
    return(Offspring_biomass)
}
