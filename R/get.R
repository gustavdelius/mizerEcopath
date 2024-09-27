#' Get somatic production
#'
#' For each species returns the rate at which somatic tissue is produced. This
#' is calculated as
#' \deqn{P_{s.i} = \int g_i(w) N_i(w) dw}
#' where \eqn{g(w)} is the somatic growth rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getEGrowth}})
#' and \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @return A named vector of somatic production for each species
#' @export
#' @family rate functions
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getSomaticProduction(params)
getSomaticProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    Ps <- as.vector((getEGrowth(params) * N) %*% dw)
    names(Ps) <- params@species_params$species
    return(Ps)
}

#' Get gonadic production
#'
#' For each species returns the rate at which gonads are produced.
#'
#' @param params A MizerParams object
#' @return A named vector of gonadic production for each species
#' @export
#' @family rate functions
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getGonadicProduction(params)
getGonadicProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    Pg <- as.vector((getERepro(params) * N) %*% dw)
    names(Pg) <- params@species_params$species
    return(Pg)
}

#' Get production
#'
#' For each species returns the rate at which biomass is produced.
#' This is calculated as
#' \deqn{P_i = \int E_{r.i}(w) N_i(w) dw}
#' where \eqn{E_{r.i}(w)} is the rate at which an individual of species \eqn{i}
#' and weight \eqn{w} allocates energy to reproduction and growth (calculated
#' with \code{\link{getEReproAndGrowth}}) and \eqn{N_i(w)} is the number density
#' of species \eqn{i} at weight \eqn{w}.
#'
#' The production rate is the sum of the gonadic and somatic production rates
#' obtained with \code{\link{getGonadicProduction}} and
#' \code{\link{getSomaticProduction}}.
#'
#' @param params A MizerParams object
#' @return A named vector of production for each species
#' @export
#' @family rate functions
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getProduction(params)
getProduction <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    P <- as.vector((getEReproAndGrowth(params) * N) %*% dw)
    names(P) <- params@species_params$species
    return(P)
}

#' Get consumption
#'
#' For each species returns the rate at which food is consumed. This is
#' calculated as
#' \deqn{Q_i = \int E_{e.i}(w) (1 - f_i(w)) N_i(w) dw}
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
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getConsumption(params)
getConsumption <- function(params, w_min = 0, w_max = Inf) {
    N <- initialN(params)
    q <- sweep(getEncounter(params) * (1 - getFeedingLevel(params)) * N,
               2, dw(params), "*")
    sel <- params@w >= w_min & params@w <= w_max
    Q <- rowSums(q[, sel])
    return(Q)
}

#' Get respiration
#'
#' For each species returns the rate at which energy is used for respiration.
#' This is calculated as
#' \deqn{R_i = \int k_i(w) N_i(w) dw}
#' where \eqn{k_i(w)} is the metabolic rate of an individual of species \eqn{i}
#' and weight \eqn{w} (calculated with \code{\link{getMetabolicRate}}) and
#' \eqn{N_i(w)} is the number density of species \eqn{i} at weight \eqn{w}.
#'
#' @param params A MizerParams object
#' @return A named vector of respiration for each species
#' @export
#' @family rate functions
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getRespiration(params)
getRespiration <- function(params) {
    N <- initialN(params)
    dw <- dw(params)
    R <- as.vector((getMetabolicRate(params) * N) %*% dw)
    names(R) <- params@species_params$species
    return(R)
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
#' @examples
#' params <- readRDS("models/Celtic_16_untuned.rds")
#' getUnassimilated(params)
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
#' \deqn{ZB_i = \int \mu_i(w) N_i(w) w dw}
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
#' \deqn{M0B = \int \mu_{ext.i}(w) N_i(w) w dw}
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
#' \deqn{M2B = \int \mu_{p.i}(w) N_i(w) w dw}
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

#' Get rate at which egg biomass is produced
#'
#' For each species the rate at which egg biomass is produced
#' is calculated as
#' \deqn{Eggs_i = R_i w_{0.i}}
#' where \eqn{R_i} is the rate of egg production for species \eqn{i}
#' and \eqn{w_{0.i}} is the egg weight of species \eqn{i}.
#'
#' @param params A MizerParams object
#' @return A named vector of egg biomass production for each species
#' @export
#' @family rate functions
getEggProduction <- function(params) {
    sp <- species_params(params)
    Egg_biomass <- getRDD(params) * sp$w_min
    return(Egg_biomass)
}

#' Get diet matrix
#'
#' Returns the diet matrix with one row for each predator species and one
#' column for each prey species and other ecosystem components. The entries
#' give the rate at which biomass flows from the prey species to the
#' predator species.
#'
#' @param params A MizerParams object
#' @param min_w_prey The minimum weight of prey species to include in the diet
#'   matrix
#' @param max_w_prey The maximum weight of prey species to include in the diet
#'   matrix
#' @return A matrix with dimnames `predator` and `prey`
#' @export
#' @family rate functions
getDietMatrix <- function(params, min_w_prey = 0, max_w_prey = Inf) {
    w_sel <- params@w >= min_w_prey & params@w <= max_w_prey
    N <- initialN(params)[, w_sel]
    dw <- dw(params)[w_sel]
    diet_matrix <- getDiet(params, proportion = FALSE)[, w_sel, ] |>
        sweep(c(1, 2), N, "*") |>
        sweep(2, dw, "*") |>
        aperm(c(1, 3, 2)) |>
        rowSums(dims = 2)
    return(diet_matrix)
}
