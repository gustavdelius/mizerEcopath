#' Set feeding levels while maintaining constant growth and reproduction
#'
#' @description
#' This function takes a tuned, non-interacting allometric model (where initially
#' \eqn{h = \infty} and \eqn{k_s = 0}) and determines the new maximum intake rate
#' \eqn{h} and metabolic loss \eqn{k_s} parameters. It does this by adjusting
#' metabolic loss, encounter rate, and maximum intake rate so that the energy
#' available for reproduction and growth remains unchanged, despite changes in
#' the feeding level and critical feeding level.
#'
#' @details
#' In order to keep the energy for growth and reproduction \eqn{E_r} unchanged,
#' the function recalculates the encounter rate \eqn{E}, metabolic loss rate
#' \eqn{k} and maximum intake rate \eqn{h} as follows:
#'
#' \deqn{E=E_r\frac{f}{\alpha(f-f_c)(1-f)}}
#' \deqn{h = E_r\frac{1}{\alpha(f-f_c)}}
#' \deqn{k_s = E_r\frac{f_c}{f-f_c}}
#'
#' This methodology allows you to set a desired feeding level for each species
#' without altering the underlying energy balance established during the
#' initial tuning phase.
#'
#' @section Mandatory Prerequisites:
#' This function will throw an error unless the following conditions are met:
#' \itemize{
#'   \item **Allometric Scaling**: The model must use allometric rates where the
#'         intake exponent (\eqn{n}) and metabolism exponent (\eqn{p}) are equal.
#'         This ensures the balance remains consistent across all fish sizes.
#'   \item **Initial State**: The species parameters \code{h} and \code{ks}
#'         must have the values \code{h = Inf} and \code{ks = 0}.
#'   \item **Non-interacting Model**: All encounter must be external encounter.
#'         This function is not compatible with models where species interact with
#'         a dynamically realized resource. If interactions exist, call
#'         \code{makeNoninteracting()} first.
#' }
#' @param params A \linkS4class{MizerParams} object.
#' @param f The target feeding level \eqn{f} for each
#'   species. Must be a vector with one value for each species or a single
#'   value used for all. Defaults to species parameter `f0` if exists or else
#'   to 0.6.
#' @param f_c The feeding level \eqn{f_c} at which
#'   growth/repro is zero. Defaults to species parameter `fc` if exists or else
#'   to 0.2.
#'
#' @return A \linkS4class{MizerParams} object with updated \eqn{h, k_s},
#'   and \code{ext_encounter} rates.
#' @export
setFeedingLevels <- function(params, f, f_c) {
    sp <- params@species_params

    # Check that params describes a non-interacting model
    if (!isTRUE(all.equal(getEncounter(params), getExtEncounter(params)))) {
        stop("This function only works for models where all encounter is external encounter. Try calling `makeNoninteracting()` first.")
    }

    # Check that n and p exponents are equal
    bad <- which(sp$n != sp$p)
    if (length(bad) > 0) {
        stop(paste(
            "Exponents n and p must be equal for species:",
            paste(rownames(sp)[bad], collapse = ", ")
        ))
    }

    # Check that h is Inf
    if ("h" %in% names(sp) && !all(is.infinite(sp$h))) {
        stop("h must be Inf before calling this function")
    }

    #Check that ks is 0
    if (!all(sp$ks == 0)) {
        stop("ks must be 0 before calling this function")
    }

    # Check that feeding level is supplied
    if (missing(f)) {
        sp <- set_species_param_default(sp, "f0", 0.6)
        f <- sp$f0
    }

    if (missing(f_c)) {
        sp <- set_species_param_default(sp, "fc", 0.2)
        f_c <- sp$fc
    }

    assert_that(is.numeric(f))
    assert_that(is.numeric(f_c))

    # If f is a single value, make it a vector
    if (length(f) == 1) {
        f <- rep(f, nrow(sp))
    }

    # If f_c is a single value, make it a vector
    if (length(f_c) == 1) {
        f_c <- rep(f_c, nrow(sp))
    }

    # Check that f is a vector with one value for each species
    if (length(f) != nrow(sp)) {
        stop(paste("The length of f vector (", length(f), ") does not match the number of species (", nrow(sp), ")."))
    }

    # Check that f_c is a vector with one value for each species
    if (length(f_c) != nrow(sp)) {
        stop(paste("The length of f_c vector (", length(f_c), ") does not match the number of species (", nrow(sp), ")."))
    }

    # Check that the feeding level is in the correct range
    if (any(f < 0 | f >= 1)) {
        stop("Feeding level must be positive and strictly less than 1.")
    }

    # Check that the critical feeding level is in the correct range
    if (any(f_c < 0 | f_c >= 1)) {
        stop("Critical feeding level must be positive and strictly less than 1.")
    }

    # Ensure critical feeding level is less than feeding level
    bad <- which(f_c >= f)
    if (length(bad) > 0) {
        stop(paste(
            "Critical feeding level must be less than the feeding level for species:",
            paste(rownames(sp)[bad], collapse = ", ")
        ))
    }

    # Check that the params object is allometric
    if (!isTRUE(isAllometric(params))) {
        stop("This function only works for models made up of allometric rates.")
    }

    # Save feeding level and f_c in species params
    # Don't use `given_species_params<-` because we do not want to trigger
    # a recalculation of anything
    params@given_species_params$f0 <- f
    sp$f0 <- f
    params@given_species_params$fc <- f_c
    sp$fc <- f_c

    alpha <- params@species_params$alpha
    n <- params@species_params$n
    E_r <- getEReproAndGrowth(params)
    # Divide out power law to get coefficient
    E_r_0 <- E_r / exp(outer(n, log(params@w)))
    # Check that we had a power law
    if (any(abs(E_r_0 - E_r_0[, 1]) > 1e-13)) {
        stop("The energy for growth and reproduction needs to be a power law with exponent ", n)
    }
    E_r_0 <- E_r_0[, 1]

    # Set New External Encounter
    ext_encounter(params) <- E_r * f / (alpha * (f - f_c) * (1 - f))

    # Set new species params
    sp$ks <- E_r_0 * f_c / (f - f_c)
    sp$h <- E_r_0 / (alpha * (f - f_c))

    # Update the params object
    params@species_params <- sp

    # Comment out so that internal calculations will be triggered in the next step
    comment(params@intake_max) <- NULL
    comment(params@metab)      <- NULL
    # Note do not comment out any other internal calculations, as these should
    # not be changed by the function.

    # Rebuild of the internal arrays based on the new sp constants
    params <- setMaxIntakeRate(params)
    params <- setMetabolicRate(params)

    return(params)
}
