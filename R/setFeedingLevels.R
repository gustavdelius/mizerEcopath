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
#' The function follows a three-step methodology to ensure that the consumption
#' rate does not change in spite of the change of feeding level:
#'
#' 1. Metabolic Loss Calculation: It uses the critical feeding level \eqn{f_c}
#' (the level where energy for growth/repro is zero) and the old encounter rate
#' \eqn{E^{old}} to find the new metabolic loss \eqn{k_{i,w}}:
#' \deqn{k_{i}(w) = \alpha_i (1 - f_{c,i}) E_{i}^{old}(w)}
#'
#' 2. Encounter Rate Adjustment: It calculates a new encounter rate (\eqn{E_{new}})
#' required to maintain the original energy available for growth and reproduction
#' (\eqn{E_{r,i}(w)}) at the target feeding level \eqn{f}:
#' \deqn{E_{i}^{new}(w) = \frac{E_{r,i}(w) + k_{i}(w)}{\alpha_i (1 - f_i)}}
#'
#' 3. Maximum Intake Back-calculation: Finally, the maximum intake rate \eqn{h}
#' is back-calculated from the new encounter rate to achieve the desired feeding level:
#' \deqn{h_{i}(w) = \frac{E_{i}^{new}(w) (1 - f_i)}{f_i}}
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
#' @param feeding_level The target feeding level \eqn{f} for each
#'   species. Must be a vector with one value for each species or a single
#'   value used for all. Defaults to species parameter `f0` if exists or else
#'   to 0.6.
#' @param critical_feeding_level The feeding level \eqn{f_c} at which
#'   growth/repro is zero. Defaults to species parameter `fc` if exists or else
#'   to 0.2.
#'
#' @return A \linkS4class{MizerParams} object with updated \eqn{h, k_s},
#'   and \code{ext_encounter} rates.
#' @export
setFeedingLevels <- function(params, feeding_level, critical_feeding_level) {
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
    if (missing(feeding_level)) {
        sp <- set_species_param_default(sp, "f0", 0.6)
        feeding_level <- sp$f0
    }

    if (missing(critical_feeding_level)) {
        sp <- set_species_param_default(sp, "fc", 0.2)
        critical_feeding_level <- sp$fc
    }

    assert_that(is.numeric(feeding_level))
    assert_that(is.numeric(critical_feeding_level))

    # If feeding_level is a single value, make it a vector
    if (length(feeding_level) == 1) {
        feeding_level <- rep(feeding_level, nrow(sp))
    }

    # If critical_feeding_level is a single value, make it a vector
    if (length(critical_feeding_level) == 1) {
        critical_feeding_level <- rep(critical_feeding_level, nrow(sp))
    }

    # Check that feeding_level is a vector with one value for each species
    if (length(feeding_level) != nrow(sp)) {
        stop(paste("The length of feeding_level vector (", length(feeding_level), ") does not match the number of species (", nrow(sp), ")."))
    }

    # Check that critical_feeding_level is a vector with one value for each species
    if (length(critical_feeding_level) != nrow(sp)) {
        stop(paste("The length of critical_feeding_level vector (", length(critical_feeding_level), ") does not match the number of species (", nrow(sp), ")."))
    }

    # Check that the feeding level is in the correct range
    if (any(feeding_level < 0 | feeding_level >= 1)) {
        stop("Feeding level must be positive and strictly less than 1.")
    }

    # Check that the critical feeding level is in the correct range
    if (any(critical_feeding_level < 0 | critical_feeding_level >= 1)) {
        stop("Critical feeding level must be positive and strictly less than 1.")
    }

    # Ensure critical feeding level is less than feeding level
    bad <- which(critical_feeding_level >= feeding_level)
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

    # Save feeding level and critical_feeding_level in species params
    # Don't use `given_species_params<-` because we do not want to trigger
    # a recalculation of anything
    params@given_species_params$f0 <- feeding_level
    params@species_params$f0 <- feeding_level
    params@given_species_params$fc <- critical_feeding_level
    params@species_params$fc <- critical_feeding_level

    w <- params@w

    Eriw <- getEReproAndGrowth(params)
    Eiw <- getEncounter(params)

    # Ensure that Eiw*alpha=Eriw
    for (i in seq(nrow(sp))) {
        if (any(Eiw[i,] * sp$alpha[i] != Eriw[i,])) {
            stop(paste("Eiw * alpha does not equal Eriw"))
        }
    }

    #Save 2 matrices of Eiw values one which will be over written with
    #new values and one which contains the old values to calculate the new values
    Eiw_old <- Eiw
    Eiw_new <- Eiw

    for (i in seq(nrow(sp))) {
        alpha <- sp$alpha[i]
        n_exp <- sp$n[i] # intake exponent
        p_exp <- sp$p[i] # metabolism exponent

        fc <- params@species_params$fc[i]
        f <- params@species_params$f0[i]

        # Calculate species specific metabolism based on the critical feeding
        # level in which Eriw. Should equal 0;
        # 0 = alpha * (1 - fc) * Eiw_old[i, ] - k_iw
        k_iw <- alpha * (1 - fc) * Eiw_old[i, ]

        # Calculate new Encounter Rate (E_new) to maintain growth at f
        # Growth = alpha * (1 - f) * E_new - k
        Eiw_new[i, ] <- (Eriw[i,] + k_iw) / (alpha * (1 - f))

        # Calculate Max Intake (h) to achieve feeding level f
        # f = E_new / (E_new + h)  => h = E_new * (1 - f) / f
        h_iw <- Eiw_new[i, ] * (1 - f) / f

        # Check that Eriw[i, ] = alpha * (1 - f) * Eiw_new[i, ] - k_iw
        if (!isTRUE(all.equal(alpha * (1 - f) * Eiw_new[i, ] - k_iw, Eriw[i, ]))) {
            stop(paste("Eriw was changed for species", rownames(sp)[i]))
        }

        sp$ks[i] <- k_iw[1] / (w[1]^p_exp)
        h <- h_iw[1] / (w[1]^n_exp)
        h[is.nan(h)] <- Inf
        sp$h[i]  <- h
        sp$Eiw[i] <- Eiw_new[i,1]/ (w[1]^n_exp)
    }

    # Set New External Encounter
    ext_encounter(params) <- Eiw_new

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
