#' Set up intake rate and metabolic loss parameters for a tuned non-interacting
#' model. Take the existing tuned encounter rates and compute the correct
#' internal component values, and sets them based on a supplied assimilation
#' efficiency, critical feeding level,  and feeding level. Defaults for these
#' are set following Andersen (2019).
#'
#'
#' @param params A MizerParams object
#' @param fc Critical feeding level. Default is 0.2
#' @param f Feeding level. Default is 0.6
#' @return A MizerParams object with the adjusted max intake rate (h) and
#' metablolic loss (k) parameters
#' @export
computeEnergeticComponents <- function(params, fc = 0.2, f = 0.6) {
    # 1. Clear "comment" protection to allow recalculation
    comment(params@intake_max) <- NULL
    comment(params@metab)      <- NULL
    comment(params@search_vol) <- NULL
    comment(params@ext_encounter) <- NULL

    sp <- params@species_params
    w <- params@w

    # 2. Get the original encounter rates (before we add h or k)
    # At this stage, growth is simply: alpha * E_old
    E_old <- getEncounter(params)

    for (i in seq(nrow(sp))) {
        alpha <- sp$alpha[i]
        n_exp <- sp$n[i] # intake exponent
        p_exp <- sp$p[i] # metabolism exponent

        # Target Growth (E_r.i) we must preserve
        E_target_growth <- alpha * E_old[i, ]

        # STEP 1: Calculate Metabolism (k) based on critical feeding level fc
        # Based on theory: 0 = alpha * (1 - fc) * E_old - k
        k_iw <- alpha * (1 - fc) * E_old[i, ]

        # STEP 2: Calculate new Encounter Rate (E_new) to maintain growth at f
        # Growth = alpha * (1 - f) * E_new - k
        E_new_iw <- (E_target_growth + k_iw) / (alpha * (1 - f))

        # STEP 3: Calculate Max Intake (h) to achieve feeding level f
        # f = E_new / (E_new + h)  => h = E_new * (1 - f) / f
        h_iw <- E_new_iw * (1 - f) / f

        # STEP 4: Update coefficients in species_params
        # We pick a reference weight (e.g., 1g) to set the allometric constants
        ref_w <- 1
        # Find index closest to 1g, or just use allometric scaling logic
        sp$ks[i] <- k_iw[1] / (w[1]^p_exp)
        sp$h[i]  <- h_iw[1] / (w[1]^n_exp)

        # STEP 5: CRITICAL - Adjust gamma so the encounter rate actually changes
        # Since E = gamma * w^q * (Prey Abundance)
        # gamma_new = gamma_old * (E_new / E_old)
        sp$gamma[i] <- sp$gamma[i] * (E_new_iw[1] / E_old[i, 1])
    }

    # 3. Update the params object
    params@species_params <- sp

    # This triggers a full rebuild of the internal arrays based on the new sp constants
    params <- setParams(params)

    return(params)
}
