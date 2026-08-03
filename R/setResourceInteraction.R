#' Set resource interaction to absorb external encounter rate
#'
#' This function sets the resource interaction for each species to absorb as
#' much of the external encounter rate as possible. It adjusts the resource
#' capacity so that the resource can sustain the extra encounter without
#' changing the current resource abundance or the resource level.
#'
#' If the initial abundances were at steady state, then they will remain at
#' steady state. This is because the total encounter rate is not modified by
#' this function and because the extra mortality on the resource is exactly
#' balanced by the changed resource capacity.
#'
#' Optionally, you can specify the resource dynamics to use (e.g.,
#' "resource_semichemostat" or "resource_constant") by setting the
#' resource_dynamics argument. This allows you to switch between dynamic and
#' constant resource regimes directly through this function. For example:
#' `params <- setResourceInteraction(params, resource_dynamics = "resource_semichemostat")`
#' will set the resource dynamics to semichemostat.
#'
#' The amount that a species can absorb is limited by the requirement that the
#' extra resource encounter rate must not exceed the external encounter rate at
#' any size, because it is subtracted from the latter. This is imposed with a
#' tolerance: the extra resource encounter rate may exceed the external
#' encounter rate by up to `tol` times the total encounter rate at that size.
#' Hence the total encounter rate, and with it the steady state, changes by at
#' most a relative amount `tol`. Without this tolerance a single size at which
#' the external encounter rate is exactly zero (as produced by
#' [makeInteracting()] when predation on fish already covers the whole intake)
#' would prevent the species from absorbing anything at all.
#'
#' Because of this, restricting the resource to small sizes with `w_pp_cutoff`
#' helps: it makes the resource encounter rate of large predators negligible
#' compared to their total encounter rate, so their largest sizes no longer
#' constrain how much they can absorb at the sizes where they do feed on the
#' resource.
#'
#' Note: The function ensures that the resource level passed to setResource()
#' is strictly between 0 and 1, as required by mizer. If any value of
#' resource_level is greater than or equal to 1, it is set to just below 1
#' (0.99999999). If any value is less than or equal to 0, it is set to a small
#' positive value (1e-8). This prevents errors and ensures biological realism.
#'
#' @param params A mizer params object
#' @param resource_dynamics Optional. Character string specifying the resource
#'   dynamics to use. If NULL (default), the current resource dynamics are kept.
#' @param tol The relative change in the total encounter rate that is
#'   tolerated. Defaults to 1e-10. See Details.
#' @return The modified mizer params object with increased `interaction_resource`
#'   species parameters and correspondingly decreased external encounter rate.
#' @export
setResourceInteraction <- function(params, resource_dynamics = NULL,
                                   tol = 1e-10) {
    # Optionally set resource dynamics if specified by the user
    if (!is.null(resource_dynamics)) {
        resource_dynamics(params) <- resource_dynamics
    }
    # Save the resource level so we can restore it later
    # The resource_level must be strictly between 0 and 1 for setResource().
    # Values >= 1 are set to just below 1, and values <= 0 are set to a small positive value.
    resource_level <- resource_level(params)
    resource_level[resource_level >= 1] <- 0.99999999
    resource_level[resource_level <= 0] <- 1e-8

    old_resource_encounter <- getResourceEncounterRate(params)
    old_interaction_resource <- species_params(params)$interaction_resource
    old_encounter <- getEncounter(params)

    # Calculate the encounter rate achieved with interaction_resource = 1.
    # We do this on a copy so that `params` keeps the old value and the change
    # we make below is detected and recorded correctly.
    p1 <- params
    p1@species_params$interaction_resource <- 1
    encounter <- getResourceEncounterRate(p1)

    # Increasing `interaction_resource` by `r` adds `r * encounter` to the
    # resource encounter rate at each size, and that has to be taken out of
    # the external encounter rate. So the binding constraint is
    #     r * encounter <= ext_encounter    at every size.
    # We allow the extra resource encounter to overshoot the external
    # encounter rate by at most `tol` times the total encounter rate. Without
    # this allowance any size at which `ext_encounter` happens to be exactly
    # zero would force `r = 0`, no matter how minute the resource encounter
    # there is. With it, sizes where the resource contributes a negligible
    # part of the diet stop being a constraint, so restricting the resource to
    # small sizes (via `w_pp_cutoff`) frees up the large predators.
    allowance <- params@ext_encounter + tol * old_encounter
    ratio <- allowance / encounter
    # ratio is NaN where both allowance and encounter are zero, and Inf where
    # encounter is zero but allowance is not. Set these to Inf so they do not
    # affect the minimum calculation below.
    ratio[!is.finite(ratio)] <- Inf
    # The minimum ratio is the maximum that can be absorbed from the external
    # encounter rate
    ratio <- apply(ratio, 1, min)
    # If a species does not encounter the resource at any size:
    ratio[ratio == Inf] <- 0

    # Increase the resource interaction. `recalculate = FALSE` records the new
    # value, so that a later recalculation does not reset it to its default,
    # without recalculating any rates. We set the rate that depends on it, the
    # external encounter rate, ourselves below.
    sp_new <- params@species_params
    sp_new$interaction_resource <- old_interaction_resource + ratio
    species_params(params, recalculate = FALSE) <- sp_new
    new_resource_encounter <- getResourceEncounterRate(params)

    # Subtract the absorbed encounter rate from the external encounter rate
    params@ext_encounter <- params@ext_encounter - new_resource_encounter +
        old_resource_encounter

    # Rounding errors could give negative values
    params@ext_encounter[params@ext_encounter < 0] <- 0

    # The construction above bounds the relative change in the total encounter
    # rate by `tol`, and `getResourceEncounterRate()` computes the resource
    # encounter rate the same way `getEncounter()` does, so this should not
    # trigger. It is kept as a guard because the bound is easy to invalidate.
    rel_diff <- abs(getEncounter(params) - old_encounter) / old_encounter
    rel_diff[!is.finite(rel_diff)] <- 0
    max_rel_diff <- apply(rel_diff, 1, max)
    affected <- max_rel_diff > 2 * tol
    if (any(affected)) {
        warning("The encounter rate changed for the following species: ",
                paste0(names(max_rel_diff)[affected], " (",
                       signif(max_rel_diff[affected], 3), ")",
                       collapse = ", "),
                ". The numbers in brackets give the maximum relative change ",
                "in the encounter rate over all sizes. If these are large ",
                "then the model may no longer be at steady state.")
    }

    # Set the resource capacity so that the the steady-state resource
    # abundance stays the same
    params <- setResource(params,
                          resource_level = resource_level)

    return(params)
}

#' Get resource encounter rate
#'
#' This function returns the part of the encounter rate that comes from feeding
#' on the resource, i.e., the contribution that `interaction_resource` scales.
#'
#' The calculation mirrors that in [mizer::mizerEncounter()] exactly, with the
#' fish prey removed from the prey spectrum, so that the result agrees with
#' [mizer::getEncounter()] to machine precision. In particular it uses the same
#' fast Fourier transform to evaluate the convolution integral, unless the
#' predation kernel has been set by hand, in which case mizer sums over the
#' predation kernel array instead and so do we.
#'
#' Using the same method matters because functions like
#' [setResourceInteraction()] subtract this rate from the external encounter
#' rate and rely on the total encounter rate being left unchanged. Note that
#' mizer's Fourier transform evaluates a circular convolution, so at the very
#' smallest sizes the result includes a small spurious contribution from prey at
#' the top of the weight grid. That artefact is reproduced here on purpose,
#' because it is present in the encounter rate that drives the model.
#'
#' @param params A mizer params object
#' @return A matrix (species x size) of resource encounter rates
#' @importFrom stats mvfft
#' @export
getResourceEncounterRate <- function(params) {
    n_pp <- initialNResource(params)
    # idx_sp are the index values of params@w_full such that
    # params@w_full[idx_sp] = params@w
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

    # If the user has set a custom predation kernel then mizer can not use the
    # fft and neither can we.
    if (!is.null(comment(params@pred_kernel))) {
        phi_resource <- params@species_params$interaction_resource *
            rowSums(sweep(
                params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
                "*",
                check.margin = FALSE
            ), dims = 2)
        return(params@search_vol * phi_resource)
    }

    prey <- outer(params@species_params$interaction_resource, n_pp)
    prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
    avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                         mvfft(base::t(prey)),
                                     inverse = TRUE))) /
        length(params@w_full)
    avail_energy <- avail_energy[, idx_sp, drop = FALSE]
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    avail_energy[avail_energy < 1e-18] <- 0

    encounter <- params@search_vol * avail_energy
    dimnames(encounter) <- dimnames(params@search_vol)
    return(encounter)
}
