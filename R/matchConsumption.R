#' Match the consumption of the model to the Ecopath consumption
#'
#' This function adjusts the consumption parameters of the model to match the
#' Ecopath consumption It does this by rescaling the external mortality and
#' encounter rates as well as the maximum intake rate and the metabolic
#' respiration rate. This means that the ratio between growth and mortality
#' is preserved and hence the size spectrum is not changed.
#'
#' Besides rescaling the rates it also rescales the associated species parameters
#' `gamma`, `h`, and `ks` and `k`.
#'
#' The `matchConsumption` function calls `matchConsumptionOnce()`
#' repeatedly until the production matches within the tolerance specified by
#' the `tol` parameter.
#' The function will return a warning if the maximum number of iterations is
#' reached without converging.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance for the match. Default is 0.1.
#' @param max_iter The maximum number of iterations. Default is 10.
#'
#' @return A MizerParams object with the consumption matched
#' @family match functions
#' @export
matchConsumption <- function(params, tol = 0.1, max_iter = 10) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- set_species_param_default(params@species_params,
                                    "ecopath_consumption", NA)

    params |> steadySingleSpecies() |>
        matchBiomasses() |> steadySingleSpecies()

    for (i in seq_len(max_iter)) {
        if (isConsumptionMatched(params, tol)) {
            break
        }
        params <- matchConsumptionOnce(params)
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Consumption converged in ", i - 1, " iterations.")
    }
    params
}

#' @rdname matchConsumption
matchConsumptionOnce <- function(params) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- set_species_param_default(params@species_params,
                                    "ecopath_consumption", NA)

    Q <- getConsumption(params)
    Qratio <- sp$ecopath_consumption / Q

    # For species without `ecopath_consumption` use age at maturity
    sel <- is.na(Qratio)
    sp <- set_species_param_default(sp, "age_mat", NA)
    # If age at maturity is not specified, calculate it from von Bertalanffy
    if (all(c("k_vb", "w_inf") %in% names(sp))) {
        sp <- set_species_param_default(sp, "age_mat", age_mat_vB(params))
    }
    Qratio[sel] <- age_mat(params)[sel] / sp$age_mat[sel]

    # Don't affect species where no info is available
    Qratio[is.na(Qratio)] <- 1

    # Adjust consumption without changing size spectrum
    search_vol(params) <- search_vol(params) * Qratio
    intake_max(params) <- intake_max(params) * Qratio
    metab(params) <- metab(params) * Qratio
    ext_mort(params) <- ext_mort(params) * Qratio
    ext_encounter(params) <- ext_encounter(params) * Qratio
    # Also adjust species params
    params@species_params$gamma <- sp$gamma * Qratio
    params@species_params$h <- sp$h * Qratio
    params@species_params$ks <- sp$ks * Qratio
    if (hasName(sp, "k")) {
        params@species_params$k <- sp$k * Qratio
    }

    # Return new steady state
    params |> steadySingleSpecies() |>
        matchBiomasses() |> steadySingleSpecies()
}

#' @rdname matchConsumption
isConsumptionMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- set_species_param_default(params@species_params,
                                    "ecopath_consumption", NA)

    Q <- getConsumption(params)
    Qratio <- sp$ecopath_consumption / Q

    # For species without `ecopath_consumption` use age at maturity
    sel <- is.na(Qratio)
    sp <- set_species_param_default(sp, "age_mat", NA)
    # If age at maturity is not specified, calculate it from von Bertalanffy
    if (all(c("k_vb", "w_inf") %in% names(sp))) {
        sp <- set_species_param_default(sp, "age_mat", age_mat_vB(params))
    }
    Qratio[sel] <- age_mat(params)[sel] / sp$age_mat[sel]

    # Ignore species where no info is available
    Qratio[is.na(Qratio)] <- 1

    return(max(abs(Qratio - 1)) < tol)
}
