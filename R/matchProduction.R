#' Match the production of the model to the Ecopath production
#'
#' This function sets the respiration rate according to the first master
#' equation of Ecopath, which is: \deqn{Q = P + R + U} where \eqn{Q} is the
#' consumption rate, \eqn{P} is the somatic production rate, \eqn{R} is the
#' respiration rate, and \eqn{U} is the unassimilated part of the consumption
#' rate, i.e., \eqn{U = (1 - \alpha) Q}, where \eqn{\alpha} is the assimilation
#' efficiency. Solving this for \eqn{R} we get: \deqn{R = Q \alpha - P}
#'
#' Ecopath does not explicitly account for the loss due to gonad production.
#' However most of the energy invested into gonad production is lost due to
#' the inefficiency of the reproductive process. Therefore in mizer we need
#' to include this loss in the respiration rate. Thus
#' \deqn{R = K + G - R_dd * w_0}
#' where \eqn{K} is the metabolic rate, \eqn{G} is the gonadic production rate,
#' \eqn{R_dd} is the rate of offspring production, and \eqn{w_0} is the weight
#' of an offspring at the size it enters the model. Thus \eqn{G - R_dd * w_0}
#' is the net loss due to gonad production.
#'
#'
#'
#' The `matchProduction` function calls `matchProductionOnce()`
#' repeatedly until the production matches within the tolerance specified by
#' the `tol` parameter.
#' The function will return a warning if the maximum number of iterations is
#' reached without converging.
#'
#' @param params A MizerParams object
#' @param tol The relative tolerance for the match. Default is 0.1.
#' @param max_iter The maximum number of iterations. Default is 10.
#'
#' @return A MizerParams object with the production matched
#' @family match functions
#' @export
matchProduction <- function(params, tol = 0.1, max_iter = 10) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }
    if (!hasName(sp, "gonad_proportion")) {
        stop("You must provide the gonad_proportion species parameter.")
    }

    for (i in seq_len(max_iter)) {
        # Break if tolerance achieved
        if (isProductionMatched(params, tol)) {
            break
        }
        params <- matchProductionOnce(params)
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Production converged in ", i - 1, " iterations.")
    }
    params
}

#' @rdname matchProduction
#' @param steady Whether to return the model to a steady state after adjusting
#'   the production. Default is TRUE.
#' @export
matchProductionOnce <- function(params, steady = TRUE) {
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

    # # Adjust gonadic production
    # Q <- getConsumption(params)
    # R_desired <- sp$alpha * Q - sp$ecopath_production
    # if (any(R_desired < 0)) {
    #     warning("Assimilated consumption not sufficient to cover production.")
    # }
    # G_desired <- sp$gonad_proportion * R_desired
    # G <- getGonadicProduction(params)
    # G_factor <- G_desired / G
    # sp <- set_species_param_default(sp, "w_repro_max", sp$w_max)
    # sp <- set_species_param_default(sp, "m", 1)
    # w_max_factor <- G_factor ^ (1 / (sp$n - sp$m))
    # sp$w_repro_max <- sp$w_repro_max * w_max_factor
    # if (any(sp$w_repro_max < sp$w_mat)) {
    #     stop("The gonadic proportion leads to a `w_repro_max` smaller than `w_mat`.")
    # }
    # species_params(params)$w_repro_max <- sp$w_repro_max

    # Adjust respiration
    Q <- getConsumption(params)
    G <- getGonadicProduction(params)
    K_desired <- sp$alpha * Q - sp$ecopath_production - G
    if (any(K_desired < 0)) {
        warning("Assimilated consumption not sufficient to cover production.")
    }
    K <- getRespiration(params)
    K_factor <- K_desired / K
    # If the metabolic rate was set manually, we need to rescale it, otherwise we
    # rescale the species parameters used to calculate the metabolic rate
    if (!is.null(comment(metab(params)))) {
        warning("This function has rescaled the metabolic rate that was set manually.")
        params@metab[] <- params@metab[] * K_factor
    } else {
        species_params(params)$ks <- sp$ks * K_factor
        if (hasName(sp, "k")) {
            species_params(params)$k <- sp$k * K_factor
        }
    }

    if (steady) {
        # Determine new steady state
        params <- params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    return(params)
}

#' @rdname matchProduction
isProductionMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    if (!hasName(sp, "ecopath_production")) {
        stop("You must provide the ecopath_production species parameter.")
    }

    # Calculate discrepancy in production
    Pratio <- sp$ecopath_production / getSomaticProduction(params)

    return(max(abs(Pratio - 1)) < tol)
}
