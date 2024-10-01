#' Match the observed yield by adjusting catchability
#'
#' This function adjusts the catchability parameters of the model to match the
#' observed yield.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' Note that adjusting the `catchability` parameter to match the yield may not
#' lead to a convergent iteration because it currently assumes that increasing
#' the catchability will increase the yield. This may not be the case if the
#' increased mortality on large individuals truncates the size spectrum too
#' much.
#'
#' The `matchYield` function calls `matchYieldOnce()` repeatedly until the
#' production matches within the tolerance specified by the `tol` parameter.
#' The function will return a warning if the maximum number of iterations is
#' reached without converging.
#'
#' @param params A MizerParams object
#' @return A MizerParams object with the catchability adjusted to match the
#'   observed yield
#' @family match functions
#' @export
matchYield <- function(params, tol = 0.1, max_iter = 10) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    for (i in seq_len(max_iter)) {
        # Break if tolerance achieved
        if (isYieldMatched(params, tol)) {
            break
        }
        params <- matchYieldOnce(params)
    }
    if (i == max_iter) {
        warning("Did not converge.")
    } else {
        message("Yield converged in ", i - 1, " iterations.")
    }
    return(params)
}

#' @rdname matchYield
#' @param steady Whether to return the model to a steady state after adjusting
#'   the production. Default is TRUE.
#' @export
matchYieldOnce <- function(params, steady = TRUE) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }

    # Adjust fishing mortality
    # TODO: I am just lazy here at the moment. Code for making yield match
    # for each gear and each species will probably require looping
    if (nrow(gp) != nrow(sp) ||
        !all(gp$species == sp$species)) {
        stop("This code currently requires a single gear for each species.")
    }
    Cratio <- gp$yield_observed / getYield(params)
    sel <- !is.na(Cratio)
    gp$catchability[sel] <- gp$catchability[sel] * Cratio
    gear_params(params)$catchability <- gp$catchability

    if (steady) {
        # Determine new steady state
        params <- params |> steadySingleSpecies() |>
            matchBiomasses() |> steadySingleSpecies()
    }
    params@time_modified <- lubridate::now()
    return(params)
}


#' @rdname matchYield
isYieldMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }
    # Calculate discrepancy in yields
    Cratio <- gp$yield_observed / getYield(params)

    return(max(abs(Cratio - 1)) < tol)
}
