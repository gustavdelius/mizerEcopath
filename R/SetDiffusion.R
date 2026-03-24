#' Set diffusion rate
#'
#' @param object A MizerParams object
#' @param diffusion Optional. An array (species x size) holding the diffusion rate.
#' @param reset If set to TRUE, the diffusion rate will be reset to the
#'   calculated value even if previously overwritten manually.
#' @param ... Unused
#'
#' @return A MizerParams object with updated diffusion rate.
#' @export
setDiffusion <- function(object, diffusion = NULL, reset = FALSE, ...) {
    UseMethod("setDiffusion")
}

#' @export
setDiffusion.MizerParams <- function(object, diffusion = NULL, reset = FALSE, ...) {
    assert_that(is.flag(reset))
    params <- object
    species_params <- params@species_params

    if (reset) {
        if (!is.null(diffusion)) {
            warning("Because you set `reset = TRUE`, the value you provided for `diffusion` will be ignored.")
            diffusion <- NULL
        }
        comment(params@diffusion) <- NULL
    }

    if (!is.null(diffusion)) {
        if (is.null(comment(diffusion))) {
            comment(diffusion) <- if (is.null(comment(params@diffusion))) "set manually" else comment(params@diffusion)
        }
        params@diffusion[] <- diffusion
        comment(params@diffusion) <- comment(diffusion)
        params@time_modified <- lubridate::now()
        return(params)
    }

    # Ensure parameter column exists
    params <- set_species_param_default(params, "d_over_g", 0)

    # CORRECTED MATH
    growth <- getEGrowth(params)
    # g_0 is growth at first size bin divided by w_min^n
    g_0 <- growth[, 1] / (params@w[1]^params@species_params$n)
    d_0 <- params@species_params$d_over_g * g_0

    # Calculate w^(n+1) matrix:
    # Base = weights, Exponent = n+1
    w_pow <- outer(params@w, params@species_params$n + 1, "^")

    # Resulting matrix: Species (rows) x Weights (cols)
    new_diffusion <- d_0 * t(w_pow)

    if (!is.null(comment(params@diffusion))) {
        if (different(new_diffusion, params@diffusion)) {
            message("The diffusion rate has been set manually and will not be recalculated.")
        }
        return(params)
    }

    params@diffusion[] <- new_diffusion
    params@time_modified <- lubridate::now()
    return(params)
}
