#' Coefficient of the external-diffusion power law for one species
#'
#' The growth-diffusion rate is modelled as a power law
#' \deqn{d(w) = D_{ext}\, w^{n+1},}
#' stored in the `ext_diffusion` slot of the [MizerParams][mizer::MizerParams]
#' object. This helper returns the coefficient `D_ext` for a single species.
#'
#' The `ext_diffusion` slot is the diffusion that the model dynamics actually
#' use, so it is treated as the source of truth: when it is non-zero the
#' coefficient is read off the stored rate. Otherwise the value is taken from the
#' `D_ext` species parameter, or, failing that, from the legacy `d_over_g`
#' parameter via \eqn{D_{ext} = (d/g)\, g_0} with \eqn{g_0 = g(w_1)/w_1^{n}}.
#' When no diffusion information is available it returns 0.
#'
#' @param params A MizerParams object.
#' @param sp_select A logical vector selecting a single species.
#' @param n The exponent of the growth power law for that species.
#' @return A single numeric value, the `D_ext` coefficient.
#' @keywords internal
diffusion_coefficient <- function(params, sp_select, n) {
    w <- params@w
    ed <- params@ext_diffusion[sp_select, ]
    if (any(ed > 0)) {
        # The slot holds a clean power law, so any positive bin recovers D_ext.
        j <- max(which(ed > 0))
        return(unname(ed[j] / w[j]^(n + 1)))
    }
    sp <- params@species_params
    if ("D_ext" %in% names(sp) && !is.na(sp[["D_ext"]][sp_select]) &&
        sp[["D_ext"]][sp_select] > 0) {
        return(sp[["D_ext"]][sp_select])
    }
    if ("d_over_g" %in% names(sp) && !is.na(sp[["d_over_g"]][sp_select])) {
        g <- getEGrowth(params)[sp_select, ]
        g0 <- g[1] / w[1]^n
        return(unname(sp[["d_over_g"]][sp_select] * g0))
    }
    0
}
