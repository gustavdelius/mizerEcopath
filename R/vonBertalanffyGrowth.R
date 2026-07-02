#' von Bertalanffy growth rate
#'
#' Return an ArraySpeciesBySize array with the von Bertalanffy growth rates
#' calculated from the `w_inf` and `k_vb` species parameters.
#'
#' The von Bertalanffy growth rate in terms of length \eqn{l} is
#' \deqn{dl/dt = k_{vb}(L_\inf - l(t)).}
#' After converting this to weight \eqn{w = al^b} this becomes
#' \deqn{dw/dt = Aw^n-Bw}
#' where \eqn{n=1-1/b}, \eqn{A=b\, k_{vb}\, w_{\inf}^{1/b}} and
#' \eqn{B=b\, k_{vb}}. The von Bertalanffy parameters \eqn{w_{\inf}} and
#' \eqn{k_{vb}} and the weight-length conversion parameters \eqn{a} and \eqn{b}
#' need to be provided in the `w_inf`, `k_vb`, `a` and `b` columns of the
#' `species_params` data frame contained in the MizerParams object. If `w_inf`
#' is missing then `w_max` is used instead and if `b` is missing then a default
#' of `b = 3` is used.
#'
#' The growth rate is evaluated on the size grid `w` of the MizerParams object
#' and is set to zero at sizes above `w_inf`, where the von Bertalanffy
#' expression would otherwise turn negative.
#'
#' @param object A `MizerParams` object.
#'
#' @return An array (species x size) with the von Bertalanffy growth rates.
#' @export
vonBertalanffyGrowth <- function(object) {
    UseMethod("vonBertalanffyGrowth")
}

#' @rdname vonBertalanffyGrowth
#' @export
vonBertalanffyGrowth.MizerParams <- function(object) {
    sp <- object@species_params
    sp <- set_species_param_default(sp, "b", 3)
    sp <- set_species_param_default(sp, "k_vb", NA)
    sp <- set_species_param_default(sp, "w_inf", sp$w_max)

    n <- 1 - 1 / sp$b
    A <- sp$b * sp$k_vb * sp$w_inf ^ (1 / sp$b)
    B <- sp$b * sp$k_vb

    w <- object@w
    # g[i, j] = A[i] * w[j]^n[i] - B[i] * w[j]
    g <- A * t(outer(w, n, "^")) - outer(B, w, "*")
    g[g < 0] <- 0

    dimnames(g) <- list(sp = sp$species, w = signif(w, 3))
    g
}
