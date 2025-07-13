#' Get diffusion rate
#'
#' The diffusion rate \eqn{d_i(w)} for species \eqn{i}
#' and weight \eqn{w} has contributions from the encounter of
#' fish prey and of resource. This is determined by summing over all prey
#' species and the resource spectrum and then integrating over all prey sizes
#' \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:
#' \deqn{
#' d_i(w) = (1-f_i(w))(\alpha_i(1-\psi_i(w)))^2\gamma_i(w) \int
#' \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) \right)
#' \phi_i(w,w_p) w_p^2 \, dw_p.
#' }{(1-f_i(w))(\alpha_i(1-\psi_i(w)))^2\gamma_i(w) \int
#' ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) )
#' \phi_i(w,w_p) w_p^2 dw_p.}
#' Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#' \eqn{N_R(w)} is the abundance density of resource.
#' The overall prefactor \eqn{\gamma_i(w)} determines the predation power of the
#' predator. It could be interpreted as a search volume and is set with the
#' [setSearchVolume()] function. The predation kernel
#' \eqn{\phi(w,w_p)} is set with the [setPredKernel()] function. The
#' species interaction matrix \eqn{\theta_{ij}} is set with [setInteraction()]
#' and the resource interaction vector \eqn{\theta_{ip}} is taken from the
#' `interaction_resource` column in `params@species_params`.
#' \eqn{f(w)} is the feeding level calculated with
#' [getFeedingLevel()]. \eqn{\psi(w)} is the proportion of the available energy
#' that is invested in reproduction instead of growth, as stored in `params@psi`.
#'
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be used, as they are meaningless.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#' @param order By setting this to integer values above 2 one can get the
#'   coefficients of higher-order terms in the expansion of the jump-growth
#'   equation. The default is 2, which gives the diffusion rate.
#' @param ... Unused
#'
#' @return A named two dimensional array (predator species x predator size) with
#'   the diffusion rates.
#' @export
#' @family mizer rate functions
getDiffusion <- function(params, n = initialN(params),
                         n_pp = initialNResource(params),
                         n_other = initialNOther(params),
                         t = 0, order = 2, ...) {

    # idx_sp are the index values of params@w_full such that
    # params@w_full[idx_sp] = params@w
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

    # If the the user has set a custom pred_kernel we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (!is.null(comment(params@pred_kernel))) {
        # n_eff_prey is the total prey abundance by size exposed to each
        # predator (prey not broken into species - here we are just working out
        # how much a predator eats - not which species are being eaten - that is
        # in the mortality calculation
        # \sum_j \theta_{ij} N_j(w_p) w_p dw_p
        n_eff_prey <- sweep(params@interaction %*% n, 2,
                            params@w^order * params@dw, "*", check.margin = FALSE)
        # pred_kernel is predator species x predator size x prey size
        # So multiply 3rd dimension of pred_kernel by the prey biomass density
        # Then sum over 3rd dimension to get consumption rate of each predator by
        # predator size
        # This line is a bottle neck
        phi_prey_species <- rowSums(sweep(
            params@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
        # Eating the background
        # This line is a bottle neck
        phi_prey_background <- params@species_params$interaction_resource *
            rowSums(sweep(
                params@pred_kernel, 3, params@dw_full * params@w_full^order * n_pp,
                "*", check.margin = FALSE), dims = 2)
        encounter <- params@search_vol * (phi_prey_species + phi_prey_background)
    } else {
        stop("Fourier transform method not implemented for diffusion rate.")
    }

    diffusion <- (1 - getFeedingLevel(params, n = n, n_pp = n_pp,
                                      n_other = n_other, t = t)) *
        (params@species_params$alpha)^order * encounter
    return(diffusion)
}
