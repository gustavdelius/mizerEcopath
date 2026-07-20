#' Create model with von Bertalanffy growth rates in the steady state
#'
#' Builds a single- or multi-species [MizerParams][mizer::MizerParams-class]
#' object in which each species grows along its von Bertalanffy curve while at
#' steady state. The von Bertalanffy parameters are supplied in the `w_inf`,
#' `k_vb` and `t0` columns of `species_params`, together with the weight-length
#' exponent `b`, the maturity size `w_mat` and the maximum size `w_max`.
#'
#' @section Growth:
#' The von Bertalanffy growth in length, \eqn{dl/dt = k_{vb}(L_\infty - l)},
#' becomes in weight \eqn{dw/dt = A\,w^n - B\,w} with \eqn{n = 1 - 1/b},
#' \eqn{A = b\,k_{vb}\,w_\infty^{1/b}} and \eqn{B = b\,k_{vb}} (see
#' [vonBertalanffyGrowth()]). The model reproduces this rate by
#' \itemize{
#'   \item switching off satiation, predation and resource interactions so that
#'     the energy income is a pure power law \eqn{\alpha E = A\,w^n},
#'   \item setting the metabolic loss rate to \eqn{B\,w} (so
#'     \eqn{\alpha E - \mathrm{metab} = A\,w^n - B\,w}), and
#'   \item switching off constant mortality and using a power-law external
#'     mortality chosen to give the juvenile biomass density a slight negative
#'     slope.
#' }
#'
#' @section Reproduction:
#' Mizer computes somatic growth as \eqn{g = (1 - \psi)(\alpha E -
#' \mathrm{metab})}, where \eqn{\psi} is the fraction of surplus energy invested
#' into reproduction. For the *somatic* growth to stay on the von Bertalanffy
#' curve the surplus energy must therefore be \eqn{g_{vb}/(1 - \psi)}, so the
#' metabolic loss is reduced by \eqn{\psi/(1 - \psi)\,g_{vb}}. Because
#' `w_repro_max` is set to `w_inf`, and mizer's maturity ogive gives
#' \eqn{\psi \propto (w/w_{\mathrm{repro\,max}})^{1-n}} above `w_mat`, this
#' reduction equals the metabolic loss across the mature range where the ogive
#' has saturated (reducing the metabolic loss to zero there, so the income funds
#' both somatic growth and reproduction) and is smaller in the maturation
#' transition just above `w_mat`; either way it never exceeds the metabolic loss,
#' which therefore stays non-negative.
#'
#' @section Negative `t0`:
#' A negative `t0` means the von Bertalanffy fish already has a positive size at
#' age 0, so mizer, which starts every fish at `w_min`, cannot follow the curve
#' below maturity. Instead a speed-up
#' \eqn{\Delta(w) = s\,w^n(1 - (w/w_{\mathrm{mat}})^{1-n})} is added to the growth
#' rate below maturity, chosen so that maturity is reached at the age the von
#' Bertalanffy curve (with its negative `t0`) predicts, namely
#' \eqn{t_{\mathrm{mat}} = t_0 - \log(1 - (w_{\mathrm{mat}}/w_\infty)^{1/b})/k_{vb}}.
#' Because \eqn{\Delta} vanishes at `w_mat`, the growth rate stays continuous
#' there, and above maturity the growth is unchanged and so coincides with the
#' von Bertalanffy curve.
#'
#' Since \eqn{g_{vb} + \Delta = (A + s)\,w^n - (B + s/w_{\mathrm{mat}}^{1/b})\,w}
#' is again of von Bertalanffy form, the age to maturity has a closed form which
#' is inverted for the amplitude \eqn{s}. Its \eqn{w^n} term is added to the
#' encounter rate (income) and its \eqn{w} term to the metabolic loss, keeping
#' consumption allometric (\eqn{\propto w^n}) and respiration linear
#' (\eqn{\propto w}). A message notes that juvenile intake is raised for the
#' affected species. An error is raised only if `t0` is so negative that
#' \eqn{t_{\mathrm{mat}} \le 0}. Species with `t0 >= 0` (or missing `t0`) are
#' left unchanged.
#'
#' @param species_params Species parameter data frame. Must contain the von
#'   Bertalanffy parameters `w_inf`, `k_vb` and `t0` as well as the usual mizer
#'   parameters such as `b`, `w_mat` and `w_max`.
#' @param no_w Number of size bins. Default 200
#' @param max_w Largest size in the model. Defaults to the largest `w_max`.
#' @return A MizerParams object with the model at steady state
#' @seealso [vonBertalanffyGrowth()]
#' @export
newVonBertalanffyParams <- function(species_params, no_w = 200, max_w = NULL) {
    sp <- validGivenSpeciesParams(species_params)

    # Impose relation between exponents
    sp <- set_species_param_default(sp, "n", 1 - 1 / sp$b)
    sp$p <- sp$n
    sp <- set_species_param_default(sp, "d", sp$n - 1)

    # Set default assimilation efficiency
    sp <- set_species_param_default(sp, "alpha", 0.8)

    # Switch off metabolic respiration
    sp$ks <- 0
    # Switch off constant mortality
    sp$z0 <- 0

    # Generate a default mizer model with the desired species We extend the
    # resource spectrum over the entire size range to ensure that all species
    # have sufficient prey throughout their life.

    if (is.null(max_w)) {max_w <- max(sp$w_max)}

    if (max_w < max(sp$w_max)) {
        warning("The maximum weight provided (max_w) is lower than the
                maximum size of the fish. The model has been generated with
                the latter")
    }

    max_w <- max(max_w, sp$w_max)

    p <- newMultispeciesParams(sp, no_w = no_w, info_level = 0,
                               # extend resource over entire size range
                               max_w = max_w,
                               w_pp_cutoff = max_w * (1 + 1e-9),
                               resource_dynamics = "resource_constant",
                               n = 0.7)
    sp <- p@species_params

    # Switch off all interactions
    interaction_matrix(p)[] <- 0
    sp$interaction_resource <- 0

    # Switch off satiation
    sp$h <- Inf
    intake_max(p)[] <- Inf

    # Set constant metabolic rate
    sp$k <- sp$b * sp$k_vb

    # Set power-law encounter rate (the coefficient will be adjusted below)
    sp$E_ext <- sp$b * sp$k_vb * sp$w_inf^(1/sp$b) / sp$alpha

    species_params(p) <- sp

    # Somatic growth in mizer is g = (1 - psi) * (alpha * E - metab), where psi
    # is the fraction of surplus energy invested into reproduction. We want the
    # somatic growth to equal a target rate g_target, which is the von
    # Bertalanffy rate g_vb plus, for species with a negative t0, a speed-up
    # Delta below maturity (built below). This requires
    # alpha * E - metab = g_target / (1 - psi), i.e. the metabolic loss is
    # reduced by the reproduction investment psi / (1 - psi) * g_target. With
    # w_repro_max = w_inf this reduction never exceeds metab, so metab stays
    # non-negative. See the "von Bertalanffy growth" vignette for the derivation.
    vb <- vonBertalanffyGrowth(p)
    psi <- p@psi
    w <- p@w

    # Speed-up for a negative von Bertalanffy birth age t0. A negative t0 means
    # the fish already has a positive size at age 0, so mizer, which starts every
    # fish at w_min, cannot follow the von Bertalanffy curve below maturity. We
    # add Delta(w) = s * w^n * (1 - (w/w_mat)^(1 - n)) to the growth rate below
    # maturity, chosen so that maturity is reached at the age the von Bertalanffy
    # curve (with its negative t0) predicts. Delta vanishes at w_mat, so the
    # growth rate stays continuous there, and above maturity the growth is
    # unchanged. Since g_vb + Delta = (A + s) w^n - (B + s/u_mat) w is again a
    # von Bertalanffy form, the age to maturity has a closed form which we invert
    # for s. Its w^n term is supplied by the encounter (income) and its w term by
    # the metabolic loss, keeping consumption allometric (~ w^n) and respiration
    # linear (~ w).
    Delta <- matrix(0, nrow(sp), length(w))
    metab_extra <- matrix(0, nrow(sp), length(w))
    enc_extra <- matrix(0, nrow(sp), length(w))
    for (i in seq_len(nrow(sp))) {
        if (is.na(sp$t0[i]) || sp$t0[i] >= 0) next
        b_i <- sp$b[i]
        n_i <- 1 - 1 / b_i
        A_i <- b_i * sp$k_vb[i] * sp$w_inf[i] ^ (1 / b_i)
        B_i <- b_i * sp$k_vb[i]
        u_min <- p@w[p@w_min_idx[i]] ^ (1 / b_i)
        u_mat <- sp$w_mat[i] ^ (1 / b_i)
        u_inf <- sp$w_inf[i] ^ (1 / b_i)
        # von Bertalanffy age at maturity with the given (negative) t0
        t_mat <- sp$t0[i] - log(1 - u_mat / u_inf) / sp$k_vb[i]
        if (t_mat <= 0) {
            stop("Species ", sp$species[i], ": t0 = ", sp$t0[i], " is too ",
                 "negative; maturity cannot be reached in positive time.")
        }
        # Age from w_min to w_mat for the modified von Bertalanffy rate
        # (A + s) w^n - (B + s/u_mat) w, as a function of the amplitude s.
        age <- function(s) {
            B_t <- B_i + s / u_mat
            u_inf_t <- (A_i + s) / B_t
            log((u_inf_t - u_min) / (u_inf_t - u_mat)) / (B_t / b_i)
        }
        # age(s) decreases monotonically from the plain von Bertalanffy transit
        # time at s = 0 towards 0. If maturity is already reached in time no
        # speed-up is needed; otherwise bracket and solve age(s) = t_mat.
        if (t_mat >= age(0)) next
        s_hi <- 1
        while (age(s_hi) > t_mat) s_hi <- s_hi * 2
        s <- stats::uniroot(function(s) age(s) - t_mat, c(0, s_hi))$root

        below <- w < sp$w_mat[i]
        wb <- w[below]
        Delta[i, below] <- s * wb ^ n_i * (1 - (wb / sp$w_mat[i]) ^ (1 - n_i))
        enc_extra[i, below] <- s * wb ^ n_i / sp$alpha[i]   # w^n term -> income
        metab_extra[i, below] <- (s / u_mat) * wb           # w term -> metab
        message("Species ", sp$species[i], ": t0 = ", sp$t0[i], " < 0; growth ",
                "below maturity sped up to reach maturity at the von Bertalanffy ",
                "age (this raises juvenile intake).")
    }

    # Reduce the metabolic loss by the reproduction investment for the target
    # growth. The pmax() guards against tiny negative values from floating-point
    # rounding near w_inf, where the reduction equals the metabolic loss.
    g_target <- vb + Delta
    repro <- psi / (1 - psi) * g_target
    repro[vb == 0] <- 0
    metab(p) <- pmax(metab(p) + metab_extra - repro, 0)
    ext_encounter(p) <- ext_encounter(p) + enc_extra

    # Set power-law mortality
    # Choose a positive coefficient so that the juvenile biomass density
    # has a slightly negative slope (of -0.2).
    sp$z_ext <- sp$alpha * sp$E_ext * (1 + 0.2 - sp$n)

    species_params(p) <- sp
    # Match Biomasses
    p <- matchBiomasses(p)
    # Set to steady state
    p <- steadySingleSpecies(p, keep = "biomass")
    p <- setBevertonHolt(p, reproduction_level = 0)

    return(p)
}
