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
#' reduction equals the metabolic loss itself at `w_inf` (all income then funds
#' reproduction) and is smaller below, so the metabolic loss stays
#' non-negative.
#'
#' @section Negative `t0`:
#' A negative `t0` means the von Bertalanffy fish already has a positive size at
#' age 0, so mizer, which starts every fish at `w_min`, cannot follow the curve
#' below maturity. Instead the growth below maturity is sped up so that maturity
#' is reached at the age the von Bertalanffy curve (with its negative `t0`)
#' predicts, namely \eqn{t_{\mathrm{mat}} = t_0 - \log(1 -
#' (w_{\mathrm{mat}}/w_\infty)^{1/b})/k_{vb}}; above maturity the growth is left
#' unchanged and so rejoins the von Bertalanffy curve. The speed-up is applied
#' in two stages:
#' \itemize{
#'   \item The metabolic loss below maturity is scaled by a constant factor
#'     \eqn{c \in (0, 1]}, which keeps the growth a von Bertalanffy form with
#'     rate \eqn{c\,k_{vb}} and asymptote \eqn{w_\infty^{1/b}/c}. The factor is
#'     found by inverting the resulting age-at-maturity for `c`.
#'   \item Removing the metabolic loss entirely (\eqn{c = 0}) gives the fastest
#'     growth this can achieve and hence a minimum reachable maturity age. If
#'     `t0` is so negative that even this is too slow, the metabolic loss is set
#'     to zero and the intake (external encounter rate) below maturity is
#'     additionally raised by the factor needed to hit the target maturity age.
#'     A message reports the factor used. Note that this raises the implied
#'     juvenile consumption above the power law for the affected species.
#' }
#' Species with `t0 >= 0` (or missing `t0`) are left unchanged.
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

    # Subtract the reproduction investment from the metabolic loss rate. Somatic
    # growth in mizer is g = (1 - psi) * (alpha * E - metab), so for g to equal
    # the von Bertalanffy rate g_vb we need alpha * E - metab = g_vb / (1 - psi).
    # The juvenile power law already gives alpha * E - metab = g_vb, so we reduce
    # metab by psi / (1 - psi) * g_vb. With w_repro_max = w_inf this reduction
    # equals metab itself at w_inf (all income then funds reproduction) and is
    # smaller below, so metab stays non-negative. The pmax() guards against tiny
    # negative values from floating-point rounding at the top of the range.
    vb <- vonBertalanffyGrowth(p)
    repro <- p@psi / (1 - p@psi) * vb
    repro[vb == 0] <- 0

    # Handle a negative von Bertalanffy birth age t0. A negative t0 means the
    # fish already has a positive size at age 0, so mizer, which starts every
    # fish at w_min, cannot follow the von Bertalanffy curve below maturity. We
    # instead speed up growth below maturity so that maturity is reached at the
    # age the von Bertalanffy curve (with its negative t0) predicts. Above
    # maturity the growth is left unchanged, so the curve rejoins the von
    # Bertalanffy curve there.
    #
    # In terms of u = w^(1/b) the von Bertalanffy rate is du/dt = k_vb(u_inf - u)
    # with u_inf = w_inf^(1/b). Scaling metab by a factor c turns this into
    # du/dt = c*k_vb(u_inf/c - u), i.e. again von Bertalanffy but with rate
    # c*k_vb and asymptote u_inf/c. The age to grow from w_min to w_mat is then
    #   age(c) = log((u_inf - c*u_min)/(u_inf - c*u_mat)) / (c*k_vb),
    # which we invert for c in (0, 1] to match the target maturity age. The
    # fastest growth reducing metab alone can achieve is at c = 0 (no metabolic
    # loss), giving the minimum maturity age floor_age below. When even that is
    # too slow (a strongly negative t0), we additionally raise the intake below
    # maturity: with metab = 0 the growth is the pure power law d*A*w^n, which
    # has constant du/dt and maturity age floor_age / d, so scaling the intake
    # by d = floor_age / t_mat hits the target.
    metab_base <- metab(p)
    for (i in seq_len(nrow(sp))) {
        if (is.na(sp$t0[i]) || sp$t0[i] >= 0) next
        u_inf <- sp$w_inf[i] ^ (1 / sp$b[i])
        u_min <- p@w[p@w_min_idx[i]] ^ (1 / sp$b[i])
        u_mat <- sp$w_mat[i] ^ (1 / sp$b[i])
        # von Bertalanffy age at maturity with the given (negative) t0
        t_mat <- sp$t0[i] - log(1 - u_mat / u_inf) / sp$k_vb[i]
        age <- function(c) {
            log((u_inf - c * u_min) / (u_inf - c * u_mat)) / (c * sp$k_vb[i])
        }
        # Smallest maturity age reducing metab alone can reach (c -> 0 limit)
        floor_age <- (u_mat - u_min) / (sp$k_vb[i] * u_inf)
        below <- p@w < sp$w_mat[i]
        if (t_mat > floor_age) {
            # Reachable by reducing the metabolic loss alone.
            c <- stats::uniroot(function(c) age(c) - t_mat,
                                interval = c(1e-8, 1))$root
            metab_base[i, below] <- c * metab_base[i, below]
        } else {
            # Zero metabolic loss is not fast enough, so also raise the intake.
            metab_base[i, below] <- 0
            d <- floor_age / t_mat
            ext_encounter(p)[i, below] <- d * ext_encounter(p)[i, below]
            message("Species ", sp$species[i], ": t0 = ", sp$t0[i],
                    " requires faster growth; ",
                    "intake below maturity raised by a factor of ",
                    signif(d, 3), " to reach maturity at the von Bertalanffy age.")
        }
        repro[i, below] <- 0
    }

    metab(p) <- pmax(metab_base - repro, 0)

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
