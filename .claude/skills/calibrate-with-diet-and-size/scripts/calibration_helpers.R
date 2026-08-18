# Calibration Helper Functions for Multi-Stage Mizer Model Fitting
#
# These helper functions support the invariance-preserving calibration workflow
# described in the `calibrate-with-diet-and-size` skill.

#' Set Feeding Level in an Interacting Model Without Moving Steady State
#'
#' Rescales search volume, external encounter, and maximum intake so that the
#' feeding level becomes `f_new` while keeping steady-state consumption,
#' predation mortality, and diet composition unchanged.
#'
#' @param params A MizerParams object.
#' @param f_new Named numeric vector of target feeding levels per species (or a scalar applied to all).
#' @return An updated MizerParams object with identical steady state.
#' @export
setFeedingLevelInteracting <- function(params, f_new) {
    f_old <- mizer::getFeedingLevel(params)
    sp <- mizer::species_params(params)$species
    if (is.null(names(f_new))) {
        f_new <- stats::setNames(rep(f_new, length(sp))[seq_along(sp)], sp)
    }
    for (s in names(f_new)) {
        i <- match(s, sp)
        if (is.na(i)) next
        cc <- (1 - f_old[i, ]) / (1 - f_new[[s]])
        sv <- mizer::search_vol(params)
        sv[i, ] <- sv[i, ] * cc
        mizer::search_vol(params) <- sv
        
        ee <- mizer::ext_encounter(params)
        ee[i, ] <- ee[i, ] * cc
        mizer::ext_encounter(params) <- ee
        
        E <- mizer::getEncounter(params)
        im <- mizer::intake_max(params)
        im[i, ] <- E[i, ] * (1 - f_new[[s]]) / f_new[[s]]
        mizer::intake_max(params) <- im
    }
    params
}

#' Reallocate Interaction Links to External Pools Without Moving Steady State
#'
#' Multiplies predator-prey interaction matrix entries by `lambda` and reallocates
#' the removed food to external encounter and removed predation to external mortality.
#' Leaves growth, mortality, spectra, and catch size distributions invariant.
#'
#' @param params A MizerParams object.
#' @param pred Character vector of predator species names.
#' @param prey Character vector of prey species names.
#' @param lambda Numeric attenuation factor (e.g. 0.01 for near-complete decoupling).
#' @return An updated MizerParams object with identical steady state.
#' @export
set_link <- function(params, pred, prey, lambda) {
    mort_old <- mizer::getPredMort(params)
    enc_old  <- mizer::getEncounter(params)
    inter <- mizer::interaction_matrix(params)
    inter[pred, prey] <- inter[pred, prey] * lambda
    mizer::interaction_matrix(params) <- inter
    
    # Encounter must be restored BEFORE mortality:
    # Predation mortality carries a factor of (1 - f), which depends on encounter rate.
    # Restoring mortality first would bake spurious external mortality into prey species.
    mizer::ext_encounter(params) <- mizer::ext_encounter(params) +
        (enc_old - mizer::getEncounter(params))
    mizer::ext_mort(params) <- mizer::ext_mort(params) +
        (mort_old - mizer::getPredMort(params))
    params
}

#' Set Beverton-Holt Reproduction Level for a Specific Species
#'
#' @param params A MizerParams object.
#' @param s Character name of the species.
#' @param rl Numeric reproduction level in (0, 1).
#' @return An updated MizerParams object.
#' @export
set_rl <- function(params, s, rl) {
    mizer::setBevertonHolt(params, reproduction_level = stats::setNames(rl, s))
}

#' Find Maximum Coherent Reproduction Level (Ensuring erepro <= 1)
#'
#' Finds the reproduction level threshold above which required reproductive
#' efficiency erepro exceeds 1.0 (physical impossibility).
#'
#' @param params A MizerParams object.
#' @param s Character name of the species.
#' @param rl_hi Numeric upper search bound.
#' @param rl_lo Numeric lower search bound.
#' @return Capped upper bound for reproduction level.
#' @export
cap_rl <- function(params, s, rl_hi = 0.9995, rl_lo = 0.005) {
    er <- function(rl) {
        q <- suppressWarnings(set_rl(params, s, rl))
        mizer::species_params(q)$erepro[mizer::species_params(q)$species == s]
    }
    if (er(rl_hi) <= 1) return(rl_hi)
    u <- stats::uniroot(function(u) log(er(1 - exp(-u))),
                        c(-log(1 - rl_lo), -log(1 - rl_hi)))$root
    1 - exp(-0.95 * u)
}

#' Compute Yield Curve Over Multiples of Target F
#'
#' @param params A MizerParams object.
#' @param s Character name of the species.
#' @param F_t Target fishing mortality.
#' @param mult Multipliers around target F (default log-spaced grid from 0.15 to 2.6).
#' @param gear Name of the commercial gear (default "commercial").
#' @return Data frame with columns F and yield.
#' @export
yield_curve <- function(params, s, F_t,
                        mult = exp(seq(log(0.15), log(2.6), length.out = 8)),
                        gear = "commercial") {
    y <- suppressMessages(
        mizerExperimental::getYieldVsF(params, s, gear = gear, F_range = F_t * mult)
    )
    y[order(y$F), ]
}

#' Calculate Peak Ratio with Quadratic Refinement in log(F)
#'
#' Safely locates the global maximum of a yield curve even if bimodal.
#'
#' @param y Data frame from `yield_curve()`.
#' @param F_t Target fishing mortality.
#' @return Named vector with `ratio` (F_peak / F_target) and `edge` flag.
#' @export
peak_ratio <- function(y, F_t) {
    i <- which.max(y$yield)
    r <- y$F[i] / F_t
    edge <- i == 1 || i == nrow(y)
    if (!edge) {
        x <- log(y$F[(i - 1):(i + 1)])
        v <- y$yield[(i - 1):(i + 1)]
        den <- v[1] - 2 * v[2] + v[3]
        if (den != 0) {
            r <- exp(x[2] + (v[1] - v[3]) / (2 * den) * (x[2] - x[1])) / F_t
        }
    }
    c(ratio = r, edge = as.numeric(edge))
}

#' Bisect Reproduction Level to Match F_MSY Target
#'
#' Uses global curve scans and bisection in u = -log(1 - rl) to match F_peak to F_t.
#'
#' @param params A MizerParams object.
#' @param s Species name.
#' @param F_t Target fishing mortality (F_MSY).
#' @param rl_lo Minimum reproduction level (default 0.005).
#' @param rl_hi Maximum reproduction level (default 0.9995).
#' @param steps Maximum bisection iterations.
#' @param tol Relative error tolerance on peak ratio.
#' @return List with `rl`, `ratio`, and `status`.
#' @export
tune_peak <- function(params, s, F_t, rl_lo = 0.005, rl_hi = 0.9995,
                      steps = 5, tol = 0.10) {
    u2rl <- function(u) 1 - exp(-u)
    rl_hi <- cap_rl(params, s, rl_hi, rl_lo)
    u_lo <- -log(1 - rl_lo)
    u_hi <- -log(1 - rl_hi)
    
    calc_r <- function(rl_val) {
        y <- yield_curve(set_rl(params, s, rl_val), s, F_t)
        peak_ratio(y, F_t)[["ratio"]]
    }
    
    r_lo <- calc_r(rl_lo)
    if (r_lo >= 1) {
        return(list(rl = rl_lo, ratio = r_lo, status = "at lower bound"))
    }
    r_hi <- calc_r(rl_hi)
    if (r_hi <= 1) {
        return(list(rl = rl_hi, ratio = r_hi, status = "at upper bound"))
    }
    
    for (i in seq_len(steps)) {
        u <- (u_lo + u_hi) / 2
        r <- calc_r(u2rl(u))
        if (abs(log(r)) < log(1 + tol)) {
            return(list(rl = u2rl(u), ratio = r, status = "converged"))
        }
        if (r < 1) u_lo <- u else u_hi <- u
    }
    u <- (u_lo + u_hi) / 2
    list(rl = u2rl(u), ratio = calc_r(u2rl(u)), status = "bisected")
}

#' Extract Target Fishing Mortalities for All Species
#'
#' Uses FMSY from species_params, falling back to commercial catchability.
#'
#' @param params A MizerParams object.
#' @return Named numeric vector of target F per species.
#' @export
F_targets <- function(params) {
    gpf <- mizer::gear_params(params)
    gpf <- gpf[gpf$gear == "commercial", ]
    Ft <- stats::setNames(mizer::species_params(params)$FMSY,
                          mizer::species_params(params)$species)
    Fc <- stats::setNames(gpf$catchability, gpf$species)
    Ft[is.na(Ft)] <- Fc[is.na(Ft)]
    Ft
}

#' Compute Excess Multinomial Negative Log-Likelihood per Fish
#'
#' Measures fit of model size distributions against raw catch length bins.
#'
#' @param params A MizerParams object.
#' @param catch Catch data frame with columns `species`, `gear`, `length`, `dl`, `catch`.
#' @return Data frame of excess negative log-likelihood per species and gear.
#' @export
catch_nll <- function(params, catch) {
    sp <- mizer::species_params(params)
    fm <- mizer::getFMortGear(params)
    out <- NULL
    for (i in seq_len(nrow(sp))) {
        s <- sp$species[i]
        w <- mizer::w(params)
        dw <- mizer::dw(params)
        for (g in unique(catch$gear[catch$species == s])) {
            obs <- catch[catch$species == s & catch$gear == g & !is.na(catch$catch), ]
            if (nrow(obs) == 0) next
            cw <- fm[g, i, ] * mizer::initialN(params)[i, ]
            edges <- sp$a[i] * c(obs$length, max(obs$length + obs$dl))^sp$b[i]
            cum <- c(0, cumsum(cw * dw))
            lo <- w[1] - dw[1] / 2
            hi <- max(w + dw / 2)
            F_at <- stats::approx(c(lo, w + dw / 2), cum,
                                  xout = pmin(pmax(edges, lo), hi), rule = 2)$y
            prob <- diff(F_at)
            prob <- prob / sum(prob)
            k <- obs$catch
            phat <- k / sum(k)
            out <- rbind(out, data.frame(
                species = s, gear = g, n = sum(k),
                excess = (-sum(k * log(pmax(prob, 1e-12))) +
                              sum(k[k > 0] * log(phat[k > 0]))) / sum(k)
            ))
        }
    }
    out
}
