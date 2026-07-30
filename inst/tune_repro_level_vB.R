# Choose the reproduction level of each species so that the maximum of its
# yield-vs-F curve sits at a target fishing mortality. A higher reproduction
# level means more density dependence in reproduction, which makes the species
# less sensitive to fishing and hence moves F_MSY up.
#
# The allometric build (`inst/tune_repro_level.R`) uses a cheap two-point test:
# compare the yield a little above and a little below the target, and the sign
# tells you which side of the peak the target is on. That is valid only for a
# unimodal curve. In the von Bertalanffy model several species have bimodal
# yield curves - Blue whiting's has a maximum at low F, a minimum near its FMSY
# and a second rise beyond - and there the sign of the local slope says nothing
# about where the global maximum is. We therefore evaluate a whole curve on a
# log grid of multiples of the target and locate its global maximum.

library(mizer)
library(mizerExperimental)
library(ggplot2)

MULT <- exp(seq(log(0.15), log(2.6), length.out = 8))

yield_curve <- function(params, s, F_t, mult = MULT) {
    y <- suppressMessages(getYieldVsF(params, s, gear = "commercial",
                                      F_range = F_t * mult))
    y[order(y$F), ]
}

# Global peak location as a multiple of the target, refined quadratically in
# log F because the grid is multiplicative. `edge` flags a maximum at an end of
# the grid, where the returned value is only a bound - but its position relative
# to 1 is still correct, which is all the bisection below needs.
peak_ratio <- function(y, F_t) {
    i <- which.max(y$yield)
    r <- y$F[i] / F_t
    edge <- i == 1 || i == nrow(y)
    if (!edge) {
        x <- log(y$F[(i - 1):(i + 1)]); v <- y$yield[(i - 1):(i + 1)]
        den <- v[1] - 2 * v[2] + v[3]
        if (den != 0) {
            r <- exp(x[2] + (v[1] - v[3]) / (2 * den) * (x[2] - x[1])) / F_t
        }
    }
    c(ratio = r, edge = as.numeric(edge))
}

gpeak <- function(params, s, F_t) peak_ratio(yield_curve(params, s, F_t), F_t)

set_rl <- function(params, s, rl) {
    setBevertonHolt(params, reproduction_level = setNames(rl, s))
}

# The largest reproduction level worth bracketing at. Above some level the
# reproductive efficiency `erepro` that mizer needs in order to sustain the
# calibrated recruitment exceeds 1, i.e. the model is being asked for more larval
# biomass than the spawning stock can physically produce, and a peak located
# there is meaningless. The crossing point is species-specific (0.96 for Herring,
# 0.99 for Megrim, up to 0.99997 for Cod) and always below the nominal upper
# bracket of 0.9995, so without this cap a species reported "at upper bound"
# would be returning a reproduction level from an incoherent region. No
# projection is involved, so this is cheap.
cap_rl <- function(params, s, rl_hi, rl_lo) {
    er <- function(rl) {
        q <- suppressWarnings(set_rl(params, s, rl))
        species_params(q)$erepro[species_params(q)$species == s]
    }
    if (er(rl_hi) <= 1) return(rl_hi)
    u <- stats::uniroot(function(u) log(er(1 - exp(-u))),
                        c(-log(1 - rl_lo), -log(1 - rl_hi)))$root
    1 - exp(-0.95 * u)
}

# Bisect the reproduction level on log(F_peak / F_target). More density
# dependence buffers reproduction against fishing and so moves the peak up, so
# the objective is increasing in the reproduction level. We bisect in
# u = -log(1 - rl) because the peak is most sensitive to rl as it approaches 1.
tune_peak <- function(params, s, F_t, rl_lo = 0.005, rl_hi = 0.9995,
                      steps = 5, tol = 0.10) {
    u2rl <- function(u) 1 - exp(-u)
    rl_hi <- cap_rl(params, s, rl_hi, rl_lo)
    u_lo <- -log(1 - rl_lo); u_hi <- -log(1 - rl_hi)
    r_lo <- gpeak(set_rl(params, s, rl_lo), s, F_t)
    if (r_lo[["ratio"]] >= 1) {
        return(list(rl = rl_lo, ratio = r_lo[["ratio"]],
                    status = "at lower bound"))
    }
    r_hi <- gpeak(set_rl(params, s, rl_hi), s, F_t)
    if (r_hi[["ratio"]] <= 1) {
        return(list(rl = rl_hi, ratio = r_hi[["ratio"]],
                    status = "at upper bound"))
    }
    for (i in seq_len(steps)) {
        u <- (u_lo + u_hi) / 2
        r <- gpeak(set_rl(params, s, u2rl(u)), s, F_t)[["ratio"]]
        if (abs(log(r)) < log(1 + tol)) {
            return(list(rl = u2rl(u), ratio = r, status = "converged"))
        }
        if (r < 1) u_lo <- u else u_hi <- u
    }
    u <- (u_lo + u_hi) / 2
    list(rl = u2rl(u),
         ratio = gpeak(set_rl(params, s, u2rl(u)), s, F_t)[["ratio"]],
         status = "bisected")
}

# The target F for each species: its FMSY where the data give one, otherwise the
# fishing mortality the model says it currently experiences. With effort = 1 the
# fully-selected F equals the fitted catchability.
F_targets <- function(params) {
    gpf <- gear_params(params)
    gpf <- gpf[gpf$gear == "commercial", ]
    Ft <- setNames(species_params(params)$FMSY, species_params(params)$species)
    Fc <- setNames(gpf$catchability, gpf$species)
    Ft[is.na(Ft)] <- Fc[is.na(Ft)]
    Ft
}

# Save a yield curve as a figure, for the record.
save_curve <- function(params, s, F_t, dir,
                       mult = exp(seq(log(0.15), log(2.6), length.out = 13))) {
    y <- yield_curve(params, s, F_t, mult)
    pk <- peak_ratio(y, F_t)
    pl <- ggplot(y, aes(x = F, y = yield)) +
        geom_line(linewidth = 0.8) +
        geom_vline(xintercept = F_t, linetype = "dashed", colour = "grey40") +
        labs(x = "Fishing mortality (1/yr)", y = "Yield",
             title = sprintf("%s  (peak at %.2f x target)", s,
                             pk[["ratio"]])) +
        theme_bw()
    ggsave(file.path(dir, paste0(gsub(" ", "_", s), ".png")), pl,
           width = 4, height = 3, dpi = 110)
    pk
}
