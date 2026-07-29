# Choose the reproduction level of each species so that the maximum of its
# yield-vs-F curve sits close to the fishing mortality the species currently
# experiences. A higher reproduction level means more density dependence in
# reproduction, which makes the species less sensitive to fishing and hence
# moves F_MSY up.
#
# Doing this by eye with plotYieldVsF() for nine species is slow, so below we
# automate it: we only need to know on which side of the peak the current F
# lies, which takes two projections instead of a whole curve, and then we
# bisect on the reproduction level.

library(mizer)
library(mizerExperimental)

# Relative difference between the yield a little above and a little below the
# current F. Positive means the peak is above F_cur, negative that it is below.
side <- function(params, sp, F_cur, d = 1.25) {
    y <- suppressMessages(getYieldVsF(params, sp, gear = "commercial",
                                      F_range = c(F_cur / d, F_cur * d)))
    lo <- y$yield[which.min(abs(y$F - F_cur / d))]
    hi <- y$yield[which.min(abs(y$F - F_cur * d))]
    (hi - lo) / max(hi, lo, 1e-30)
}

set_rl <- function(params, sp, rl) {
    setBevertonHolt(params, reproduction_level = setNames(rl, sp))
}

# Bisect the reproduction level until the peak sits at F_cur. We bisect in
# u = -log(1 - rl) because F_MSY is most sensitive to the reproduction level
# as it approaches 1.
tune <- function(params, sp, F_cur, rl_lo = 0.01, rl_hi = 0.99,
                 steps = 6, tol = 0.02) {
    u2rl <- function(u) 1 - exp(-u)
    u_lo <- -log(1 - rl_lo)
    u_hi <- -log(1 - rl_hi)
    s_lo <- side(set_rl(params, sp, rl_lo), sp, F_cur)
    if (s_lo >= 0) return(list(rl = rl_lo, side = s_lo, status = "at lower bound"))
    s_hi <- side(set_rl(params, sp, rl_hi), sp, F_cur)
    if (s_hi <= 0) return(list(rl = rl_hi, side = s_hi, status = "at upper bound"))
    for (i in seq_len(steps)) {
        u <- (u_lo + u_hi) / 2
        s <- side(set_rl(params, sp, u2rl(u)), sp, F_cur)
        if (abs(s) < tol) return(list(rl = u2rl(u), side = s, status = "converged"))
        if (s < 0) u_lo <- u else u_hi <- u
    }
    u <- (u_lo + u_hi) / 2
    list(rl = u2rl(u), side = side(set_rl(params, sp, u2rl(u)), sp, F_cur),
         status = "bisected")
}

# Location of the peak, for checking the result. The quadratic refinement is
# in log F because the grid is multiplicative.
peak <- function(params, sp, F_cur,
                 mult = seq(1/4, 1.5, 1/8)) {
    y <- suppressMessages(getYieldVsF(params, sp,  gear = "commercial",,
                                      F_range = F_cur * mult))
    y <- y[order(y$F), ]
    pl <- ggplot(y, aes(x = F, y = yield)) +
        geom_line() +
        xlab("Fishing mortality (1/yr)") +
        ylab("Yield") +
        ggtitle(sp) +
        geom_vline(xintercept = F_cur, linetype = "dashed", colour = "grey")
    print(pl)
    i <- which.max(y$yield)
    if (i > 1 && i < nrow(y)) {
        x <- log(y$F[(i - 1):(i + 1)])
        v <- y$yield[(i - 1):(i + 1)]
        d <- (v[1] - v[3]) / (2 * (v[1] - 2 * v[2] + v[3]))
        pk <- exp(x[2] + d * (x[2] - x[1]))
    } else {
        pk <- y$F[i]  # peak outside the range
    }
    c(F_peak = pk, ratio = pk / F_cur)
}


