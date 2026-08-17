# Shared smoother for noisy quantities estimated by length class.
#
# Used by `catchability_Richard.R` (fishing mortality) and
# `survey_selectivity_Richard.R` (survey selectivity), which differ only in the
# link function and in which part of the model they write to.

#' Smooth a noisy quantity given by length class
#'
#' Fits a weighted local linear regression to `link(y)` against length. The
#' width of the local window is set by two constraints, whichever is tighter:
#'
#' * a bandwidth measured in *fish* rather than in cm. The scatter is driven by
#'   sample size, and the sparse length classes all sit at the large-length end,
#'   so a window holding a fixed number of fish is narrow through the
#'   well-sampled rise and wide over the noisy tail. A window of fixed length
#'   cannot do both.
#' * a cap in cm. Where a species is sparsely sampled at *small* lengths (Cod
#'   has only 20 fish below 26 cm) the fish-count window would otherwise reach
#'   tens of cm up into the plateau and drag the fit at the lower boundary up
#'   by more than an order of magnitude.
#'
#' @param length_cm Lengths (cm) at which `y` was estimated.
#' @param y The quantity to smooth, strictly positive (and strictly below 1 for
#'   the logit link).
#' @param n Number of fish behind each estimate. Sets both the regression
#'   weights and the cumulative-fish axis on which the window is measured.
#' @param link,inv_link The scale on which to smooth. `log`/`exp` for a positive
#'   quantity such as a rate, `qlogis`/`plogis` for a proportion.
#' @param bw_fish Bandwidth in number of fish. Capped at `sum(n) / 5` so that
#'   species with few fish in total (Cod, Herring) do not collapse to a single
#'   global straight line.
#' @param bw_cm Cap on the bandwidth in cm.
#' @param fall_off Factor by which the quantity drops per cm below the smallest
#'   observed length. The data start at 20 cm, but fish smaller than that are
#'   not retained by the net at all, so below the data the curve must fall off
#'   much faster than a continuation of the fitted trend would suggest.
#' @return A function giving the smoothed quantity at any length. Outside the
#'   range of the data it drops off steeply below the smallest observed length
#'   and is held constant above the largest observed length. Holding it constant
#'   freezes a curve that is usually still descending there, which is deliberate:
#'   the decline at large sizes is real and must not be smoothed back up.
smooth_at_length <- function(length_cm, y, n, link = log, inv_link = exp,
                             bw_fish = 2000, bw_cm = 6, fall_off = 10) {
    o <- order(length_cm)
    length_cm <- length_cm[o]
    y <- y[o]
    n <- n[o]
    ly <- link(y)

    # Position of each length class on the cumulative-fish axis
    u <- cumsum(n) - n / 2
    bw_fish <- min(bw_fish, sum(n) / 5)

    # Weighted local linear fit at one length, returning the fitted value on the
    # link scale and the local slope with respect to length.
    local_fit <- function(l) {
        u0 <- approx(length_cm, u, xout = l, rule = 2)$y
        k <- n * dnorm(u - u0, sd = bw_fish) * dnorm(length_cm - l, sd = bw_cm)
        co <- lm.wfit(cbind(1, length_cm - l), ly, k)$coefficients
        # The slope is unidentified if the window collapses onto one class
        c(value = co[[1]], slope = if (is.na(co[[2]])) 0 else co[[2]])
    }

    l_min <- min(length_cm)
    l_max <- max(length_cm)
    edge_lo <- local_fit(l_min)
    edge_hi <- local_fit(l_max)
    # Below the data, use the steeper of the fitted slope at the lower edge and
    # the imposed net-retention fall-off.
    slope_lo <- max(edge_lo[["slope"]], log(fall_off))

    function(l) {
        out <- numeric(length(l))
        below <- l < l_min
        above <- l > l_max
        inside <- !below & !above
        out[inside] <- vapply(l[inside], function(x) local_fit(x)[["value"]],
                              numeric(1))
        out[below] <- edge_lo[["value"]] + slope_lo * (l[below] - l_min)
        out[above] <- edge_hi[["value"]]
        inv_link(out)
    }
}
