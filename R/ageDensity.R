# Code from https://chatgpt.com/share/68a34a22-6470-8007-8348-51c0ca02575f

# ===== Helper functions (circular arithmetic on [0,1)) =====
wrap01 <- function(x) x - floor(x)                      # wrap to [0,1)
circ_diff01 <- function(a, b) {                         # shortest signed diff a-b in (-0.5,0.5]
    d <- a - b
    d - round(d)
}

# ===== von Mises density on [0,1) via angle 2*pi*s =====
dvm01 <- function(s, mu, kappa) {
    ang <- 2*pi*wrap01(s - mu)
    (exp(kappa*cos(ang))) / (2*pi*I0(kappa))
}

# Modified Bessel I0 wrapper (base R has besselI)
I0 <- function(kappa) besselI(kappa, nu = 0, expon.scaled = FALSE)

# ===== Edge readability function v(δ): probability newest annulus is counted =====
v_edge <- function(delta, alpha, sigma) {
    if (sigma <= 0) return(as.numeric(circ_diff01(delta, alpha) >= 0))  # hard step
    plogis(circ_diff01(delta, alpha) / sigma)
}

# ===== P(K = k | R = r, δ) given survey phase theta and edge v(δ) =====
pK_given_r_delta <- function(k, r, delta, theta, alpha, sigma) {
    # Deterministic baseline ring count K* at survey:
    # After α in the survey year, the "new" annulus is present for everyone.
    Kstar <- if (theta >= alpha) r else max(r - 1, 0L)

    v <- v_edge(delta, alpha, sigma)
    # With probability v, we count one more than K* (new annulus read);
    # with 1-v, we read K*.
    p <- (1 - v) * as.numeric(k == Kstar) + v * as.numeric(k == (Kstar + 1L))

    # Clip impossible negatives (defensive)
    if (k < 0) p <- 0
    p
}

# ===== Main: density of A given K = k =====
# Returns a data.frame with columns: a (age), dens (pdf), plus posterior over R.
age_density_given_k <- function(
        k,
        theta,               # survey phase in [0,1)
        mu, kappa,           # von Mises spawning params on [0,1)
        alpha,               # annulus visibility phase in [0,1)
        sigma = 0.03,        # width of the edge transition; 0 => hard step
        R_max = max(10, k + 5),  # maximum completed calendar years to consider
        w_R   = NULL,        # optional length R_max+1 vector of prior weights over R=0..R_max
        n_delta = 2048,      # integration grid over δ
        n_per_year = 512     # resolution for the final A density within each [r, r+1)
) {
    stopifnot(theta >= 0 && theta < 1, alpha >= 0 && alpha < 1)
    R_vals <- 0:R_max
    if (is.null(w_R)) w_R <- rep(1, length(R_vals))
    if (length(w_R) != length(R_vals)) stop("w_R must have length R_max+1.")
    w_R <- pmax(w_R, 0); if (sum(w_R) == 0) stop("All w_R are zero.")
    w_R <- w_R / sum(w_R)

    # Grid for δ in [0,1)
    delta_grid <- seq(0, 1, length.out = n_delta + 1); delta_grid <- delta_grid[-length(delta_grid)]
    d_delta <- 1 / n_delta

    # f_Δ(δ) = f_S((θ - δ) mod 1) where S ~ von Mises(μ, κ) on [0,1)
    f_delta <- dvm01(wrap01(theta - delta_grid), mu = mu, kappa = kappa)
    # Normalize (numerical safety)
    f_delta <- f_delta / (sum(f_delta) * d_delta)

    # For each R, compute Pr(K=k | R)
    pK_given_R <- vapply(
        R_vals,
        function(r) {
            p_vec <- pK_given_r_delta(k, r, delta_grid, theta, alpha, sigma)
            sum(p_vec * f_delta) * d_delta
        },
        numeric(1)
    )

    # Posterior over R given K=k
    numer_R <- w_R * pK_given_R
    post_R <- if (sum(numer_R) > 0) numer_R / sum(numer_R) else rep(0, length(numer_R))

    # Now build the conditional density of A on a grid by mixing over R.
    # For each r, the support is a in [r, r+1) with density ∝ Pr(K=k | r, δ=a-r) f_Δ(δ)
    build_r_piece <- function(r, weight_r) {
        if (weight_r == 0) return(NULL)
        # Grid for a in [r, r+1)
        delta_r <- seq(0, 1, length.out = n_per_year + 1); delta_r <- delta_r[-length(delta_r)]
        # Interpolate f_delta at these delta_r points (use nearest-neighbor from precomputed grid)
        idx <- pmin(n_delta, pmax(1L, round(delta_r * n_delta)))
        fdel <- f_delta[idx]
        pKr <- pK_given_r_delta(k, r, delta_r, theta, alpha, sigma)
        g_unnorm <- pKr * fdel
        area <- sum(g_unnorm) * (1 / n_per_year)
        if (area == 0) return(NULL)
        g <- g_unnorm / area
        data.frame(a = r + delta_r, dens = weight_r * g, r = r)
    }

    pieces <- do.call(rbind, Map(build_r_piece, R_vals, post_R))
    if (is.null(pieces) || nrow(pieces) == 0) {
        warning("Degenerate result: no posterior mass; check parameters.")
        return(list(density = data.frame(a = numeric(0), dens = numeric(0)),
                    post_R = data.frame(R = R_vals, post = post_R),
                    pK_given_R = data.frame(R = R_vals, p = pK_given_R)))
    }

    # Normalize the overall mixture density (guard small numerical drift)
    total_area <- sum(pieces$dens) * (1 / n_per_year)
    pieces$dens <- pieces$dens / total_area

    list(
        density = pieces[, c("a", "dens")],
        post_R  = data.frame(R = R_vals, post = post_R),
        pK_given_R = data.frame(R = R_vals, p = pK_given_R)
    )
}
