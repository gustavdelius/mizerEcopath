#' Solve steady-state equation including diffusion
#'
#' This function solves the steady-state equation for a single species
#' including diffusion. It uses an upwind difference scheme to solve the
#' second-order ODE.
#'
#' For more details on the upwind difference scheme, see the vignette
#' \code{vignette("diffusion")}.
#'
#' @param params A MizerParams object
#' @param species The species for which to solve the steady-state equation.
#' @return A vector of steady-state abundances for the specified species
#' @export
steady_diffusion <- function(params, species) {
    params <- validParams(params)
    species <- valid_species_arg(params, species = species,
                                 error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Currently this species deals with a single species at a time.")
    }
    sps <- species_params(params)[species, ]
    n <- sps$n
    d_over_g <- sps$d_over_g
    w <- w(params)
    x <- log(w / w[1])
    h <- x[2] - x[1]
    N <- length(w) - 2 # Number of interior points
    # TODO: select only the relevant sizes from w_min to w_max

    # Get mortality and growth rates
    mu <- getMort(params)[species, ]
    g <- getEGrowth(params)[species, ]

    # Calculate diffusion rate as a power law
    g_0 <- g[1] / w[1]^n
    d_0 <- d_over_g * g_0
    d <- d_0 * w^(n + 1)
    # Transform to logarithmic space
    dtilde <- d / w^2
    dtilde_prime <- d_0 * (n - 1) * w^(n-1)
    gtilde <- g / w - 0.5 * dtilde
    # Transform to standard form for diffusion term
    ghat <- gtilde - dtilde_prime / 2

    # Set initial abundance at smallest size
    n0 <- initialN(params)[species, 1] * w[1]

    # Solve for steady state using the upwind scheme
    n_steady <- solve_diffusion_ode(dtilde, ghat, mu, n0, h) / w

    return(n_steady)
}

#' Solve 2nd order ODE
#'
#' Solve ODE of the form
#' \deqn{(d(x) n'(x))'-(g(x)n(x))'-\mu(x)n(x)=0}
#' with boundary conditions \eqn{n(0)=n_0} and \eqn{n(N+1)=0}
#' using an upwind difference scheme.
#'
#' For more details on the upwind difference scheme, see the vignette
#' \code{vignette("diffusion")}. This function is used by [steady_diffusion()].
#'
#' @param d A vector of diffusion coefficients \eqn{d_i} for each grid point
#' @param g A vector of advection coefficients \eqn{g_i} for each grid point
#' @param mu A vector of mortality coefficients \eqn{\mu_i} for each grid point
#' @param n0 The value of \eqn{n} at the left boundary
#' @param h The grid spacing
#' @return A vector of solutions \eqn{n_i} for each grid point
#' @examples
#' d <- c(0, 1, 1, 0)
#' g <- c(0, 1, -1, 0)
#' mu <- c(0, 0.1, 0.1, 0)
#' n0 <- 1
#' h <- 1
#' solve_diffusion_steady(d, g, mu, n0, h)
#' @export
solve_diffusion_ode <- function(d, g, mu, n0, h) {
    # Number of interior points (N in the maths, but R indices run 1 to N+2)
    N <- length(d) - 2

    # Calculate midpoint values for diffusion
    # \tilde{d}_{i+1/2} = (\tilde{d}_i+\tilde{d}_{i+1})/2 for i = 0,...,N
    # In R:
    # d_half[1] = \tilde{d}_{1/2}, ..., d_half[N+1] = \tilde{d}_{N+1/2}
    d_half <- (d[1:(N+1)] + d[2:(N+2)]) / 2

    # Prepare diagonals for upwind scheme
    abs_g <- abs(g)
    g_plus <- (abs_g + g) / 2  # g^+
    g_minus <- (abs_g - g) / 2 # g^-

    # Lower diagonal L_i = (d_{i-1/2}/2) + h * g^+_{i-1}
    # In R: L[i] = (d_half[i]/2) + h * g_plus[i],
    # corresponds to L_i for i = 1,...,N
    L <- (d_half[1:N] / 2) + h * g_plus[1:N]
    # Upper diagonal U_i = (d_{i+1/2}/2) + h * g^-_{i+1}
    # In R: U[i] = (d_half[i+1]/2) + h * g_minus[i+1],
    # corresponds to U_i for i = 1,...,N
    U <- (d_half[2:(N+1)] / 2) + h * g_minus[3:(N+2)]
    # Main diagonal (D): -((d_{i+1/2} + d_{i-1/2})/2 + h*|g_i| + h^2*mu_i)
    # In R: D[i] = -((d_half[i+1] + d_half[i])/2 + h * abs_g[i] + h^2 * mu[i]),
    # corresponds to D_i for i = 1,...,N
    D <- -((d_half[2:(N+1)] + d_half[1:N]) / 2 +
               h * abs_g[2:(N+1)] + h^2 * mu[2:(N+1)])

    # Check stability conditions
    if (any(is.na(U)) || any(is.na(L)) || any(is.na(D))) {
        stop("NA values detected in diagonals. Check input lengths and indexing.")
    }
    if (any(U <= 0) || any(L <= 0) || any(D >= 0)) {
        warning("Stability conditions not met.")
    }
    if (any(U + L > -D)) {
        warning("Stability condition U + L <= -D not met.")
    }

    # Right-hand side: n0 at left boundary, zeros elsewhere
    b <- numeric(N)
    b[1] <- -L[1] * n0
    # The right boundary is handled by the system structure (last value is zero)

    # Solve using double sweep method (Thomas algorithm)
    n <- solve_double_sweep(U, L, D, b)

    # Add boundary values: n0 at left, 0 at right
    solution <- numeric(N + 2)
    solution[1] <- n0
    solution[2:(N + 1)] <- n
    solution[N + 2] <- 0
    return(solution)
}
