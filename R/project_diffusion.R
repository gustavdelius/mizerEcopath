#' Project forward in time with diffusion
#'
#' This function solves the PDE for a single species
#' including diffusion. It uses an upwind difference scheme.
#'
#' For more details  see the vignette
#' \code{vignette("diffusion")}.
#'
#' @param params A MizerParams object
#' @param species The species for which to solve the steady-state equation.
#' @param dt Time step size.
#' @param nsteps Number of time steps to compute.
#' @return A matrix (time x size) with the solution of the PDE.
#' @export
project_diffusion <- function(params, species, dt = 0.05, nsteps = 200) {
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
    N <- length(w) # Number of grid points
    # TODO: select only the relevant sizes from w_min to w_max

    # Get mortality and growth rates
    mu <- getMort(params)[species, ]
    g <- getEGrowth(params)[species, ]
    emigration <- emigration(params)[species, ]

    # Calculate diffusion rate as a power law
    g_0 <- g[1] / w[1]^n
    d_0 <- d_over_g * g_0
    d <- d_0 * w^(n + 1)
    # Transform to logarithmic space
    dtilde <- d / w^2
    dtilde_prime <- d_0 * (n - 1) * w^(n-1)
    gtilde <- g / w - 0.5 * dtilde
    etilde <- emigration * w
    # Transform to standard form for diffusion term
    ghat <- gtilde - dtilde_prime / 2

    # Set initial abundance at smallest size
    n_init <- initialN(params)[species, ]

    # Solve the PDE
    n_hist <- solve_diffusion_pde(dtilde, ghat, mu, etilde,
                                  n_init, h, dt, nsteps)
    # Convert to density in w
    n_hist <- sweep(n_hist, 2, w, "/")

    return(n_hist)
}

#' Function to solve a reaction-diffusion-advection PDE
#'
#' This function uses the implicit Euler method to solve a PDE of the form:
#' \deqn{\dot{n}(x,t) = (d(x) n'(x,t))'-(g(x)n(x,t))'-\mu(x)n(x,t)}
#' @param d A vector of diffusion rates.
#' @param g A vector of advection rates.
#' @param mu A vector of reaction rates.
#' @param emigration A vector of emigration rates.
#' @param n_init Initial condition for the solution.
#' @param h Spatial step size.
#' @param dt Time step size.
#' @param nsteps Number of time steps to compute.
#' @return A matrix (time x size) with the solution of the PDE.
#' @export
solve_diffusion_pde <- function(d, g, mu, emigration, n_init, h, dt, nsteps) {
    # Because of the way we calculate the derivative of the diffusion we
    # we can solve the equation only at interior points
    N <- length(d) - 2
    interior <- 2:(N + 1)
    d_half <- (d[1:(N+1)] + d[2:(N+2)]) / 2
    abs_g <- abs(g)
    g_plus <- (abs_g + g) / 2
    g_minus <- (abs_g - g) / 2
    L <- (d_half[1:N] / 2) + h * g_plus[1:N]
    U <- (d_half[2:(N+1)] / 2) + h * g_minus[3:(N+2)]
    D <- -((d_half[2:(N+1)] + d_half[1:N]) / 2 +
               h * abs_g[2:(N+1)] + h^2 * mu[2:(N+1)])

    # For implicit Euler: (I - dt*A) n^{k+1} = n^k
    D_new <- 1 - dt / h^2 * D
    L_new <- -dt / h^2 * L
    U_new <- -dt / h^2 * U

    n <- n_init[interior]
    n_hist <- matrix(0, nrow = nsteps + 1, ncol = length(d))
    n_hist[1, interior] <- n

    b <- n - dt * emigration[interior]  # Add emigration term

    for (step in 1:nsteps) {
        # Solve (I - dt*A) n_new = n_old
        n_new <- solve_double_sweep(U_new, L_new, D_new, b)
        n_new[length(n_new)] <- 0  # Enforce boundary at large size
        n <- n_new
        n_hist[step + 1, interior] <- n
    }
    return(n_hist)
}
