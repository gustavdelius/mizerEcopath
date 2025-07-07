#' Set initial abundances to single-species steady state abundances
#'
#' `r lifecycle::badge("experimental")`
#' This first calculates growth and death rates that arise from the current
#' initial abundances. Then it uses these growth and death rates to
#' determine the steady-state abundances of the selected species.
#'
#' The result of applying this function is of course not a multi-species steady
#' state, because after changing the abundances of the selected species the
#' growth and death rates will have changed.
#'
#' @param params A MizerParams object
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param d_over_g The ratio of the diffusion rate to the growth rate at the
#'   size of offspring.
#' @param keep A string determining which quantity is to be kept constant. The
#'   choices are "egg" which keeps the egg density constant, "biomass" which
#'   keeps the total biomass of the species constant and "number" which keeps
#'   the total number of individuals constant.
#' @return A MizerParams object in which the initial abundances of the selected
#'   species are changed to their single-species steady state abundances.
#' @export
steadySingleSpeciesDiffusion <-
    function(params, species = NULL, d_over_g = 0.15,
             keep = c("egg", "biomass", "number")) {
    species <- valid_species_arg(params, species)
    keep <- match.arg(keep)

    biomass <- getBiomass(params)
    number <- getN(params)

    # Use growth and mortality from current abundances
    growth_all <- getEGrowth(params)
    mort_all <- getMort(params)

    # Loop through all species and calculate their steady state abundances
    # using the current growth and mortality rates
    for (sp in species) {
        growth <- growth_all[sp, ]
        mort <- mort_all[sp, ]

        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
        idx <- w_min_idx:(w_max_idx - 1)

        # Check that species can grow to maturity at least
        w_mat_idx <- sum(params@w <= params@species_params[sp, "w_mat"])
        if (any(growth[w_min_idx:w_mat_idx] == 0)) {
            stop(sp, " cannot grow to maturity")
        }

        # Keep egg density constant
        n0 <- params@initial_n[sp, w_min_idx]
        # Steady state solution of the upwind-difference scheme used in project
        solve_ode <-
            solve_ode_steady_state(growth[idx], mort[idx], n0, d_over_g)
        params@initial_n[sp, ] <- 0
        params@initial_n[sp, w_min_idx:w_max_idx] <-
            n0 * c(1, cumprod(growth[idx] /
                                  ((growth + mort * params@dw)[idx + 1])))
    }

    if (any(is.infinite(params@initial_n))) {
        stop("Candidate steady state holds infinities")
    }
    if (any(is.na(params@initial_n) | is.nan(params@initial_n))) {
        stop("Candidate steady state holds non-numeric values")
    }

    if (keep == "biomass") {
        factor <- biomass / getBiomass(params)
        params@initial_n <- params@initial_n * factor
    }
    if (keep == "number") {
        factor <- number / getN(params)
        params@initial_n <- params@initial_n * factor
    }

    params@time_modified <- lubridate::now()
    params
    }

# Helper function to solve steady state ODE
solve_ode_steady_state <- function(growth, mort, d_over_g, n0, w, N0, h) {
    N <- length(growth) - 1  # Number of internal points
    if (length(mort) != N + 1 || length(growth) != N + 1) {
        stop("Growth and mortality vectors must have the same length.")
    }
    g_0 <- growth[1] / w[1]^n
    d_0 <- d_over_g * g_0
    d <- d_0 * w^(n + 1)
    dtilde <- d / w^2
    gtilde <- g / w - 0.5 * dtilde

    # Coefficients for the tridiagonal matrix
    U <- (dtilde / 2)[3:(N+2)]    # Upper diagonal
    L <- (dtilde / 2 + h * gtilde)[1:N]     # Lower diagonal
    D <- (-dtilde - h * gtilde - h^2 * mort)[2:(N+1)]  # Main diagonal

    # Solve using double sweep method
    return(solve_ode_double_sweep(U, L, D, N0))
}

# Helper function for double sweep method
solve_ode_double_sweep <- function(U, L, D, n0) {
    N <- length(D)  # Number of internal points
    if (length(U) != N || length(L) != N) {
        stop("U, L, and D must have the same length.")
    }
    # Initialize arrays for alpha and beta coefficients
    alpha <- numeric(N + 1)
    beta <- numeric(N + 1)

    # Initialize n array for solution
    n <- numeric(N + 1)
    n[N + 1] <- 0  # Boundary condition at x_max

    # Initial condition
    alpha[1] <- 0
    beta[1] <- n0

    # Forward sweep - calculate alpha and beta coefficients
    for (i in 1:N) {
        denom <- alpha[i] * L[i] + D[i]
        alpha[i + 1] <- -U[i] / denom
        beta[i + 1] <- -(beta[i] * L[i]) / denom
    }

    # Work backwards to get interior n values
    for (i in N:1) {
        n[i] <- alpha[i + 1] * n[i + 1] + beta[i + 1]
    }

    # Add initial value
    solution <- numeric(N + 2)
    solution[1] <- n0  # Set the first value to n0
    solution[2:(N + 2)] <- n  # Fill in the other points
    return(solution)
}
