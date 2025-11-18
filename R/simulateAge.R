#' Simulate cohort evolution and calculate matrix G
#'
#' This function simulates the evolution of a cohort over time and calculates
#' the matrix G which represents the cohort's size distribution over time.
#'
#' @param params A MizerParams object.
#' @param species Name of the species to simulate.
#' @param l_survey Survey length‑bin edges (cm)
#' @param t_max Maximum time for simulation (years).
#' @param dt Time step for the model simulation (years). Default is 0.05.
#'
#' @return A list containing:
#'   - G: Matrix representing cohort size distribution over time
#'   - age: Age grid [years]
#' @export
simulateCohort <- function(params, species, l_survey, t_max, dt = 0.05) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    ## Set up initial pulse ----
    sps <- species_params(params)[species, ]

    # model length‑bin edges [cm]
    w <- w(params)
    l <- (w/sps$a)^(1/sps$b)
    m <- length(l) - 1  # number of model bins

    # Set initial condition
    n_init <- rep(0, length(l))
    n_init[2] <- 1  # Point source at small size
    initialN(params)[species, ] <- n_init

    ## Solve the PDE ----

    # fine‑age grid (centres) and widths  [years]
    age   <- seq(dt, t_max,  by = dt)
    nsteps <- length(age)

    n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

    ## Convert to numbers from densities ----
    G <- n_hist[-1, 1:m]
    G <- sweep(G, 2, dw(params)[1:m], "*")

    # re-bin from survey bins to mizer bins ----
    B <- length_rebinning_matrix(l_model = l, l_survey = l_survey)
    G <- G %*% B

    return(list(
        G = G,
        age = age
    ))
}

#' Get log likelihood of age observations
#'
#' @param params A MizerParams object.
#' @param species Name of the species.
#' @param surveys A data frame with survey age-at-length observations with
#'   columns `survey_date`, `Length`, `K`, and `count`.
#' @param dt Time step for the model simulation (years). Default is 0.05.
#'
#' @return A data frame with
#' @export
getLogLik <- function(params, species, surveys, dt = 0.05) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)
    sps <- params@species_params[species, ]

    # survey length‑bin edges [cm]
    l_survey <- min(surveys$Length):(max(surveys$Length) + 1)
    s <- length(l_survey) - 1

    # Evolve a cohort over time in model
    t_max <- max(surveys$K) + 2
    cohort_result <- simulateCohort(params, species, l_survey, t_max, dt)
    G <- cohort_result$G
    age <- cohort_result$age

    calculate_and_aggregate_likelihood(
        surveys, G, a = age, l = l_survey[1:s],
        mu = sps$spawning_mu, kappa = sps$spawning_kappa,
        annuli_date = sps$annuli_date, annuli_min_age = sps$annuli_min_age)
}
