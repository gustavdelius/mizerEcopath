#' Simulate age-at-length data
#'
#' This function compares the mean age-at-length predicted by the model to
#' observed survey data for a given species. It aggregates observed
#' age-at-length data into length bins, simulates a cohort using the model, and
#' plots the mean age-at-length with confidence intervals for the observed data.
#'
#' @param params A MizerParams object.
#' @param species Name of the species to plot.
#' @param observed_df A data frame with columns `survey_date`, `Length`, `K`,
#'   and `count`
#' @param dt Time step for the model simulation (years). Default is 0.05.
#'
#' @return A data frame with
#' @export
simulateAge <- function(params, species, observed_df, dt = 0.05) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # highest survey age
    K_max <- max(observed_df$K)
    # survey length‑bin edges [cm]
    l_survey <- min(observed_df$Length):(max(observed_df$Length) + 1)
    s <- length(l_survey) - 1

    # Evolve a cohort over time in model ----

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
    t_max <- K_max + 2
    age   <- seq(dt, t_max,  by = dt)
    nsteps <- length(age)

    n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

    ## Convert to numbers from densities ----
    G <- n_hist[-1, 1:m]
    G <- sweep(G, 2, dw(params)[1:m], "*")

    # re-bin from survey bins to mizer bins
    B <- length_rebinning_matrix(l_model = l, l_survey = l_survey)
    G <- G %*% B

    # Split the data frame by unique survey date
    surveys <- split(observed_df, observed_df$survey_date)
    simulated_surveys <- list()

    # Loop through each survey to generate simulated data ----
    for (survey_date_str in names(surveys)) {
        survey_date_current <- as.numeric(survey_date_str)

        # Get the observations for this specific survey
        current_obs_df <- surveys[[survey_date_str]]

        # 1. Generate the model's predictions (proportions) for this specific date
        P_model <- generate_model_predictions_for_date(
            survey_date_current, G, a = age, l = l_survey[1:s],
            mu = sps$spawning_mu, kappa = sps$spawning_kappa,
            t_r = sps$t_r, a_min = sps$a_min
        )

        # 2. Simulate a sample from the model that mimics the real sampling effort
        simulated_K_values <- simulate_sample_from_model(P_model, current_obs_df)

        # 3. Store the results
        sim_df <- current_obs_df
        sim_df$K <- simulated_K_values
        simulated_surveys[[survey_date_str]] <- sim_df
    }

    # Combine the lists of data frames back into single data frames
    do.call(rbind, simulated_surveys)
}

#' Get log likelihood of age observations
#'
#' @param params A MizerParams object.
#' @param species Name of the species.
#' @param observed_df A data frame with columns `survey_date`, `Length`, `K`,
#'   and `count`
#' @param dt Time step for the model simulation (years). Default is 0.05.
#'
#' @return A data frame with
#' @export
getLogLik <- function(params, species, observed_df, dt = 0.05) {
    # TODO: remove code duplication from previous function
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # highest survey age
    K_max <- max(observed_df$K)
    # survey length‑bin edges [cm]
    l_survey <- min(observed_df$Length):(max(observed_df$Length) + 1)
    s <- length(l_survey) - 1

    # Evolve a cohort over time in model ----

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
    t_max <- K_max + 2
    age   <- seq(dt, t_max,  by = dt)
    nsteps <- length(age)

    n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

    ## Convert to numbers from densities ----
    G <- n_hist[-1, 1:m]
    G <- sweep(G, 2, dw(params)[1:m], "*")

    # re-bin from survey bins to mizer bins
    B <- length_rebinning_matrix(l_model = l, l_survey = l_survey)
    G <- G %*% B

    # Split the data frame by unique survey date
    surveys <- split(observed_df, observed_df$survey_date)

    calculate_and_aggregate_likelihood(
        surveys, G, a = age, l = l_survey[1:s],
        mu = sps$spawning_mu, kappa = sps$spawning_kappa,
        t_r = sps$t_r, a_min = sps$a_min)
}
