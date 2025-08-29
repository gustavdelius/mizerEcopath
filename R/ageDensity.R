#' von Mises Probability Density Function
#' A lightweight implementation used to model circular seasonality
#' (e.g., spawning day within a year) without adding dependencies.
#' @param x Angle in radians.
#' @param mu Mean direction in radians.
#' @param kappa Concentration parameter (higher means more concentrated).
#' @return The probability density evaluated at `x`.
#' @keywords internal
#' @examples
#' # Density at angle pi when centered at pi with moderate concentration
#' von_mises_pdf(pi, mu = pi, kappa = 2)
von_mises_pdf <- function(x, mu, kappa) {
    denominator <- 2 * pi * besselI(kappa, 0)
    numerator <- exp(kappa * cos(x - mu))
    return(numerator / denominator)
}

#' Spawning Density Function S(d)
#' Calculates relative spawning intensity for numeric dates using a
#' von Mises density on the unit circle of the year.
#' @param numeric_dates A vector of numeric dates (e.g., 2023.45). Only the
#'   fractional part is used to determine day-of-year.
#' @param mu The mean spawning date as a fraction of a year in \[0, 1). For
#'   example, 0.5 is mid-year.
#' @param kappa The concentration parameter for the spawning distribution.
#' @return A numeric vector of relative spawning intensities with the same
#'   length as `numeric_dates`.
#' @keywords internal
#' @examples
#' # Peak at mid-year with moderate spread
#' spawning_density(numeric_dates = c(2020.45, 2020.50, 2020.55),
#'                  mu = 0.5, kappa = 4)
spawning_density <- function(numeric_dates, mu, kappa) {
    day_fraction <- numeric_dates %% 1
    day_rad <- day_fraction * 2 * pi
    mu_rad <- mu * 2 * pi
    density <- von_mises_pdf(day_rad, mu = mu_rad, kappa = kappa)
    return(density)
}


#' Age-to-Ring Mapping Function Calculate_K(a)
#' Deterministically maps true age to the expected annuli count, given
#' a survey date, an annual ring-formation day, and a minimum age threshold.
#' @param age_in_years A numeric vector of true ages in years.
#' @param survey_date The survey date as a numeric year (e.g., 2023.25).
#' @param annuli_date Ring formation day as fraction of a year in \[0, 1).
#' @param annuli_min_age The minimum age (years) required to form the first ring.
#' @return An integer vector with the calculated number of rings (K) for each age.
#' @keywords internal
#' @examples
#' calculate_K(age_in_years = c(0.3, 1.2, 2.7), survey_date = 2023.5,
#'             annuli_date = 0.25, annuli_min_age = 0.5)
calculate_K <- function(age_in_years, survey_date, annuli_date, annuli_min_age) {
    sapply(age_in_years, function(age) {
        if (age < 0) return(0)
        birth_date <- survey_date - age
        birth_year <- floor(birth_date)
        next_ring_date <- birth_year + annuli_date
        if (next_ring_date <= birth_date) {
            next_ring_date <- (birth_year + 1) + annuli_date
        }
        k_count <- 0
        while (next_ring_date < survey_date) {
            age_at_ring_date <- next_ring_date - birth_date
            if (age_at_ring_date >= annuli_min_age) {
                k_count <- k_count + 1
            }
            next_ring_date <- next_ring_date + 1
        }
        return(k_count)
    })
}

#' Generate model predictions for a specific survey date
#' Encapsulates the prediction pipeline to obtain P(K | length) for one survey.
#' @param survey_date Numeric survey date (e.g., 2023.25).
#' @param G Impulse-response matrix from the single-cohort simulation
#'   (rows = ages, cols = length classes).
#' @param a Numeric vector of high-resolution ages corresponding to rows of `G`.
#' @param l Numeric vector of length-class centers corresponding to columns of `G`.
#' @param mu Mean spawning date as a fraction of a year in \[0, 1).
#' @param kappa Spawning concentration parameter.
#' @param annuli_date Ring formation day as fraction of a year in \[0, 1).
#' @param annuli_min_age Minimum age (years) at which the first ring can form.
#' @return A matrix of proportions with rows named by `l` and columns by K bins.
#' @export
#' @examples
#' # Minimal schematic example (using toy inputs)
#' a <- seq(0, 3, length.out = 5)
#' l <- seq(10, 30, length.out = 3)
#' G <- matrix(abs(sin(outer(a, l, "+"))), nrow = length(a))
#' generate_model_predictions_for_date(2023.5, G, a, l, mu = 0.5, kappa = 3,
#'                                     annuli_date = 0.25, annuli_min_age = 0.5)
generate_model_predictions_for_date <- function(
        survey_date, G, a, l, mu, kappa, annuli_date, annuli_min_age) {
    # Population Convolution
    birth_dates <- survey_date - a
    spawning_weights <- spawning_density(birth_dates, mu, kappa)
    N_pop <- diag(spawning_weights) %*% G

    # Observation Convolution
    k_for_each_age <- calculate_K(a, survey_date, annuli_date, annuli_min_age)
    max_K <- max(k_for_each_age)
    k_bins <- 0:max_K
    N_model <- matrix(0, nrow = length(l), ncol = length(k_bins))
    rownames(N_model) <- l
    colnames(N_model) <- k_bins

    for (k_val in k_bins) {
        age_indices <- which(k_for_each_age == k_val)
        if (length(age_indices) > 0) {
            pop_subset <- N_pop[age_indices, , drop = FALSE]
            N_model[, k_val + 1] <- colSums(pop_subset)
        }
    }

    # Step D: Convert to proportions
    P_model_K_given_l <- prop.table(N_model, margin = 1)
    P_model_K_given_l[is.nan(P_model_K_given_l)] <- 0
    return(P_model_K_given_l)
}

#' Calculate and Aggregate Signed Log-Likelihood Contributions
#' Loops through surveys, calculates signed NLL for each, and aggregates the results.
#' @param surveys A list of data frames, split by survey date.
#' @inheritParams generate_model_predictions_for_date
#' @return A data frame with the total NLL contribution for each Length-K bin.
calculate_and_aggregate_likelihood <- function(surveys, G, a, l, mu, kappa,
                                               annuli_date, annuli_min_age) {

    log_lik_contributions <- list()

    for (survey_date_str in names(surveys)) {
        survey_date_current <- as.numeric(survey_date_str)

        current_obs_df <- surveys[[survey_date_str]]

        # 1. Generate model predictions for this date
        P_model <- generate_model_predictions_for_date(
            survey_date_current, G, a, l, mu, kappa, annuli_date, annuli_min_age
        )

        # 2. Get sample sizes per length for this survey

        sample_sizes <- current_obs_df %>%
            group_by(Length) %>%
            summarise(N = sum(count, na.rm = TRUE), .groups = 'drop') %>%
            mutate(Length = as.factor(Length))

        # 3. Convert model probabilities to a long data frame for joining
        P_model_df <- as.data.frame.table(P_model)
        colnames(P_model_df) <- c("Length", "K", "Prob")

        # 4. Join all data together
        likelihood_df <- current_obs_df %>%
            left_join(sample_sizes, by = "Length") %>%
            left_join(P_model_df, by = c("Length", "K"))

        # 5. Calculate the signed negative log-likelihood contribution
        epsilon <- 1e-9 # To prevent log(0)
        likelihood_df <- likelihood_df %>%
            mutate(
                Expected = N * Prob,
                NegLogLik = - (count * log(Prob + epsilon)),
                SignedNegLogLik = sign(count - Expected) * NegLogLik
            )

        log_lik_contributions[[survey_date_str]] <- likelihood_df
    }

    # Aggregate contributions across all surveys
    all_contributions_df <- do.call(rbind, log_lik_contributions)

    total_contributions <- all_contributions_df %>%
        group_by(Length, K) %>%
        summarise(
            TotalObserved = sum(count, na.rm = TRUE),
            TotalExpected = sum(Expected, na.rm = TRUE),
            TotalNegLogLik = sum(NegLogLik, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        mutate(
            # The final signed NLL is the total misfit, with the sign determined
            # by the overall difference between observed and expected counts.
            SignedNegLogLik = sign(TotalObserved - TotalExpected) * TotalNegLogLik
        )

    return(total_contributions)
}

#' Simulate a sample from model predictions
#' Draws K values using multinomial sampling with probabilities P(K | length)
#' for each length observed in a given survey.
#' @param P_model_K_given_l Matrix of predicted proportions P(K | length).
#' @param survey_obs A data frame for one survey with at least a `Length` column.
#' @return An integer vector of simulated K values aligned with `survey_obs` rows.
#' @keywords internal
#' @examples
#' set.seed(1)
#' P <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' rownames(P) <- c("20", "30"); colnames(P) <- c("0", "1")
#' survey_obs <- data.frame(Length = c(20, 20, 30, 30, 30))
#' simulate_sample_from_model(P, survey_obs)
simulate_sample_from_model <- function(P_model_K_given_l, survey_obs) {
    simulated_K <- integer(nrow(survey_obs))

    # Get the unique length classes that were actually sampled in this survey
    unique_lengths <- unique(survey_obs$Length)

    # Get the possible K values from the column names of the proportion matrix.
    k_values <- as.numeric(colnames(P_model_K_given_l))

    # For each unique length class...
    for (len_val in unique_lengths) {
        indices_to_fill <- which(survey_obs$Length == len_val)
        n_fish <- length(indices_to_fill)
        len_char <- as.character(len_val)

        # Get the model's predicted proportions of K for this length
        k_proportions <- P_model_K_given_l[len_char, , drop = TRUE]

        # Check if the length class exists in the model predictions and has valid probabilities
        if (len_char %in% rownames(P_model_K_given_l) && sum(k_proportions, na.rm = TRUE) > 0) {

            # Perform a multinomial random sample
            sim_counts <- rmultinom(1, size = n_fish, prob = k_proportions)

            # Create a vector of the simulated K values
            simulated_k_vector <- rep(k_values, sim_counts)

            # Assign these simulated K's to the correct rows in the output.
            # The length of simulated_k_vector is guaranteed to match n_fish.
            if (length(simulated_k_vector) > 0) {
                simulated_K[indices_to_fill] <- simulated_k_vector
            }
        }
    }
    return(simulated_K)
}



#' Pre-process length-at-age data frame
#'
#' Filters and aggregates age-at-length data for a species, rounding lengths down
#' to the nearest cm and aggregating counts by quarter. Returns a tidy data frame
#' with columns `survey_date`, `Length`, `K`, and `count`.
#' @param params A `mizer::MizerParams` object.
#' @param species Species name as in `species_params(params)$species`.
#' @param age_at_length Data frame with at least the columns `Scientific_name`,
#'   `Quarter`, `LngtClass`, `Age`, and `CANoAtLngt`.
#' @return A data frame with the columns described above, one row per
#'   quarter-length-age combination.
#' @export
#' @examples
#' # Using package data would typically look like:
#' # df <- preprocess_length_at_age(params, species = "Cod", age_at_length = your_df)
preprocess_length_at_age <- function(params, species, age_at_length) {
    sci_name <- species_params(params)[species, "SciName"]
    survey_dates <- c(0.125, 0.375, 0.625, 0.875)
    age_at_length |>
        filter(.data$Scientific_name == sci_name) |>
        # remove rows with NA in any column
        filter(!is.na(.data$LngtClass) & !is.na(.data$Age) & !is.na(.data$CANoAtLngt) &
                   !is.na(.data$Quarter)) |>
        # round down to cm
        mutate(LngtClass = floor(.data$LngtClass)) |>
        # aggregate counts by quarter
        group_by(.data$Quarter, .data$LngtClass, .data$Age) |>
        summarise(CANoAtLngt = sum(.data$CANoAtLngt, na.rm = TRUE), .groups = "drop") |>
        transmute(survey_date = survey_dates[.data$Quarter],
                  Length = as.integer(.data$LngtClass),
                  K = as.integer(.data$Age),
                  count = .data$CANoAtLngt)
}

#' Build a length rebinning matrix
#' Computes the fraction of each model length bin that overlaps each survey
#' length bin, producing a matrix suitable for rebinning/aggregation.
#' @param l_model Numeric vector of model bin edges (strictly increasing).
#' @param l_survey Numeric vector of survey bin edges (strictly increasing).
#' @return A matrix of dimension `(length(l_model) - 1) x (length(l_survey) - 1)`
#'   where each entry is the fraction of the model bin width overlapping a
#'   survey bin.
#' @keywords internal
length_rebinning_matrix <- function(l_model, l_survey) {
    ## LENGTH aggregation B : (surveyLen × modelLen) ----
    # model length bins are defined by l_model
    low_L  <- l_model[-length(l_model)]
    high_L <- l_model[-1]
    # survey length bins are defined by l_survey
    low_S  <- l_survey[-length(l_survey)]
    high_S <- l_survey[-1]

    B <- matrix(0, nrow = length(high_L), ncol = length(high_S))
    for (l in seq_along(high_L)) {
        for (j in seq_along(high_S)) {
            overlap <- max(0, min(high_L[l], high_S[j]) - max(low_L[l], low_S[j]))
            B[l, j] <- overlap / (high_L[l] - low_L[l])   # fraction of model bin
        }
    }
    return(B)
}
