#' von Mises Probability Density Function
#' A custom implementation to avoid external package dependencies.
#' @param x Angle in radians.
#' @param mu Mean direction in radians.
#' @param kappa Concentration parameter.
#' @return The probability density.
von_mises_pdf <- function(x, mu, kappa) {
    denominator <- 2 * pi * besselI(kappa, 0)
    numerator <- exp(kappa * cos(x - mu))
    return(numerator / denominator)
}

#' Spawning Density Function S(d)
#' Calculates the relative spawning intensity for any given numeric date.
#' @param numeric_dates A vector of numeric dates (e.g., 2023.45).
#' @param mu The mean spawning date.
#' @param kappa The concentration parameter for the spawning distribution.
#' @return A numeric vector of relative spawning intensities.
spawning_density <- function(numeric_dates, mu, kappa) {
    day_fraction <- numeric_dates %% 1
    day_rad <- day_fraction * 2 * pi
    mu_rad <- mu * 2 * pi
    density <- von_mises_pdf(day_rad, mu = mu_rad, kappa = kappa)
    return(density)
}


#' Age-to-Ring Mapping Function Calculate_K(a)
#' For a given true age, calculates the deterministic number of otolith rings.
#' @param age_in_years A numeric vector of true ages.
#' @param survey_date The numeric representation of the survey date.
#' @param t_r The numeric representation of the ring formation day.
#' @param a_min The minimum age for a fish to form its first ring.
#' @return An integer vector of the same length, with the calculated K for each age.
calculate_K <- function(age_in_years, survey_date, t_r, a_min) {
    sapply(age_in_years, function(age) {
        if (age < 0) return(0)
        birth_date <- survey_date - age
        birth_year <- floor(birth_date)
        next_ring_date <- birth_year + t_r
        if (next_ring_date <= birth_date) {
            next_ring_date <- (birth_year + 1) + t_r
        }
        k_count <- 0
        while (next_ring_date < survey_date) {
            age_at_ring_date <- next_ring_date - birth_date
            if (age_at_ring_date >= a_min) {
                k_count <- k_count + 1
            }
            next_ring_date <- next_ring_date + 1
        }
        return(k_count)
    })
}

#' Generate Model Predictions for a Specific Survey Date
#' This function encapsulates the entire prediction pipeline for one survey.
#' @param survey_date The numeric survey date.
#' @param G The impulse response matrix from the single cohort simulation.
#' @param a The vector of high-resolution ages.
#' @param l The vector of length classes.
#' @param mu Mean spawning date.
#' @param kappa Spawning concentration parameter.
#' @param t_r
#' @param a_min Minimum age at which first ring can form
#' @return A matrix of proportions, P(K|s), for the given survey date.
#' @export
generate_model_predictions_for_date <- function(
        survey_date, G, a, l, mu, kappa, t_r, a_min) {
    # Population Convolution
    birth_dates <- survey_date - a
    spawning_weights <- spawning_density(birth_dates, mu, kappa)
    N_pop <- diag(spawning_weights) %*% G

    # Observation Convolution
    k_for_each_age <- calculate_K(a, survey_date, t_r, a_min)
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

#' Simulate a Sample from Model Predictions
#' Uses the model proportions and observed sample sizes to generate a simulated dataset.
#' @param P_model_K_given_l The predicted proportions from the model.
#' @param survey_obs A data frame of observations for a single survey.
#' @return A vector of simulated K values for the given lengths.
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



#'Pre-process length-at-age data frame
preprocess_length_at_age <- function(params, species, length_at_age) {
    sci_name <- species_params(params)[species, "SciName"]
    survey_dates <- c(0.125, 0.375, 0.625, 0.875)
    age_at_length |>
        filter(Scientific_name == sci_name) |>
        # remove rows with NA in any column
        filter(!is.na(LngtClass) & !is.na(Age) & !is.na(CANoAtLngt) &
                   !is.na(Quarter)) |>
        # round down to cm
        mutate(LngtClass = floor(LngtClass)) |>
        # aggregate counts by quarter
        group_by(Quarter, LngtClass, Age) |>
        summarise(CANoAtLngt = sum(CANoAtLngt, na.rm = TRUE), .groups = "drop") |>
        transmute(survey_date = survey_dates[Quarter],

                  Length = as.integer(LngtClass),
                  K = as.integer(Age),
                  count = CANoAtLngt)
}

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
