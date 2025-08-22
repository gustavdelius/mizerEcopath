#' Plot Negative Log-Likelihood Contributions
#' Creates a heatmap of the contribution of each cell to the NLL.
#' @param contributions_df A data frame with Length, K, and TotalNegLogLik.
#' @return A ggplot object.
plot_log_likelihood <- function(contributions_df) {

    # Convert factors to numeric for plotting
    contributions_df$Length <- as.numeric(as.character(contributions_df$Length))
    contributions_df$K <- as.numeric(as.character(contributions_df$K))

    # Create the plot
    p <- ggplot(contributions_df, aes(x = Length, y = K, fill = TotalNegLogLik)) +
        geom_tile() +
        # Using a sequential color scale is more appropriate for NLL (always >= 0)
        scale_fill_viridis_c(
            name = "Neg Log-Lik\nContribution",
            option = "plasma", # A visually appealing color scale
            direction = -1 # Puts darker colors for higher values
        ) +
        labs(
            title = "Model Fit Diagnostic: Negative Log-Likelihood",
            subtitle = "Darker colors indicate cells with poorer model fit (higher 'surprise')",
            x = "Fish Length (cm)",
            y = "Otolith Ring Count (K)"
        ) +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9)
        )

    return(p)
}

#' Calculate and plot Pearson residuals
#' Creates a heatmap of Pearson residuals to visualize model fit in K-by-length
#' space, with overlaid average K lines for observed and simulated data.
#' @param observed_df Data frame of observed data with columns `Length` and `K`.
#' @param simulated_df Data frame of simulated data with columns `Length` and `K`.
#' @return A `ggplot2` object.
#' @keywords internal
#' @examples
#' # Toy example
#' obs <- data.frame(Length = c(20, 20, 25, 30), K = c(0, 1, 1, 2))
#' sim <- data.frame(Length = c(20, 25, 25, 30, 30), K = c(0, 0, 1, 1, 2))
#' calculate_and_plot_residuals(obs, sim)
calculate_and_plot_residuals <- function(observed_df, simulated_df) {

    # --- 1. Create contingency tables (counts of fish by Length and K) ---
    obs_counts <- as.data.frame(table(observed_df$Length, observed_df$K))
    sim_counts <- as.data.frame(table(simulated_df$Length, simulated_df$K))

    colnames(obs_counts) <- c("Length", "K", "Observed")
    colnames(sim_counts) <- c("Length", "K", "Simulated")

    # --- 2. Merge the two tables ---
    # A full join ensures that we keep all cells, even if one of the
    # tables has a zero count for that cell.
    all_counts <- full_join(obs_counts, sim_counts, by = c("Length", "K"))

    # Replace NA values with 0 for counts
    all_counts$Observed[is.na(all_counts$Observed)] <- 0
    all_counts$Simulated[is.na(all_counts$Simulated)] <- 0

    # --- 3. Calculate Pearson residuals ---
    # Formula: (Observed - Simulated) / sqrt(Simulated)
    # We add a small epsilon to the denominator to avoid division by zero.
    epsilon <- 1e-9
    all_counts$Residual <- (all_counts$Observed - all_counts$Simulated) /
        sqrt(all_counts$Simulated + epsilon)

    # --- 4. Prepare for plotting ---
    # Convert factor columns to numeric for plotting
    all_counts$Length <- as.numeric(as.character(all_counts$Length))
    all_counts$K <- as.numeric(as.character(all_counts$K))

    # Cap residuals for better color scaling, as extreme outliers can
    # wash out the interesting patterns. A cap of +/- 4 is common.
    res_limit <- 4
    all_counts$Residual_capped <- pmax(-res_limit, pmin(res_limit, all_counts$Residual))

    # --- 5. Calculate average K for each length bin ---
    # Calculate average K for observed data
    obs_avg_k <- observed_df %>%
        group_by(.data$Length) %>%
        summarise(
            avg_k = weighted.mean(.data$K, w = rep(1, dplyr::n()), na.rm = TRUE),
            .groups = "drop"
        )

    # Calculate average K for simulated data
    sim_avg_k <- simulated_df %>%
        group_by(.data$Length) %>%
        summarise(
            avg_k = weighted.mean(.data$K, w = rep(1, dplyr::n()), na.rm = TRUE),
            .groups = "drop"
        )

    # --- 6. Create the plot using ggplot2 ---
    p <- ggplot(all_counts, aes(x = .data$Length, y = .data$K, fill = .data$Residual_capped)) +
        geom_tile() + # This creates the heatmap tiles
        # Add lines for average K
        geom_line(data = obs_avg_k, aes(x = .data$Length, y = .data$avg_k, color = "Observed Avg K"),
                  linewidth = 2, inherit.aes = FALSE) +
        geom_line(data = sim_avg_k, aes(x = .data$Length, y = .data$avg_k, color = "Simulated Avg K"),
                  linewidth = 1, inherit.aes = FALSE) +
        scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = 0,
            limit = c(-res_limit, res_limit),
            name = "Pearson\nResidual"
        ) +
        scale_color_manual(
            name = "Average K",
            values = c("Observed Avg K" = "#000000", "Simulated Avg K" = "#ff7f0e")
        ) +
        labs(
            title = "Pearson Residuals of Model Fit",
            subtitle = "Blue = Model Underestimates, Red = Model Overestimates",
            x = "Fish Length (cm)",
            y = "Otolith Ring Count (K)"
        ) +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right"
        )

    return(p)
}

#' Plot observed vs simulated ring counts by length
#' Generates a residual heatmap and overlays average K-by-length for both
#' observed and simulated data, for a single species.
#' @param params A `mizer::MizerParams` object.
#' @param species Species name as in `species_params(params)$species`.
#' @param age_at_length Data frame of raw age-at-length observations; will be
#'   preprocessed internally by `preprocess_length_at_age()`.
#' @return A `ggplot2` object suitable for display in Shiny or saving.
#' @export
#' @examples
#' # In practice provide a real `age_at_length` table for the species
#' # p <- plotAge(params, species = "Cod", age_at_length = df)
plotSimulatedAge <- function(params, species, age_at_length) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # Preprocess observed data
    observed_df <- preprocess_length_at_age(params, species, age_at_length)

    # Simulate age data
    simulated_df <- simulateAge(params, species, observed_df)

    # Calculate and plot residuals
    p <- calculate_and_plot_residuals(observed_df, simulated_df)

    return(p)
}

#' Plot the likelihood of the observed age data
#'
#' @param params A `mizer::MizerParams` object.
#' @param species Species name as in `species_params(params)$species`.
#' @param age_at_length Data frame of raw age-at-length observations; will be
#'   preprocessed internally by `preprocess_length_at_age()`.
#' @return A `ggplot2` object suitable for display in Shiny or saving.
#' @export
#' @examples
#' # In practice provide a real `age_at_length` table for the species
#' # p <- plotAge(params, species = "Cod", age_at_length = df)
plotAgeLikelihood <- function(params, species, age_at_length) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # Preprocess observed data
    observed_df <- preprocess_length_at_age(params, species, age_at_length)

    # Simulate age data
    logLik <- getLogLik(params, species, observed_df)

    # Calculate and plot residuals
    plot_log_likelihood(logLik)
}
