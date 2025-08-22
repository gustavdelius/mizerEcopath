#' Calculate and Plot Pearson Residuals
#' Creates a heatmap of Pearson residuals to visualize model fit.
#' @param observed_df The data frame of observed data (Length, K).
#' @param simulated_df The data frame of simulated data (Length, K).
#' @return A ggplot object.
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
        group_by(Length) %>%
        summarise(
            avg_k = weighted.mean(K, w = rep(1, n()), na.rm = TRUE),
            .groups = "drop"
        )

    # Calculate average K for simulated data
    sim_avg_k <- simulated_df %>%
        group_by(Length) %>%
        summarise(
            avg_k = weighted.mean(K, w = rep(1, n()), na.rm = TRUE),
            .groups = "drop"
        )

    # --- 6. Create the plot using ggplot2 ---
    p <- ggplot(all_counts, aes(x = Length, y = K, fill = Residual_capped)) +
        geom_tile() + # This creates the heatmap tiles
        # Add lines for average K
        geom_line(data = obs_avg_k, aes(x = Length, y = avg_k, color = "Observed Avg K"),
                  linewidth = 2, inherit.aes = FALSE) +
        geom_line(data = sim_avg_k, aes(x = Length, y = avg_k, color = "Simulated Avg K"),
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

plotAge <- function(params, species, age_at_length) {
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
