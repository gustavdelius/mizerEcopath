#' Plot Signed Negative Log-Likelihood Contributions with Mean K Lines
#' Creates a heatmap of the signed contribution of each cell to the NLL,
#' and overlays the mean observed and expected K values.
#' @param contributions_df A data frame from calculate_and_aggregate_likelihood.
#' @return A ggplot object.
plot_log_likelihood <- function(contributions_df) {

    # Convert factors to numeric for plotting
    contributions_df$Length <- as.numeric(as.character(contributions_df$Length))
    contributions_df$K <- as.numeric(as.character(contributions_df$K))

    # Calculate mean K lines ---
    mean_lines_df <- contributions_df %>%
        group_by(Length) %>%
        summarise(
            # Weighted mean of K by the number of observed/expected fish
            Observed = sum(K * TotalObserved) / sum(TotalObserved),
            Model = sum(K * TotalExpected) / sum(TotalExpected),
            .groups = 'drop'
        ) %>%
        # Reshape the data for easy plotting with ggplot
        pivot_longer(
            cols = c("Observed", "Model"),
            names_to = "Source",
            values_to = "MeanK"
        )

    # Determine symmetric color scale limits for the heatmap
    max_abs_val <- max(abs(contributions_df$SignedNegLogLik), na.rm = TRUE)

    # Create the plot
    p <- ggplot(contributions_df, aes(x = Length, y = K)) +
        # Heatmap layer for the misfit
        geom_tile(aes(fill = SignedNegLogLik)) +
        # Line layers for the mean K values
        geom_line(
            data = mean_lines_df,
            aes(y = MeanK, color = Source, linetype = Source),
            size = 1
        ) +
        # --- Define scales and labels ---
        scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = 0,
            limit = c(-max_abs_val, max_abs_val),
            name = "Direction & Magnitude\nof Misfit (NLL)"
        ) +
        scale_color_manual(
            name = "Mean K",
            values = c("Observed" = "black", "Model" = "darkgreen")
        ) +
        scale_linetype_manual(
            name = "Mean K",
            values = c("Observed" = "solid", "Model" = "solid")
        ) +
        labs(
            title = "Model Fit Diagnostic",
            subtitle = "Color shows direction (Blue: Obs > Exp, Red: Obs < Exp). Intensity shows magnitude of misfit.",
            x = "Fish Length (cm)",
            y = "Annuli Count (K)"
        ) +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 18),
            legend.position = "right"
        )

    return(p)
}

#' Plot the likelihood of the observed age at length data
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
    if (is.null(age_at_length)) {
        return(NULL)
    }
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # Preprocess observed data
    observed_df <- preprocess_length_at_age(params, species, age_at_length)

    # Simulate age data
    logLik <- getLogLik(params, species, observed_df)

    # Calculate and plot residuals
    plot_log_likelihood(logLik)
}
