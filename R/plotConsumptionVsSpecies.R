#' Plot the Consumption against species
#'
#' @param params A MizerParams object
#' @export
plotConsumptionVsSpecies <- function(params) {
    no_sp <- length(params@species_params$species)
    observed <- params@species_params$consumption_observed
    if (is.null(observed)) observed <- rep(NA, no_sp)

    model <- getConsumption(params)

    species <- factor(params@species_params$species)
    df <- rbind(
        data.frame(Species = species,
                   Type = "Observation",
                   Consumption = observed,
                   other = model),
        data.frame(Species = species,
                   Type = "Model",
                   Consumption = model,
                   other = observed)
    )
    # Get rid of unobserved entries
    df <- df[df$Consumption > 0 & !is.na(df$Consumption), ]

    df_lr <- df[!is.na(df$other) & df$other > 0, ]
    p <- ggplot(df, aes(x = Species, y = Consumption)) +
        geom_point(aes(shape = Type), size = 4)
    if (nrow(df_lr) > 0) {
        p <- p + geom_linerange(data = df_lr,
                           aes(ymin = Consumption, ymax = other, colour = Species))
    }
    p +
        scale_y_continuous(name = "Consumption [g/year]", trans = "log10",
                           breaks = log_breaks()) +
        scale_colour_manual(values = getColours(params)) +
        scale_shape_manual(values = c(Model = 1, Observation = 15)) +
        guides(colour = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
