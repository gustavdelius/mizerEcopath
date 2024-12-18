#' Plot the production against species
#'
#' @param params A MizerParams object
#' @export
plotProductionVsSpecies <- function(params) {
    no_sp <- length(params@species_params$species)
    observed <- params@species_params$production_observed
    if (is.null(observed)) observed <- rep(NA, no_sp)

    model <- getProduction(params)

    species <- factor(params@species_params$species)
    df <- rbind(
        data.frame(Species = species,
                   Type = "Observation",
                   Production = observed,
                   other = model),
        data.frame(Species = species,
                   Type = "Model",
                   Production = model,
                   other = observed)
    )
    # Get rid of unobserved entries
    df <- df[df$Production > 0 & !is.na(df$Production), ]

    ggplot(df, aes(x = Species, y = Production)) +
        geom_point(aes(shape = Type), size = 4) +
        geom_linerange(aes(ymin = Production, ymax = other, colour = Species)) +
        scale_y_continuous(name = "Production [g/year]", trans = "log10",
                           breaks = log_breaks()) +
        scale_colour_manual(values = getColours(params)) +
        scale_shape_manual(values = c(Model = 1, Observation = 15)) +
        guides(colour = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
