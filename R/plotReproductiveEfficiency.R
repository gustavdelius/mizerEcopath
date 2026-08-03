#' Plot the efficiency of reproduction
#'
#' Shows, for each species, which fraction of the energy that adults invest
#' into reproduction is actually turned into egg biomass in the model. A
#' value larger than 1 (above the red line) is energetically impossible and
#' signals that the model is asking reproduction to create biomass out of
#' nothing.
#'
#' @section The quantity that is plotted:
#'
#' Mature individuals divert a part of their income into reproduction. Summing
#' this over the population gives the total rate \eqn{E_R} (grams per year) at
#' which the species invests into reproduction. Mizer assumes that half of that
#' investment is made by females, so \eqn{E_R / 2} is the energy that is
#' available for making eggs.
#'
#' Not all of it arrives in the size spectrum as eggs. Mizer removes energy in
#' two steps:
#'
#' 1. The reproductive efficiency `erepro` \eqn{= \epsilon} accounts for all
#'    losses on the way from the gonads to a viable egg (metabolic cost of
#'    producing gametes, unfertilised or unviable eggs, ...). It gives the
#'    density-independent rate of egg production
#'    \deqn{R_{di} = \frac{\epsilon\, E_R}{2\, w_{min}},}
#'    see [mizer::getRDI()].
#' 2. Density dependence, usually a Beverton-Holt stock-recruitment curve,
#'    reduces this further to the realised rate \eqn{R_{dd}}, see
#'    [mizer::getRDD()]. The surviving fraction is \eqn{R_{dd} / R_{di}}, which
#'    is the same as `1 - reproduction_level`, see
#'    [mizer::setBevertonHolt()].
#'
#' The eggs enter the spectrum at the egg weight \eqn{w_{min}}, so the model
#' produces egg biomass at the rate \eqn{R_{dd}\, w_{min}}. The
#' **reproductive efficiency** plotted by this function is that biomass rate
#' expressed as a fraction of the energy that was available for it:
#' \deqn{\frac{R_{dd}\, w_{min}}{E_R / 2} = \epsilon \, \frac{R_{dd}}{R_{di}}.}
#'
#' @section How to read the plot:
#'
#' For each species the plot shows two points joined by a vertical line:
#'
#' * **Potential**: the efficiency `erepro` alone, i.e. what the species would
#'   achieve if density dependence removed no eggs at all.
#' * **Realised**: the efficiency after density dependence, i.e. the quantity
#'   defined above.
#'
#' The length of the line is therefore the amount of reproduction that density
#' dependence removes, and the position of the lower point tells you how close
#' the species is to the energetic limit.
#'
#' The red line at 1 is that limit: it marks the case where every gram of
#' female reproductive investment ends up as egg biomass. Species sitting
#' close below it have no room left --- any further loss of reproductive
#' output would have to be compensated by an impossible efficiency. Species
#' sitting far below it are reproducing well within their energy budget,
#' either because `erepro` is small or because density dependence is strong.
#' Species above it indicate a mis-specified model: the observed abundance of
#' small individuals cannot be sustained by the reproductive investment that
#' the model's adults are able to make, so either their food intake, their
#' investment into reproduction, or the mortality of the small individuals
#' needs to be revisited.
#'
#' @param params A MizerParams object
#' @param species The species to be included. Optional. By default all
#'   foreground species are included. A vector of species names, or a numeric
#'   vector of species indices, or a logical vector indicating for each species
#'   whether it is to be included (length must equal the number of species in
#'   the model).
#' @param return_data A boolean value that determines whether the data frame
#'   underlying the plot is returned instead of the plot itself. Defaults to
#'   FALSE. The data frame has one row per species and the columns `Species`,
#'   `erepro` (the potential efficiency), `dd_survival` (the fraction
#'   \eqn{R_{dd}/R_{di}} of eggs surviving density dependence) and `efficiency`
#'   (the realised efficiency, the product of the previous two).
#'
#' @return A ggplot2 object, or a data frame if `return_data = TRUE`.
#' @export
#' @seealso [getReproductiveEfficiency()] returns the realised efficiency
#'   without plotting it, also for a MizerSim object.
#' @examples
#' plotReproductiveEfficiency(celtic_params)
plotReproductiveEfficiency <- function(params, species = NULL,
                                       return_data = FALSE) {
    params <- validParams(params)
    sel <- valid_species_arg(params, species, return.logical = TRUE,
                             error_on_empty = TRUE)
    sp <- params@species_params[sel, ]

    efficiency <- unname(getReproductiveEfficiency(params)[sel])

    df <- data.frame(Species = factor(sp$species, levels = sp$species),
                     erepro = sp$erepro,
                     # Fraction of eggs surviving density dependence,
                     # = 1 - reproduction_level
                     dd_survival = efficiency / sp$erepro,
                     efficiency = efficiency)
    rownames(df) <- NULL
    if (return_data) return(df)

    dfl <- rbind(
        data.frame(Species = df$Species, Type = "Potential",
                   value = df$erepro),
        data.frame(Species = df$Species, Type = "Realised",
                   value = df$efficiency)
    )
    dfl$Type <- factor(dfl$Type, levels = c("Potential", "Realised"))
    dfl <- dfl[!is.na(dfl$value) & dfl$value > 0, ]

    ggplot(dfl, aes(x = Species, y = value)) +
        geom_hline(yintercept = 1, colour = "red") +
        geom_segment(data = df,
                     aes(x = Species, xend = Species,
                         y = erepro, yend = efficiency, colour = Species),
                     inherit.aes = FALSE, linewidth = 1) +
        geom_point(aes(shape = Type), size = 3) +
        scale_y_continuous(name = "Reproductive efficiency",
                           trans = "log10", breaks = log_breaks()) +
        scale_shape_manual(values = c(Potential = 1, Realised = 16),
                           name = NULL) +
        scale_colour_manual(values = getColours(params)) +
        guides(colour = "none") +
        theme(text = element_text(size = 12),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
