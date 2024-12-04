# Plot observed catches against modelled catches
#' @export
plot_catch <- function(params, species = 1, catch) {
    params <- validParams(params)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Only one species can be updated at a time.")
    }
    sp <- species_params(params)
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]

    # Validate catch data frame
    catch <- valid_catch(catch, species)

    heights <- catch$count / catch$dl / sum(catch$count)

    # Create an empty plot with the correct x and y limits
    plot(NULL, xlim = c(min(catch$length), max(catch$length + catch$dl)),
         ylim = c(0, max(heights)), xlab = "Length [cm]", ylab = "Catch Density",
         main = species)

    # Draw the bars using rect()
    rect(catch$length, 0, catch$length + catch$dl, heights,
         col = "blue", border = "blue")

    # We use lengths instead of weights in the catch data
    lengths <- (params@w / sps$a)^(1/sps$b)
    # Catch density as a function of weight
    model_catch <- params@initial_n[sp_select, ] * getFMort(params)[sp_select, ]
    # Normalise the density
    model_catch <- model_catch  / sum(model_catch * params@dw)
    # Convert to density as a function of length
    model_catch <- model_catch * sps$b * params@w / lengths

    # Add the fitted PDF to the plot
    lines(lengths, model_catch, col = 'red', lwd = 2)
    legend('topright', legend = c('Observed', 'Modelled'),
           col = c('blue', 'red'), lwd = 2)
}
