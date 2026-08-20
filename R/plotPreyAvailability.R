#' Plot the prey available to a predator of a given size
#'
#' Shows, for a predator of a given species and size, how the prey that it
#' encounters is distributed over prey sizes. Four curves are shown, each
#' normalised to integrate to 1 over log prey weight:
#'
#' * **Feeding kernel**: the predation kernel \eqn{\phi(w_p/w)} alone, i.e.
#'   the size preference of the predator, irrespective of what is there.
#' * **Number density**: the density of available prey individuals in prey
#'   weight, i.e. the kernel times the abundance of the prey community.
#' * **Number density in log w**: the same, but as a density in log weight,
#'   so that the area under the curve over a size range is proportional to
#'   the number of available prey in that range.
#' * **Biomass density in log w**: the corresponding biomass density, so that
#'   the area under the curve is proportional to the available prey biomass.
#'
#' The prey community is the sum of the resource spectrum and the fish
#' spectra, each weighted by the predator's interaction strength with it,
#' see [mizer::setInteraction()]. The size of the predator is marked by a
#' point on the x-axis.
#'
#' @param params A MizerParams object
#' @param species The name of the predator species. Optional if the model has
#'   only a single species.
#' @param pred_size The weight of the predator in grams. Defaults to the
#'   maturity weight `w_mat` of the species.
#' @param return_data A boolean value that determines whether the data frame
#'   underlying the plot is returned instead of the plot itself. Defaults to
#'   FALSE. The data frame has the variables `w`, `value` and `Type`.
#'
#' @return A ggplot2 object, or a data frame if `return_data = TRUE`. The
#'   `plotly` version returns a plotly object.
#' @param ... Additional arguments passed to `plotPreyAvailability()` by the
#'   plotly wrapper.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions], [mizer::plotDiet()]
#' @examples
#' plotPreyAvailability(celtic_params, species = "Cod", pred_size = 1000)
plotPreyAvailability <- function(params, species = NULL, pred_size = NULL,
                                 return_data = FALSE) {
    params <- validParams(params)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("You can only plot the prey availability for ",
             "one species at a time.")
    }
    sp_idx <- which(params@species_params$species == species)
    if (is.null(pred_size)) {
        pred_size <- params@species_params$w_mat[[sp_idx]]
    }
    assert_that(is.number(pred_size))
    # Index of the size bin holding the predator
    w_idx <- max(1, sum(params@w <= pred_size))

    dx <- log(params@w_full[2]) - log(params@w_full[1])
    no_w_full <- length(params@w_full)
    fish_idx <- (no_w_full - length(params@w) + 1):no_w_full
    # Abundance of all prey, weighted by how strongly the predator interacts
    # with it
    total_n <- params@initial_n_pp *
        params@species_params$interaction_resource[[sp_idx]]
    total_n[fish_idx] <- total_n[fish_idx] +
        params@interaction[sp_idx, ] %*% params@initial_n

    kernel <- pred_kernel(params)[sp_idx, w_idx, ]
    n <- total_n * kernel
    n_logw <- n * params@w_full
    b_logw <- n_logw * params@w_full
    # Turn each curve into a probability density in log w
    normalise <- function(value) value / sum(value * dx)

    plot_dat <- rbind(
        data.frame(w = params@w_full, value = normalise(kernel),
                   Type = "Feeding kernel"),
        data.frame(w = params@w_full, value = normalise(n),
                   Type = "Number density"),
        data.frame(w = params@w_full, value = normalise(n_logw),
                   Type = "Number density in log w"),
        data.frame(w = params@w_full, value = normalise(b_logw),
                   Type = "Biomass density in log w")
    )
    if (return_data) return(plot_dat)

    params <- setPreyAvailabilityColours(params)
    plotDataFrame(plot_dat, params,
                  xlab = "Prey weight [g]", ylab = "Density",
                  xtrans = "log10") +
        annotate("point", x = pred_size, y = 0, size = 4, colour = "blue")
}

#' @rdname plotPreyAvailability
#' @export
plotlyPreyAvailability <- function(params, species = NULL, pred_size = NULL,
                                   ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotPreyAvailability", argg),
             tooltip = c("w", "value", "Type"))
}

#' Set the colours used for the curves in `plotPreyAvailability()`
#'
#' Only sets the colour of those curves for which no colour has been chosen
#' already, so that the user can override them with [mizer::setColours()].
#' @param params A MizerParams object
#' @return The MizerParams object with the colours set
#' @noRd
setPreyAvailabilityColours <- function(params) {
    colours <- c("Feeding kernel" = "black",
                 "Number density" = "#E69F00",
                 "Number density in log w" = "#56B4E9",
                 "Biomass density in log w" = "#009E73")
    missing <- setdiff(names(colours), names(getColours(params)))
    if (length(missing) > 0) {
        params <- setColours(params, colours[missing])
    }
    missing <- setdiff(names(colours), names(getLinetypes(params)))
    if (length(missing) > 0) {
        params <- setLinetypes(params, setNames(rep("solid", length(missing)),
                                                missing))
    }
    params
}
