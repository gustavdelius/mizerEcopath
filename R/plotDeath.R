# This is a copy of the code from mizerExperimental, included here for easier
# modification before it is merged back into mizerExperimental later.

#' Plot the sources of external, predation and fishing mortality
#' per species and size

#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The name of the predator species for which to plot the
#'   mortality.
#' @param proportion A boolean value that determines whether values should be
#'   displayed as proportions from 0 to 1 or with their actual values. Default
#'   is TRUE.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#' @param ... Other arguments (currently unused)
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'w', 'value', 'Cause', 'Species' is returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' plotDeath(NS_params, species = "Cod")
#'
#' # Returning the data frame
#' fr <- plotDeath(NS_params, species = "Cod", return_data = TRUE)
#' str(fr)
#' }
plotDeathX <- function(object, species = NULL, proportion = TRUE,
                      return_data = FALSE,
                      xtrans = c("identity", "log10"),
                      xvar = c("Length", "Weight")) {
    if (is(object, "MizerSim")) {
        params <- object@params
        params <- setInitialValues(params, object)
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
    }
    if (!"External" %in% names(getColours(params))) {
        params <- setColours(params, c("External" = "grey"))
    }
    xtrans <- match.arg(xtrans)
    xvar <- match.arg(xvar)
    species <- valid_species_arg(params, species)

    pred_rate <- getPredRate(params)
    f_mort <- getFMort(params)
    mort <- getMort(params)
    plot_dat <- NULL
    for (iSpecies in species) {
        sps <- species_params(params)[iSpecies, ]
        fish_idx_full <- (params@w_full >= sps$w_min) &
            (params@w_full <= sps$w_max)
        fish_idx <- (params@w >= sps$w_min) &
            (params@w <= sps$w_max)
        w <- params@w[fish_idx]
        predation <- params@interaction[, iSpecies] *
            pred_rate[, fish_idx_full]
        fishing <- f_mort[iSpecies, fish_idx]
        external <- ext_mort(params)[iSpecies, fish_idx]
        total <- mort[iSpecies, fish_idx]
        ylab <- "Death rate [1/year]"
        if (proportion) {
            predation <- predation / rep(total, each = dim(predation)[[1]])
            external <- external / total
            fishing <- fishing / total
            ylab <- "Proportion of all death"
        }
        # Make data.frame for plot
        if (xvar == "Length") {
            size <- (w / sps$a) ^ (1 / sps$b)
        } else {
            size <- w
        }
        plot_dat <-
            rbind(
                plot_dat,
                data.frame(
                    size = size,
                    value = external,
                    Cause = "External",
                    Prey = iSpecies
                ),
                data.frame(
                    size = size,
                    value = fishing,
                    Cause = "Fishing",
                    Prey = iSpecies
                ),
                data.frame(
                    size = rep(size,
                        each = dim(predation)[[1]]
                    ),
                    value = c(predation),
                    Cause = params@species_params$species,
                    Prey = iSpecies
                )
            )
    }

    if (return_data) {
        return(plot_dat)
    }

    plotDataFrame(plot_dat, params,
        style = "area",
        wrap_var = "Prey", wrap_scale = "free_x",
        xlab = ifelse(xvar == "Length", "Size [cm]", "Size [g]"),
        ylab = ylab,
        xtrans = xtrans
    )
}


#' @rdname plotDeathX
#' @export
plotlyDeathX <- function(object,
                        species = NULL,
                        proportion = TRUE,
                        ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotDeathX", argg),
        tooltip = c("value", "Cause", "size")
    )
}
