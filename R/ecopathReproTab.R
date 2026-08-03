#' Ecopath reproduction tab
#'
#' This tab displays parameters related to species maturity and
#' reproductive success. It includes:
#' *   **Reproductive efficiency plot**: Shows which fraction of the energy
#'     that the species invests into reproduction ends up as egg biomass,
#'     both before and after density dependence. Values above the red line
#'     at 1 are energetically impossible. See
#'     [plotReproductiveEfficiency()].
#' *   **Maturity ogive plot**: Shows the proportion of energy allocated
#'     to reproduction (psi) vs size for the selected species.
#'
#' @inheritParams ecopathDeathTab
#' @family gadget tabs
#' @export
ecopathReproTab <- function(input, output, session, params, logs, ...) {
    # erepro plot ----
    output$plot_erepro <- renderPlotly({
        plotReproductiveEfficiency(params())
    })

    # Plot psi ----
    output$plot_psi <- renderPlotly({
        p <- params()
        species <- input$sp
        plotHover(psi(p), species = species,
                  size_axis = "l", llim = c(1, NA), log = "")
    })
}

#' @rdname ecopathReproTab
ecopathReproTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_erepro"),
        plotlyOutput("plot_psi")
    )
}

ecopathReproTabTitle <- "Repro"
