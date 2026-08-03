#' Reproduction tab for tuning gadget
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
#' @inheritParams biomassTab
#' @family gadget tabs
#' @export
reproTab <- function(input, output, session, params, logs, ...) {
    # Reproductive efficiency plot ----
    output$plot_erepro <- renderPlotly({
        plotReproductiveEfficiency(params())
    })

    # Plot psi ----
    output$plot_psi <- renderPlotly({
        plotHover(psi(params()), species = req(input$sp),
                  wlim = c(1, NA), log = "")
    })
}

#' @rdname reproTab
reproTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_erepro"),
        plotlyOutput("plot_psi")
    )
}
