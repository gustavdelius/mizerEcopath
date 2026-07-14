#' Ecopath reproduction tab
#'
#' This tab displays parameters related to species maturity and
#' reproductive success. It includes:
#' *   **Reproductive success plot**: Compares the realized recruitment
#'     (RDD) with the potential recruitment (RDI) for each species.
#'     A value of 1 (red line) indicates that recruitment is
#'     primarily density-independent.
#' *   **Maturity ogive plot**: Shows the proportion of energy allocated
#'     to reproduction (psi) vs size for the selected species.
#'
#' @inheritParams ecopathDeathTab
#' @family gadget tabs
#' @export
ecopathReproTab <- function(input, output, session, params, logs, ...) {
    # erepro plot ----
    output$plot_erepro <- renderPlotly({
        p <- params()
        foreground <- !is.na(p@A)
        rdi <- getRDI(p)[foreground]
        rdd <- getRDD(p)[foreground]
        repro_success <- p@species_params$erepro[foreground] * rdd / rdi
        df <- data.frame(Species = factor(p@species_params$species[foreground],
                                          levels = p@species_params$species[foreground]),
                         value = repro_success)
        ggplot(df, aes(x = Species, y = value)) +
            geom_col() + geom_hline(yintercept = 1, color = "red") +
            scale_y_log10(name = "Reproductive success") +
            theme(text = element_text(size = 12)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
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
