#' Ecopath spectra tab
#'
#' This tab displays the size spectra of the forground and background
#' species. It provides controls for:
#' *   **Scaling the background**: Adjusting the resource and background
#'     species' biomass.
#' *   **Adjusting/removing background species**: Fine-tuning the
#'     background community to align with the foreground.
#'
#' @inheritParams ecopathDeathTab
#' @family gadget tabs
#' @export
ecopathSpectraTab <- function(input, output, session,
                       params, logs, trigger_update, ...) {

    ## Plot spectra ####
    output$plotSpectra <- renderPlotly({
        power <- 2
        plot <- plotSpectra(params(), power = power, highlight = input$sp,
                            total = FALSE,
                            resource = FALSE) +
            theme(text = element_text(size = 12))
        ggplotly(plot, tooltip = c("Species", "w", "value"))
    })

    ## Scale ####
    observeEvent(input$scale_bkgrd_by, {
        if (input$scale_bkgrd_by == 1) return(1)
        tryCatch({
            p <- scaleDownBackground(params(), input$scale_bkgrd_by)
            tuneParams_update_params(p, params)
        }, error = error_fun)
        updateSliderInput(session, "scale_bkgrd_by", value = 1)
    })

    ## Retune background ####
    observeEvent(input$retune_background, {
        p <- adjustBackgroundSpecies(params())
        tuneParams_update_params(p, params)
    })

    ## Remove background ####
    observeEvent(input$remove_background, {
        p <- removeBackgroundSpecies(params())
        tuneParams_update_params(p, params)
    })
}

#' @rdname ecopathSpectraTab
ecopathSpectraTabUI <- function(params, help = TRUE, ...) {
    p <- isolate(params())

    tl <- tagList(plotlyOutput("plotSpectra"),
                  div(style = "display:inline-block;vertical-align:middle; width: 300px;",
                      popify(sliderInput("scale_bkgrd_by",
                                         "Scale background down by a factor of:",
                                         value = 1, min = 0.5, max = 2, step = 0.1),
                             title = "Scaling the background",
                             content = "You can scale down the background in which the fish find themselves (the resource and any background species that your model may contain).")),
    )

    if (anyNA(p@A)) {
        tl <- tagList(tl,
                      popify(actionButton("retune_background", "Adj bs"),
                             title = "Adjust background species",
                             content = "Adjust the biomasses of the background species."),
                      popify(actionButton("remove_background", "Rem bs"),
                             title = "Remove background species",
                             content = "Remove all background species."))
    }
    if (help) {
        tl <- tagList(tl,
                      h1("Size spectra"),
                      p("This tab shows the biomass size spectra of the individual fish species and of the resource."))
    }
    tl
}

ecopathSpectraTabTitle <- "Spectra"

