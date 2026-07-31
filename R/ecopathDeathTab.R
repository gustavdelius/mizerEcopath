#' Ecopath death tab
#'
#' This tab provides diagnostic plots for mortality and production. It
#' allows the user to explore:
#' *   **Mortality density**: Breakdown of various mortality sources by size
#'     and species.
#' *   **Production inputs**: Numeric fields for editing the observed production
#'     (`production_observed`) and the production weight (`production_weight`)
#'     for the selected species. The `production_weight` sets how strongly the
#'     deviation of the model production from the observed production is
#'     penalised in `matchCatch()`.
#' *   **Production comparison**: The model production for the selected species,
#'     alongside the observed production and the resulting contribution to the
#'     loss function.
#' *   **Production vs Species**: Comparison of biomass produced across
#'     all species in the model.
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param logs Environment holding the log of steady states.
#' @param diet A diet matrix to match the diet of the model to.
#' @param ... Unused
#' @family gadget tabs
#' @export
ecopathDeathTab <- function(input, output, session, params, logs,
                     diet = NULL, ...) {
    # Plot mortality
    output$plot_mort <- renderPlotly({
        req(input$sp)
        p <- params()
        if (!is.null(diet)) {
            # p <- setFeedingLevels(params = p, f = 0.6, f_c = 0.2)
            p <- matchDiet(p, diet)
        }
        plotlyDeathX(p, species = input$sp,
                     proportion = input$death_prop == "Proportion",
                     xtrans = input$death_xtrans,
                     xvar = input$death_xvar)
    })

    # Input fields for observed production and production weight ----
    output$ecopath_production_inputs <- renderUI({
        p <- isolate(params())
        sp <- req(input$sp)
        spp <- p@species_params
        if (!sp %in% rownames(spp)) return(NULL)
        prod_obs <- if ("production_observed" %in% names(spp)) {
            spp[sp, "production_observed"]
        } else {
            NA_real_
        }
        prod_wt <- if ("production_weight" %in% names(spp)) {
            spp[sp, "production_weight"]
        } else {
            1
        }
        tagList(
            div(style = "display:inline-block; margin-right: 20px;",
                numericInput("production_observed",
                             paste0("Observed production for ", sp,
                                    " [g/year]"),
                             value = prod_obs)),
            div(style = "display:inline-block",
                numericInput("production_weight",
                             paste0("Production weight for ", sp),
                             value = prod_wt, min = 0, step = 0.1))
        )
    })

    # Adjust observed production ----
    observeEvent(input$production_observed, {
        p <- params()
        sp <- req(input$sp)
        if (!sp %in% rownames(p@species_params)) return()
        # `recalculate = FALSE` records the observation, so that a later
        # recalculation does not undo it, without recalculating any rates.
        sp_new <- p@species_params
        sp_new[sp, "production_observed"] <- input$production_observed
        species_params(p, recalculate = FALSE) <- sp_new
        tuneParams_update_params(p, params)
    }, ignoreInit = TRUE)

    # Adjust production weight ----
    observeEvent(input$production_weight, {
        p <- params()
        sp <- req(input$sp)
        if (!sp %in% rownames(p@species_params)) return()
        # `recalculate = FALSE` records the new weight, so that a later
        # recalculation does not undo it, without recalculating any rates.
        sp_new <- p@species_params
        if (!"production_weight" %in% names(sp_new)) {
            sp_new$production_weight <- 1
        }
        sp_new[sp, "production_weight"] <- input$production_weight
        species_params(p, recalculate = FALSE) <- sp_new
        tuneParams_update_params(p, params)
    }, ignoreInit = TRUE)

    # Comparison of model production, observed production and loss contribution
    output$ecopath_production_compare <- renderText({
        p <- params()
        sp <- req(input$sp)
        spp <- p@species_params
        if (!sp %in% rownames(spp)) return("")
        # The objective penalises the somatic production, so use that here.
        model_production <- getSomaticProduction(p)[[sp]]
        parts <- paste0("Model production: ", signif(model_production, 4),
                        " g/year")
        prod_obs <- if ("production_observed" %in% names(spp)) {
            spp[sp, "production_observed"]
        } else {
            NA_real_
        }
        if (!is.na(prod_obs) && prod_obs > 0) {
            prod_wt <- if ("production_weight" %in% names(spp)) {
                spp[sp, "production_weight"]
            } else {
                NA_real_
            }
            contrib <- prod_wt * (log(model_production / prod_obs))^2
            parts <- paste0(parts,
                            "   |   Observed production: ",
                            signif(prod_obs, 4),
                            " g/year   |   Loss contribution: ",
                            signif(contrib, 4))
        }
        parts
    })

    # Plot Production
    output$plot_prod <- renderPlotly({
        plotProductionVsSpecies(params())
    })
}

#' @rdname ecopathDeathTab
ecopathDeathTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_mort"),
        fluidRow(
            column(
                width = 4,
                radioButtons("death_prop", "Show:",
                             choices = c("Rate", "Proportion"),
                             selected = "Rate",
                             inline = TRUE)
            ),
            column(
                width = 4,
                radioButtons("death_xvar", "Show Size as:",
                             choices = c("Length", "Weight"),
                             selected = "Length", inline = TRUE)
            ),
            column(
                width = 4,
                radioButtons("death_xtrans", "x-axis scale:",
                             choices = c("log10", "identity"),
                             selected = "identity", inline = TRUE)
            )
        ),
        uiOutput("ecopath_production_inputs"),
        textOutput("ecopath_production_compare"),
        plotlyOutput("plot_prod")
    )
}

ecopathDeathTabTitle <- "Death"
