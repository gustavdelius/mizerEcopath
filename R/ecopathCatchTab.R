#' Ecopath catch tab
#'
#' This tab provides diagnostic plots for comparison with catch data:
#' *   **Catch size distribution**: Detailed size breakdown of
#'     yield for the selected species and gear.
#' *   **Size distribution weight**: A slider for the `catch_dist_weight` gear
#'     parameter of the selected species and gear. It sets how strongly the
#'     mismatch between the observed and the modelled catch size distribution is
#'     penalised in `matchCatch()`.
#' *   **Yield inputs**: Numeric fields for editing the observed yield
#'     (`yield_observed`) and the yield weight (`yield_weight`) for the selected
#'     species and gear. The `yield_weight` sets how strongly the deviation of
#'     the model yield from the observed yield is penalised in `matchCatch()`.
#' *   **Yield comparison**: The model yield for the selected species and gear,
#'     alongside the observed yield and the resulting contribution to the loss
#'     function.
#' *   **Total yield plot**: A bar-like plot showing circles (model)
#'     and squares (observed) for each species caught by the gear.
#'
#' @inheritParams ecopathDeathTab
#' @param catch Data frame holding binned observed catch data.
#' @family gadget tabs
#' @export
ecopathCatchTab <- function(input, output, session, params, logs,
                     catch = NULL, ...) {

    if (!is.null(catch)) {
        assert_that(
            is.data.frame(catch),
            "catch" %in% names(catch),
            "species" %in% names(catch),
            "gear" %in% names(catch),
            all(c("length", "dl") %in% names(catch)) |
                all(c("weight", "dw") %in% names(catch))
        )
    }

    # Catch density for selected species ----
    output$plotCatchDist <- renderPlotly({
        plotlyYieldVsSize(params(), species = req(input$sp),
                          gear = input$gear,
                          catch = catch, x_var = "Length")
    })

    # Slider for the weight of the catch size distribution penalty ----
    output$ecopath_catch_dist_weight <- renderUI({
        p <- isolate(params())
        sp <- req(input$sp)
        gear <- req(input$gear)
        spgear <- paste(sp, gear, sep = ", ")
        gp <- p@gear_params
        if (!spgear %in% rownames(gp)) return(NULL)
        wt <- if ("catch_dist_weight" %in% names(gp)) {
            gp[spgear, "catch_dist_weight"]
        } else {
            1
        }
        if (is.na(wt)) wt <- 0
        sliderInput("catch_dist_weight",
                    paste0("Size distribution weight for ", spgear),
                    value = wt, min = 0,
                    max = signif(max(5, wt * 2), 2), step = 0.05)
    })

    # Adjust the weight of the catch size distribution penalty ----
    observeEvent(input$catch_dist_weight, {
        p <- params()
        spgear <- paste(req(input$sp), req(input$gear), sep = ", ")
        if (!spgear %in% rownames(p@gear_params)) return()
        if (!"catch_dist_weight" %in% names(p@gear_params)) {
            p@gear_params$catch_dist_weight <- 1
        }
        p@gear_params[spgear, "catch_dist_weight"] <- input$catch_dist_weight
        tuneParams_update_params(p, params)
    }, ignoreInit = TRUE)

    # Input fields for observed yield and yield weight ----
    output$ecopath_yield_inputs <- renderUI({
        p <- isolate(params())
        sp <- req(input$sp)
        gear <- req(input$gear)
        spgear <- paste(sp, gear, sep = ", ")
        gp <- p@gear_params
        if (!spgear %in% rownames(gp)) return(NULL)
        yield_obs <- if ("yield_observed" %in% names(gp)) {
            gp[spgear, "yield_observed"]
        } else {
            NA_real_
        }
        yield_wt <- if ("yield_weight" %in% names(gp)) {
            gp[spgear, "yield_weight"]
        } else {
            1
        }
        tagList(
            div(style = "display:inline-block; margin-right: 20px;",
                numericInput("yield_observed",
                             paste0("Observed yield for ", spgear, " [g/year]"),
                             value = yield_obs)),
            div(style = "display:inline-block",
                numericInput("yield_weight",
                             paste0("Yield weight for ", spgear),
                             value = yield_wt, min = 0, step = 0.1))
        )
    })

    # Adjust observed yield ----
    observeEvent(input$yield_observed, {
        p <- params()
        spgear <- paste(req(input$sp), req(input$gear), sep = ", ")
        if (!spgear %in% rownames(p@gear_params)) return()
        p@gear_params[spgear, "yield_observed"] <- input$yield_observed
        tuneParams_update_params(p, params)
    }, ignoreInit = TRUE)

    # Adjust yield weight ----
    observeEvent(input$yield_weight, {
        p <- params()
        spgear <- paste(req(input$sp), req(input$gear), sep = ", ")
        if (!spgear %in% rownames(p@gear_params)) return()
        if (!"yield_weight" %in% names(p@gear_params)) {
            p@gear_params$yield_weight <- 1
        }
        p@gear_params[spgear, "yield_weight"] <- input$yield_weight
        tuneParams_update_params(p, params)
    }, ignoreInit = TRUE)

    # Comparison of model yield, observed yield and loss contribution ----
    output$ecopath_yield_compare <- renderText({
        p <- params()
        sp <- req(input$sp)
        gear <- req(input$gear)
        spgear <- paste(sp, gear, sep = ", ")
        gp <- p@gear_params
        if (!spgear %in% rownames(gp)) return("")
        model_yield <- sum(p@initial_n[sp, ] * p@w * p@dw *
                               getFMortGear(p)[gear, sp, ])
        parts <- paste0("Model yield: ", signif(model_yield, 4), " g/year")
        yield_obs <- if ("yield_observed" %in% names(gp)) {
            gp[spgear, "yield_observed"]
        } else {
            NA_real_
        }
        if (!is.na(yield_obs) && yield_obs > 0) {
            yield_wt <- if ("yield_weight" %in% names(gp)) {
                gp[spgear, "yield_weight"]
            } else {
                NA_real_
            }
            contrib <- yield_wt * (log(model_yield / yield_obs))^2
            parts <- paste0(parts,
                            "   |   Observed yield: ", signif(yield_obs, 4),
                            " g/year   |   Loss contribution: ",
                            signif(contrib, 4))
        }
        parts
    })

    # Total yield by species ----
    output$plotTotalYield <- renderPlot({
        plotYieldVsSpecies(params(), gear = input$gear) +
            theme(text = element_text(size = 18))
    })
}

#' @rdname ecopathCatchTab
ecopathCatchTabUI <- function(...) {
    tagList(
        plotlyOutput("plotCatchDist"),
        uiOutput("ecopath_catch_dist_weight"),
        uiOutput("ecopath_yield_inputs"),
        textOutput("ecopath_yield_compare"),
        plotOutput("plotTotalYield")
    )
}

ecopathCatchTabTitle <- "Catch"
