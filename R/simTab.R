#' @rdname simTab
simTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_sim")
    )
}

#' Simulation tab for tuning gadget
#'
#' This tab shows the results of running the model to steady state. It
#' displays the biomass dynamics during the simulation, helping to
#' identify if the system reaches a stable equilibrium.
#'
#' @inheritParams spectraTab
#' @family gadget tabs
#' @export
simTab <- function(input, output, session, params, params_old, logs, ...) {
    
    ## Plot run to steady ####
    output$plot_sim <- renderPlotly({
        sim <- tuneParams_run_steady(params(), return_sim = TRUE,
                                     params = params, params_old = params_old,
                                     logs = logs,
                                     session = session, input = input)
        plotlyBiomass(sim)
    })
}
