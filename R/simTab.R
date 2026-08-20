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
#' @param params_old Reactive value holding the MizerParams object from before
#'   the most recent update.
#' @param method The numerical method used for the projection to steady state,
#'   passed on to [mizer::project()]. Supplied by [tuneParams()].
#' @family gadget tabs
#' @export
simTab <- function(input, output, session, params, params_old, logs,
                   method = "euler", ...) {

    ## Plot run to steady ####
    output$plot_sim <- renderPlotly({
        sim <- tuneParams_run_steady(params(), return_sim = TRUE,
                                     params = params, params_old = params_old,
                                     logs = logs, method = method,
                                     session = session, input = input)
        plotlyBiomass(sim)
    })
}
