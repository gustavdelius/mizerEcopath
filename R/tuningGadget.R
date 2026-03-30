#' Launch a customisable shiny gadget for tuning mizer parameters
#'
#' This is the shared framework used by [tuneParams()] and, in the
#' mizerEcopath package, by `tuneEcopath()`. Advanced users can also call it
#' directly to build custom tuning gadgets.
#'
#' The function opens a shiny gadget with a sidebar of parameter controls and a
#' main panel of diagnostic tabs.
#'
#' # Undo functionality
#'
#' The gadget keeps a log of all steady states you create while working with the
#' gadget. You can go back to the last steady state by hitting the "Undo"
#' button. You can go back an arbitrary number of states and also go forward
#' again. There is also a button to go right back to the initial steady state.
#'
#' When you leave the gadget by hitting the "Return" button, this log is
#' cleared. If you stop the gadget from RStudio by hitting the "Stop" button,
#' then the log is left behind. You can then restart the gadget by calling
#' `tuneParams()` without a `params` argument and it will re-instate the states
#' from the log.
#'
#' # Customisation
#'
#' You can customise which functionality is included in the app via the
#' `controls` and `tabs` arguments. You can remove some of the controls and tabs
#' by providing shorter lists to those arguments. You can also add your own
#' controls and tabs.
#'
#' For an entry "foo" in the `controls` list there needs to be a function
#' "fooControlUI" that defines the input elements and a function "fooControl"
#' that processes those inputs to change the params object.
#'
#' For any entry "Foo" or "foo" in the `tabs` list there needs to be a function
#' "fooTabUI" that defines the tab layout and a function "fooTab" that
#' calculates the outputs to be displayed on the tab.
#'
#' @param params MizerParams object to tune. If missing, the gadget tries to
#'   recover information from log files left over from aborted previous runs.
#' @param controls A character vector of names of input parameter control
#'   sections that should be displayed in the sidebar.
#' @param tabs A character vector of names of the tabs that should be displayed
#'   in the main section.
#' @param match_choices A character vector of the available match options that
#'   appear as checkboxes.
#' @param match_selected A character vector indicating which match options are
#'   selected by default.
#' @param action_label The HTML label for the main action button (e.g.,
#'   `HTML("<u>s</u>teady")`).
#' @param action_tooltip Tooltip text for the main action button.
#' @param action_key JavaScript key code for the keyboard shortcut that
#'   triggers the action button (e.g., 83 for 's', 77 for 'm').
#' @param action_handler A function that is called when the action button is
#'   pressed. It receives the arguments `(p, input, session, params,
#'   params_old, logs, flags, trigger_update, ...)`.
#' @param prepare_params_hook An optional function `(p) -> p` called after the
#'   shared parameter preparation to apply package-specific defaults.
#' @param finalise_params_hook An optional function `(p) -> p` called after the
#'   shared parameter finalisation to apply package-specific cleanup.
#' @param update_species_hook An optional function `(sp, p) -> p` called inside
#'   [tuneParams_update_species()] after `steadySingleSpecies()` to apply
#'   package-specific post-processing (e.g., `matchBiomasses()`).
#' @param tab_extra_args A named list of extra arguments to pass to tab server
#'   functions (e.g., `list(catch = catch, diet = diet)`).
#' @param preserve Specifies whether the `reproduction_level` should be
#'   preserved or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro` (Default). See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#' @param return_app Boolean. For testing purposes only. When set to TRUE will
#'   return a shinyApp object instead of running the gadget.
#' @param ... Other params needed by individual tabs.
#'
#' @return The tuned MizerParams object
#' @md
#' @import shinyBS
#' @export
tuningGadget <- function(params,
                         controls,
                         tabs,
                         match_choices,
                         match_selected,
                         action_label,
                         action_tooltip,
                         action_key,
                         action_handler,
                         prepare_params_hook = NULL,
                         finalise_params_hook = NULL,
                         update_species_hook = NULL,
                         tab_extra_args = list(),
                         preserve = c("erepro", "reproduction_level", "R_max"),
                         return_app = FALSE,
                         ...) {

    # Define some local variables to avoid "no visible bindings for global
    # variable" warnings in CMD check
    wpredator <- wprey <- Nprey <- weight_kernel <- L_inf <-
        Legend <- w_mat <- erepro <- Type <- Abundance <- Catch <-
        Kernel <- Numbers <- Cause <- psi <- Predator <- Density <- NULL

    # Flags to skip certain observers ----
    flags <- new.env()

    # Prepare logs for undo/redo functionality ----
    logs <- new.env()
    logs$files <- vector(mode = "character")
    logs$idx <- 0

    if (missing(params)) {
        # Try to recover old log files ----
        logs$files <- sort(list.files(path = tempdir(),
                                pattern = "mizer_params_...._.._.._at_.._.._..\\.rds",
                                full.names = TRUE))
        logs$idx <- length(logs$files)
        if (logs$idx == 0) {
            stop("You need to specify a MizerParams object. ",
                 "There are no temporary parameter files to recover.")
        }
        p <- readRDS(logs$files[logs$idx])
    } else {
        p <- params # because I use `params` as a reactive value later
        validObject(p)

        # Add the info that should be preserved to the species_params for later
        # recall. This is not needed when a params object is restored from the
        # logs because its already included there.
        preserve <- match.arg(preserve)
        if (preserve == "reproduction_level") {
            p@species_params$tuneParams_old_repro_level <-
                getReproductionLevel(p)
        }
        if (preserve == "R_max") {
            p@species_params$tuneParams_old_R_max <- p@species_params$R_max
        } else {
            p@species_params$tuneParams_old_erepro <- p@species_params$erepro
        }
        p <- prepare_params(p, prepare_params_hook)
    }

    # Store the update_species_hook in an environment accessible to helpers
    tuning_env$update_species_hook <- update_species_hook

    # Store the finalise_params_hook
    local_finalise <- function(p) {
        p <- finalise_params(p)
        if (!is.null(finalise_params_hook)) {
            p <- finalise_params_hook(p)
        }
        p
    }

    # User interface ----
    ui <- fluidPage(
        theme = bslib::bs_theme(version = 4, bootswatch = "cerulean"),
        shinyjs::useShinyjs(),
        prompter::use_prompt(),
        introjsUI(),
        tags$script(HTML(paste0("$(function(){
          $(document).keydown(function(e) {
          if (e.which == ", action_key, ") {
            $('#sp_steady').click()
          }
          if (e.which == 78) {
            $('#next_sp').click()
          }
          if (e.which == 80) {
            $('#previous_sp').click()
          }
        });});"))),
        tags$head(
            tags$style(HTML(".center{float:center;}"))),

        sidebarLayout(

            ## Sidebar ####
            sidebarPanel(
                width = 3,
                introBox(
                    prompter::add_prompt(
                        actionButton("help", "Help"),
                        message = "Start the introductory instructions",
                        position = "right"),
                    prompter::add_prompt(
                        downloadButton("download_params", ""),
                        message = "Download the current params object",
                        position = "right"),
                    prompter::add_prompt(
                        actionButton("done", "Return", icon = icon("check"),
                                     onclick = "setTimeout(function(){window.close();},500);"),
                        message = "Return the current params objects to R"),
                    data.step = 8,
                    data.position = "right",
                    data.intro = "At any point you can press the download button to save the current state of the params object. When you press the 'Return' button, the gadget will close and the current params object will be returned. The undo log will be cleared."
                ),
                introBox(
                    prompter::add_prompt(
                        actionButton("sp_steady", action_label),
                        message = action_tooltip,
                        position = "right"),
                    prompter::add_prompt(
                        actionButton("undo_all", "", icon = icon("angles-left")),
                        message = "Undo all changes"),
                    prompter::add_prompt(
                        actionButton("undo", "", icon = icon("angle-left")),
                        message = "Go back to previous steady state"),
                    prompter::add_prompt(
                        actionButton("redo", "", icon = icon("angle-right")),
                        message = "Go forward to next steady state"),
                    data.step = 5,
                    data.intro = "Each time you change a parameter, the spectrum of the selected species is immediately recalculated. However to calculate the true multi-species steady state you have to press the action button or use the keyboard shortcut. Do this frequently, before changing the parameters too much. Otherwise there is the risk that the steady state can not be found any more. You can go backwards and forwards among the previously calculated steady states with the 'Undo All', 'Undo' and 'Redo' buttons.",
                    data.position = "right"
                ),

                introBox(
                    prompter::add_prompt(
                        checkboxGroupInput("match", "Match:",
                                           choices = match_choices,
                                           selected = match_selected,
                                           inline = TRUE),
                        message = "Choose quantities to match to observations automatically"),
                    data.step = 6,
                    data.position = "right",
                    data.intro = "Here you can specify that each time you hit the action button the selected quantities are matched to their observed values. This does of course not mean that a perfect match will be achieved in the steady state. But usually each time you hit the button the match will improve."
                ),

                introBox(
                    prompter::add_prompt(
                        uiOutput("sp_sel"),
                        message = "Select target species",
                        position = "right"),
                    data.step = 2,
                    data.position = "right",
                    data.intro = "Here you select the species whose parameters you want to change or whose properties you want to concentrate on."
                ),
                introBox(
                    introBox(
                        # Add links to input sections
                        lapply(controls, function(section) {
                            list("->",
                                 tags$a(control_title(section),
                                        href = paste0("#", section)))
                        }),
                        data.step = 4,
                        data.position = "right",
                        data.intro = "There are many parameters, organised into sections. To avoid too much scrolling you can click on a link to jump to a section."),
                    tags$br(),
                    tags$div(id = "params",
                             uiOutput("sp_params")
                    ),
                    tags$head(tags$style(
                        type = 'text/css',
                        '#params { max-height: 60vh; overflow-y: auto; }'
                    )),
                    data.step = 3,
                    data.intro = "Here you find controls for changing model parameters. The controls for species-specific parameters are for the species you have chosen above. Many of the controls are sliders that you can move by dragging or by clicking. As you change parameters, the plots in the main panel will update immediately."
                )
            ),  # endsidebarpanel

            ## Main panel ####
            mainPanel(
                width = 9,
                introBox(uiOutput("tabs"),
                         data.step = 1,
                         data.intro = "This main panel has tabs that display various aspects of the steady state of your model. At the bottom of each tab you find text explanations for that tab. You may need to scroll down in the tab to see them. Individual components may show tooltips when you hover over them."
                )
            )  # end mainpanel
        )  # end sidebarlayout
    )

    server <- function(input, output, session) {
        hintjs(session)
        ## Create params object as a reactive value ####
        params <- reactiveVal()
        if (logs$idx == 0) {
            # Adding it to the logs allows us to get back to the initial state
            tuneParams_add_to_logs(logs, p, params)
        } else {
            # Otherwise p comes from recovered logs and only needs to be saved
            # in reactive value
            params(p)
        }

        # The abundances in params() will get updated using
        # `singleSpeciesSteady()` each time a parameter gets changed. It is
        # important that this updating always starts again from the previous
        # abundances. So we need to preserve those previous abundances. That
        # is what `params_old()` will be for.
        # The abundances in params_old() will only get changed when we run to
        # steady state or take a steady state from the logs. They will not be
        # changed in `tuneParams_update_species()`.
        params_old <- reactiveVal(p)

        if (logs$idx == length(logs$files)) {
            shinyjs::disable("redo")
        }
        if (logs$idx <= 1) {
            # There is nothing to undo if there is nothing older in the logs
            shinyjs::disable("undo")
            shinyjs::disable("undo_all")
        }

        # Define a reactive value for triggering an update of species sliders
        trigger_update <- reactiveVal(0)

        ## UI for side bar ####
        # Drop-down menu for selecting active species
        output$sp_sel <- renderUI({
            p <- isolate(params())
            species <- as.character(p@species_params$species[!is.na(p@A)])
            tagList(
                selectInput("sp", "Species to tune:", species),
                prompter::add_prompt(
                    actionButton("previous_sp", HTML("<u>p</u>revious")),
                    message = "Select previous species. Keyboard shortcut: p",
                    position = "right"),
                prompter::add_prompt(
                    actionButton("next_sp", HTML("<u>n</u>ext")),
                    message = "Select next species. Keyboard shortcut: n",
                    position = "right"))
            })
        # Sliders for the species parameters
        output$sp_params <- renderUI({
            # The parameter sliders get updated whenever the species selector
            # changes
            req(input$sp)
            # or when the trigger is set somewhere
            trigger_update()
            # but not each time the params change
            p <- isolate(params())

            lapply(controls,
                   function(section) {
                       do.call(paste0(section, "ControlUI"),
                               list(p = p, input = input))
                   })
        })

        # Serve controls ####
        for (section in controls) {
            fun <- paste0(section, "Control")
            substr(fun, 1, 1) <- tolower(substr(fun, 1, 1))
            do.call(fun, list(input = input,
                              output = output,
                              session = session,
                              params = params,
                              params_old = params_old,
                              flags = flags,
                              trigger_update = trigger_update))
        }

        ## UI for tabs ####
        output$tabs <- renderUI({
            tablist <- lapply(tabs, function(tab) {
                tabname <- tab_name(tab)
                tab_content <- div(
                    style = "max-height: 94vh; overflow-y: auto; overflow-x: hidden;",
                    do.call(paste0(tabname, "TabUI"),
                            list(params = params)))
                tabPanel(tab_title(tab), tab_content)
            })
            args <- c(id = "mainTabs", type = "tabs", tablist)
            do.call(tabsetPanel, args)
        })

        ## Serve tabs ####
        tab_args <- c(list(input = input,
                           output = output,
                           session = session,
                           params = params,
                           params_old = params_old,
                           logs = logs,
                           trigger_update = trigger_update),
                      tab_extra_args,
                      list(...))
        for (tab in tabs) {
            tabname <- tab_name(tab)
            fun <- paste0(tabname, "Tab")
            do.call(fun, tab_args)
        }

        # Help button ----
        observeEvent(
            input$help,
            introjs(session)
        )

        ## Action button handler ####
        # triggered by the main action button in sidebar
        observeEvent(input$sp_steady, {
            action_handler(p = params(),
                           input = input,
                           session = session,
                           params = params,
                           params_old = params_old,
                           logs = logs,
                           flags = flags,
                           trigger_update = trigger_update)
        })

        ## Previous ####
        observeEvent(input$previous_sp, {
            p <- params()
            all_species <- species_params(p)$species[!is.na(p@A)]
            no_sp <- length(all_species)
            idx <- which(all_species == input$sp)
            prev <- ((idx - 2) %% no_sp) + 1
            updateSelectInput(session, "sp",
                              selected = all_species[[prev]])
        })
        ## Next ####
        observeEvent(input$next_sp, {
            p <- params()
            all_species <- species_params(p)$species[!is.na(p@A)]
            no_sp <- length(all_species)
            idx <- which(all_species == input$sp)
            updateSelectInput(session, "sp",
                              selected = all_species[[(idx %% no_sp) + 1]])
        })

        ## Undo ####
        observeEvent(input$undo, {
            if (logs$idx < 1) stop("This should never happen")
            p_new <- readRDS(logs$files[logs$idx])
            p_old <- params()
            # if the params have not changed, go to the previous one,
            # if it exists
            if (is.null(attr(p_old, "changes")) && logs$idx > 1) {
                logs$idx <- logs$idx - 1
                shinyjs::enable("redo")
                p_new <- readRDS(logs$files[logs$idx])
            } else {
                shinyjs::disable("redo")
            }
            if (logs$idx <= 1) {
                shinyjs::disable("undo")
                shinyjs::disable("undo_all")
            }
            params(p_new)
            params_old(p_new)
            # Trigger an update of sliders
            # First clear all flags so that update of sliders does not trigger
            # any controls.
            rm(list = ls(flags), pos = flags)
            trigger_update(runif(1))
        })
        ## Redo ####
        observeEvent(input$redo, {
            if (logs$idx >= length(logs$files)) return()
            logs$idx <- logs$idx + 1
            p <- readRDS(logs$files[logs$idx])
            params(p)
            params_old(p)
            # Trigger an update of sliders
            rm(list = ls(flags), pos = flags)
            trigger_update(runif(1))
            shinyjs::enable("undo")
            shinyjs::enable("undo_all")
            if (logs$idx == length(logs$files)) shinyjs::disable("redo")
        })
        ## Undo All ####
        observeEvent(input$undo_all, {
            if (logs$idx > 1) shinyjs::enable("redo")
            shinyjs::disable("undo")
            shinyjs::disable("undo_all")
            logs$idx <- 1
            p <- readRDS(logs$files[logs$idx])
            params(p)
            params_old(p)
            # Trigger an update of sliders
            rm(list = ls(flags), pos = flags)
            trigger_update(runif(1))
        })

        ## Prepare for download of params object ####
        output$download_params <- downloadHandler(
            filename = "tuned_params.rds",
            content = function(file) {
                saveRDS(local_finalise(params()), file = file)
        })

        ## Return ####
        # When the user hits the "Return" button we want to clear the logs and
        # return with the latest params object
        observeEvent(input$done, {
            file.remove(logs$files)
            stopApp(local_finalise(params()))
        })

    } #the server

    if (return_app) {
        return(shinyApp(ui, server))
    }

    runGadget(ui, server, viewer = browserViewer())
}
