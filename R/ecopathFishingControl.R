#' Ecopath-specific fishing control
#'
#' This control allows for fine-tuning of fishing parameters for each gear
#' capturing a given species. It allows adjusting:
#' *   **Catchability**: A global scaling factor for the fishing mortality of
#'     the gear.
#' *   **Selectivity**: Depending on the selection function (knife-edge,
#'     sigmoid, or double-sigmoid), users can adjust the size at which
#'     individuals become vulnerable to the gear.
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
#' @family gadget controls
#' @export
ecopathFishingControl <- function(input, output, session, params, params_old,
                           flags, ...) {

    # Persistent storage: species x gear x parameter
    gear_vals <- reactiveValues()

    get_val <- function(sp, gear, name, default) {
        if (!is.null(gear_vals[[sp]]) &&
            !is.null(gear_vals[[sp]][[gear]]) &&
            !is.null(gear_vals[[sp]][[gear]][[name]])) {
            gear_vals[[sp]][[gear]][[name]]
        } else {
            default
        }
    }

    set_val <- function(sp, gear, name, value) {
        if (is.null(gear_vals[[sp]])) gear_vals[[sp]] <- list()
        if (is.null(gear_vals[[sp]][[gear]])) gear_vals[[sp]][[gear]] <- list()
        gear_vals[[sp]][[gear]][[name]] <- value
    }


    # RESTORE SLIDERS when species or gear changes
    observeEvent(
        list(input$sp, input$gear),
        {
            req(input$sp, input$gear)

            sp   <- input$sp
            gear <- input$gear
            p    <- params()

            gp_idx <- which(p@gear_params$species == sp &
                                p@gear_params$gear == gear)
            if (length(gp_idx) == 0) return()

            sel_func <- p@gear_params[gp_idx, "sel_func"]

            # ---- Catchability ----
            catch_log <- get_val(
                sp, gear, "catchability",
                log10(p@gear_params[gp_idx, "catchability"])
            )

            updateSliderInput(
                session, "catchability",
                value = catch_log
            )

            # ---- Knife edge ----
            if (sel_func == "knife_edge") {
                ke <- get_val(
                    sp, gear, "knife_edge_size",
                    p@gear_params[gp_idx, "knife_edge_size"]
                )
                updateSliderInput(
                    session, "knife_edge_size",
                    value = ke,
                    max = signif(ke * 2, 2)
                )
            }

            # ---- Sigmoid ----
            if (sel_func %in% c("sigmoid_length", "double_sigmoid_length")) {
                l50 <- get_val(sp, gear, "l50", p@gear_params[gp_idx, "l50"])
                ldiff <- get_val(
                    sp, gear, "ldiff",
                    p@gear_params[gp_idx, "l50"] -
                        p@gear_params[gp_idx, "l25"]
                )

                updateSliderInput(session, "l50",
                                  value = l50,
                                  max = signif(l50 * 2, 2))
                updateSliderInput(session, "ldiff",
                                  value = ldiff,
                                  max = signif(ldiff * 2, 2))
            }

            # ---- Double sigmoid (right) ----
            if (sel_func == "double_sigmoid_length") {
                l50r <- get_val(
                    sp, gear, "l50_right",
                    p@gear_params[gp_idx, "l50_right"]
                )
                ldiffr <- get_val(
                    sp, gear, "ldiff_right",
                    p@gear_params[gp_idx, "l25_right"] -
                        p@gear_params[gp_idx, "l50_right"]
                )

                updateSliderInput(session, "l50_right",
                                  value = l50r,
                                  max = signif(l50r * 2, 2))
                updateSliderInput(session, "ldiff_right",
                                  value = ldiffr,
                                  max = signif(ldiffr * 2, 2))
            }
        },
        ignoreInit = TRUE
    )


    # WRITE TO PARAMS when sliders change
    observeEvent(
        list(
            input$catchability,
            input$knife_edge_size,
            input$l50,
            input$ldiff,
            input$l50_right,
            input$ldiff_right
        ),
        {
            req(input$sp, input$gear)

            sp   <- input$sp
            gear <- input$gear
            p    <- isolate(params())

            gp_idx <- which(p@gear_params$species == sp &
                                p@gear_params$gear == gear)
            if (length(gp_idx) == 0) return()

            sel_func <- p@gear_params[gp_idx, "sel_func"]

            # ---- Catchability ----
            set_val(sp, gear, "catchability", input$catchability)
            p@gear_params[gp_idx, "catchability"] <- 10^input$catchability

            # ---- Knife edge ----
            if (sel_func == "knife_edge" && !is.null(input$knife_edge_size)) {
                set_val(sp, gear, "knife_edge_size", input$knife_edge_size)
                p@gear_params[gp_idx, "knife_edge_size"] <- input$knife_edge_size
            }

            # ---- Sigmoid ----
            if (sel_func %in% c("sigmoid_length", "double_sigmoid_length") &&
                !is.null(input$l50) && !is.null(input$ldiff)) {

                set_val(sp, gear, "l50", input$l50)
                set_val(sp, gear, "ldiff", input$ldiff)

                p@gear_params[gp_idx, "l50"] <- input$l50
                p@gear_params[gp_idx, "l25"] <- input$l50 - input$ldiff
            }

            # ---- Double sigmoid (right) ----
            if (sel_func == "double_sigmoid_length" &&
                !is.null(input$l50_right) && !is.null(input$ldiff_right)) {

                set_val(sp, gear, "l50_right", input$l50_right)
                set_val(sp, gear, "ldiff_right", input$ldiff_right)

                p@gear_params[gp_idx, "l50_right"] <- input$l50_right
                p@gear_params[gp_idx, "l25_right"] <- input$l50_right + input$ldiff_right
            }

            p <- setFishing(p)
            tuneParams_update_species(sp, p, params, params_old)
        },
        ignoreInit = TRUE
    )
}

#' @rdname ecopathFishingControl
#'
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the gear parameters
ecopathFishingControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    gp <- params@gear_params[params@gear_params$species == sp$species, ]
    if (nrow(gp) == 0) { # Species not selected by any gears
        return(tagList())
    }
    gears <- as.character(gp$gear)
    if (is.null(input$gear) || !(input$gear %in% gears)) {
        gear <- gears[[1]]
    } else {
        gear <- input$gear
    }

    gp <- gp[gp$gear == gear, ]
    l1 <- list(tags$h3(tags$a(id = "fishing"), "Fishing"),
               selectInput("gear", "Gear to tune:", gears,
                           selected = gear),
               sliderInput("catchability", "Catchability",
                           value = log10(gp$catchability),
                           min = -13,
                           max = 4,
                           step=0.01)
    )

    if (gp$sel_func == "knife_edge") {
        l1 <- c(l1, list(
            sliderInput("knife_edge_size", "knife_edge_size",
                        value = gp$knife_edge_size,
                        min = 1,
                        max = signif(gp$knife_edge_size * 2, 2),
                        step = 0.1)))
    } else if (gp$sel_func == "sigmoid_length") {
        l1 <- c(l1, list(
            sliderInput("l50", "L50",
                        value = gp$l50,
                        min = 1,
                        max = signif(gp$l50 * 2, 2),
                        step = 0.1),
            sliderInput("ldiff", "L50-L25",
                        value = gp$l50 - gp$l25,
                        min = 0.1,
                        max = signif((gp$l50 - gp$l25) * 2, 2),
                        step = 0.1)))
    } else if (gp$sel_func == "double_sigmoid_length") {
        l1 <- c(l1, list(
            sliderInput("l50", "L50",
                        value = gp$l50,
                        min = 1,
                        max = signif(gp$l50 * 2, 2),
                        step = 0.1),
            sliderInput("ldiff", "L50-L25",
                        value = gp$l50 - gp$l25,
                        min = 0.1,
                        max = signif((gp$l50 - gp$l25) * 2, 2),
                        step = 0.1),
            sliderInput("l50_right", "L50 right",
                        value = gp$l50_right,
                        min = 1,
                        max = signif(gp$l50_right * 2, 2),
                        step = 0.1),
            sliderInput("ldiff_right", "L25-L50 right",
                        value = gp$l25_right - gp$l50_right,
                        min = 0.1,
                        max = signif((gp$l25_right - gp$l50_right) * 2, 2),
                        step = 0.1)
        ))
    }
    l1
}
