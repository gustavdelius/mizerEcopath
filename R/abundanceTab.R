#' Abundance tab for tuning gadget
#'
#' The Abundance tab combines the [biomassTab()] and the [numberTab()]
#' on a single tab, providing a comprehensive view of both biomass and
#' individual numbers density across size bins and species.
#'
#' @inheritParams spectraTab
#' @family gadget tabs
#' @export
abundanceTab <- function(input,  ...) {
    biomassTab(input, ...)
    numberTab(input, ...)
}

#' @rdname abundanceTab
abundanceTabUI <- function(params, ...) {
    tagList(h1("Biomasses"),
            biomassTabUI(params, ...),
            h1("Numbers"),
            numberTabUI(params, ...))
}
