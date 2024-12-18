#' Get diet matrix
#'
#' Returns the diet matrix with one row for each predator species and one
#' column for each prey species and other ecosystem components. The entries
#' give the rate at which biomass flows from the prey species to the
#' predator species.
#'
#' Note that the diet matrix has the absolute rates of biomass flow, not the
#' proportions of the diet.
#'
#' @param params A MizerParams object
#' @param min_w_prey The minimum weight of prey species to include in the diet
#'   matrix
#' @param max_w_prey The maximum weight of prey species to include in the diet
#'   matrix
#' @return A matrix with dimnames `predator` and `prey` containing the rates of
#'   biomass flow from a prey species to a predator species.
#' @export
#' @family rate functions
#' @examples
#' dimnames(getDietMatrix(NS_params))
#' getDietMatrix(NS_params)["Cod", ]
getDietMatrix <- function(params, min_w_prey = 0, max_w_prey = Inf) {
    w_sel <- params@w >= min_w_prey & params@w <= max_w_prey
    N <- initialN(params)[, w_sel]
    dw <- dw(params)[w_sel]
    diet_matrix <- getDiet(params, proportion = FALSE)[, w_sel, ] |>
        sweep(c(1, 2), N, "*") |>
        sweep(2, dw, "*") |>
        aperm(c(1, 3, 2)) |>
        rowSums(dims = 2)
    return(diet_matrix)
}
