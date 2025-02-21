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
#' If `min_w_pred` is set, then the diet matrix only represents the flow of
#' biomass into the part of the predator species population that is larger than
#' `min_w_pred`. This may make the returned diet matrix correspond more closely
#' to diet matrices found in the literature, which are often estimated from
#' observation of large individuals only and ignore that the diet of small
#' individuals will be different.
#'
#' @param params A MizerParams object
#' @param min_w_pred The minimum weight of predators to include in the diet
#'   matrix. Default is 0.
#' @return A matrix with dimnames `predator` and `prey` containing the rates of
#'   biomass flow from a prey species to a predator species (or that part of
#'   the predator species population that is larger than `min_w_pred`).
#' @export
#' @family rate functions
#' @examples
#' dimnames(getDietMatrix(NS_params))
#' getDietMatrix(NS_params)["Cod", ]
getDietMatrix <- function(params, min_w_pred = 0) {
    w_sel <- params@w >= min_w_pred
    N <- initialN(params)[, w_sel]
    dw <- dw(params)[w_sel]
    diet_matrix <- getDiet(params, proportion = FALSE)[, w_sel, ] |>
        sweep(c(1, 2), N, "*") |>
        sweep(2, dw, "*") |>
        aperm(c(1, 3, 2)) |>
        rowSums(dims = 2)
    return(diet_matrix)
}
