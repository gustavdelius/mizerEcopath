#' Align resource abundance to fish abundance
#'
#' This function sets the initial resource abundance density $N_R(w)$ to a power
#' law \deqn{N_R(w)=\kappa w^{-\lambda}}, where $lambda$ is taken from
#' `resource_params(params)$lambda` The coefficient $\kappa=$`kappa` is chosen
#' such that the resource abundance power law is tangent to the abundance
#' density $N_c(w)$ of the fish community in one point. It then truncates the
#' resource abundance power law at `resource_params(params)$w_pp_cutoff`. The
#' function also sets the resource carrying capacity equal to this initial
#' resource abundance and updates the resource parameter `kappa` in the
#' MizerParams object.
#'
#' @param params A MizerParams object
#' @return A MizerParams object with the resource abundance scaled up or down
#'   to match the fish abundance
#' @export
alignResource <- function(params) {
    # Set resource to a power law
    N_resource <- params@w_full^(-params@resource_params$lambda)
    # Set resource to be in line with fish
    total <- colSums(initialN(params))
    fish_sel <- params@w_full >= params@w[1]
    ratio <- max(total / initialNResource(params)[fish_sel])
    N_resource <- N_resource * ratio
    # Truncate resource
    N_resource[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    # Set initial resource abundance
    initialNResource(params) <- N_resource
    # Set resource carrying capacity
    resource_capacity(params) <- N_resource
    # Update resource parameter
    kappa <- N_resource[1] * params@w_full[1]^params@resource_params$lambda
    params@resource_params$kappa <- kappa
    return(params)
}
