#' Align resource abundance to fish abundance
#'
#' This function sets the initial resource abundance density $N_R(w)$ to a power
#' law \deqn{N_R(w)=\kappa w^{-\lambda}}, where $lambda$ is taken from
#' `resource_params(params)$lambda` The coefficient `kappa` is chosen
#' such that the resource abundance power law is tangent to the abundance
#' density \eqn{N_c(w)} of the fish community in one point. It then truncates the
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
    N_R <- params@w_full^(-params@resource_params$lambda)
    # Set resource to be in line with fish
    total <- colSums(initialN(params))
    fish_sel <- params@w_full >= params@w[1]
    # Find the amount by which the log of the resource abundance needs to be
    # decreased to be tangent to the log of the fish abundance.
    # log(N_R) <- log(N_R) - dist_log
    dist_log <- min(log(N_R[fish_sel]) - log(total))
    # So N_R <- N_R * exp(-dist_log)
    ratio <- exp(-dist_log)
    N_R <- N_R * ratio
    # Truncate resource
    N_R[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    # Set initial resource abundance
    initialNResource(params) <- N_R
    # Set resource carrying capacity
    resource_capacity(params) <- N_R
    # Update resource parameter
    # The initial coefficient of the N_R power law was 1, so after rescaling
    # by `ration` the new coefficient is just `ratio`
    params@resource_params$kappa <- ratio
    return(params)
}
