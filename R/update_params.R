#' Prepare a TMB Objective Function for Optimising Model Parameters
#'
#' This function returns a list with the data to be passed to the TMB objective
#' function. The data includes the observed catch data, the model parameters,
#' and some precomputed values that are used in the likelihood calculation.
#' The main preprocessing makes sure that we have a comprehensive set of bins
#' that cover the entire size range, even though there will not be observations
#' at all sizes. Missing observations should be interpreted as a 0 count.
#'
#' @param params A MizerParams object
#' @param species The species for which the data is to be prepared. By default
#'   the first species in the model.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#'
#' @return A list with the data to be passed to the TMB objective function. If
#'   there is no catch data for the species, the function returns NULL.
#' @export
update_params <- function(params, species = 1, pars, data) {
    params <- validParams(params)
    sp <- species_params(params)

    # Need to add mu_mat column if it does not exist
    # TODO: change once we have introduced a standard species param for this
    mat_idx <- colSums(outer(params@w, sp$w_mat, "<"))
    mu_mat <- ext_mort(params)[cbind(seq_len(nrow(sp)), mat_idx)]
    params <- set_species_param_default(params, "mu_mat", mu_mat)

    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Only one species can be updated at a time.")
    }

    sp_select <- sp$species == species
    sps <- sp[sp_select, ]

    gp <- params@gear_params
    gp_select <- gp$species == species
    gps <- gp[gp_select, ]
    # if (nrow(gps) > 1) {
    #   stop("The code currently assumes that there is only a single gear for each species.")
    # }   # not true now

    # Update the gear parameters
    gplist <- list()
    gpnames <- c( 'logit_l50', 'log_ratio_left', 'log_l50_right_offset', 'log_ratio_right',
                  'log_catchability')

    for (i in gpnames) gplist[[i]] <- as.numeric(pars[[i]])

    l50 <- min(data$l) + (max(data$l) - min(data$l)) * plogis(gplist$logit_l50)
    l25 <- l50 * (1 - plogis(gplist$log_ratio_left))
    l50_right <- l50 + exp(gplist$log_l50_right_offset)
    l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
    catchability <- exp(gplist$log_catchability)

    gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)

    gps[,'l50'] <- gp_res$l50
    gps[,'l25'] <- gp_res$l25
    gps[,'l50_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l50_right, NA)
    gps[,'l25_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l25_right, NA)
    gps[,'catchability'] <- gp_res$catchability

    gear_params(params)[gp_select, ] <- gps

    # recalculate the power-law mortality rate
    sps$mu_mat <- pars[["mu_mat"]]
    sps$m <- pars[["m"]]
    # Note that `mu_mat` is the mortality at the w just below w_mat
    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    ext_mort(params)[sp_select, ] <-
        pars[["mu_mat"]] * (params@w / w_mat)^sps$d

    params@species_params[sp_select, ] <- sps
    params <- setReproduction(params)

    # Calculate the new steady state

    params <- mizer::steadySingleSpecies(params, species = species)
    # Rescale it to get the observed biomass
    params <- matchBiomasses(params, species = species)
    # Set the reproduction level to zero
    rl <- numeric(length(species))
    names(rl) <- species
    params <- setBevertonHolt(params, reproduction_level = rl)

    return(params)
}
