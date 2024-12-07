#' Match the observed yield by adjusting catchability
#'
#' This function adjusts the gear selectivity parameters and the catchability
#' for the selected species so that the model in steady state reproduces the
#' observed catch size distribution and the observed total observed yield.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' @param params A MizerParams object
#' @param species The name of the selected species. By default the first
#'   species in the model.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @return A MizerParams object with the adjusted gear selectivity and
#'   catchability for the selected species
#' @family match functions
#' @export
matchCatch <- function(params, species = 1, catch, yield_lambda = 1) {
    data <- prepare_data(params, species = species, catch,
                         yield_lambda = yield_lambda)
    species <- valid_species_arg(params, species = species)
    sp <- species_params(params)
    gp <- gear_params(params)
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]
    gps <- gp[gp$species == species, ]

    if (!"mu_mat" %in% names(sps) || is.na(sps$mu_mat)) {
        # determine external mortality at maturity
        mat_idx <- sum(params@w < sps$w_mat)
        mu_mat <- ext_mort(params)[sp_select, mat_idx]
    } else {
        mu_mat <- sps$mu_mat
    }

    # Steepness of maturity ogive
    U <- log(3) / log(sps$w_mat / sps$w_mat25)

    # Initial parameter estimates
    initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50,
                        mu_mat = mu_mat, U = U,
                        catchability = gps$catchability)

    # Prepare the objective function.
    obj <- MakeADFun(data = data,
                     parameters = initial_params,
                     DLL = "mizerEcopath",
                     silent = TRUE)

    # Set parameter bounds
    lower_bounds <- c(l50 = 5, ratio = 0.1, mu_mat = 0, U = 1,
                      catchability = 0)
    upper_bounds <- c(l50 = Inf, ratio = 0.99, mu_mat = Inf, U = 20,
                      catchability = Inf)

    # Perform the optimization.
    optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                           lower = lower_bounds, upper = upper_bounds,
                           control = list(trace = 0))

    # Set model to use the optimal parameters
    w_select <- w(params) %in% data$w
    optimal_params <- update_params(params, species, optim_result$par,
                                    data$biomass, w_select)

    return(optimal_params)
}

#' @rdname matchCatch
isCatchMatched <- function(params, tol = 0.1) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    sp <- params@species_params
    gp <- params@gear_params
    if (!hasName(gp, "yield_observed")) {
        stop("You must provide the yield_observed gear parameter.")
    }
    # Calculate discrepancy in yields
    Cratio <- gp$yield_observed / getYield(params)

    return(max(abs(Cratio - 1)) < tol)
}
