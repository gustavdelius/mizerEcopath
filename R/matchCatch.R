#' Match the observed catch and yield
#'
#' This function adjusts various model parameters for the selected species so
#' that the model in steady state reproduces the observed catch size
#' distribution, the observed yield and the observed production, if available.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' The function sets new values for the following parameters:
#' * `l50`: The size at which the gear selectivity is 50%.
#' * `l25`: The size at which the gear selectivity is 25%.
#' * `catchability`: The catchability of the gear.
#' * `mu_mat`: The external mortality at maturity.
#'
#' Only the parameters of the selected species are adjusted. The function then
#' recalculates the corresponding rate arrays in the params object. It sets the
#' initial size spectrum to the steady state size spectrum. The total biomass of
#' each species remains unchanged.
#'
#' The function estimates these parameters by minimizing an objective function.
#' The objective function is the negative log likelihood of the observed catch
#' size distribution given the probabilities predicted by the model plus the sum
#' of squares difference between the log of the observed yield and the log of
#' the predicted yield, multiplied by `yield_lambda`, as well as the sum of
#' squares difference between the log of the observed production and the log of
#' the predicted production, multiplied by `production_lambda`.
#'
#' The function deals with missing data in the following way, for each species
#' individually:
#'
#' -  If the observed yield is not available, the function will only match the
#' observed catch size distribution and the observed production.
#'
#' -  If the observed production is not available, the function will only match the
#' observed catch size distribution and the observed yield.
#'
#' -  If the observed catch size distribution is not available, the function will
#' only match the observed yield and the observed production.
#'
#' -  If neither the observed yield nor the observed production are available, the
#' function raises an error.
#'
#' The catch predicted by the model is calculated by integrating the
#' catchability and the gear selectivity over the size distribution of the
#' species. The size distribution itself is shaped by the model through the
#' interplay between growth and mortality. This is why the gear selectivity and
#' catchability need to be adjusted together with the other parameters that
#' shape the size distribution.
#'
#' The objective function is coded in C++ and the TMB package is used to compile
#' the objective function and to create functions for automatically calculating
#' the gradients. These are then passed to the `nlminb` function to minimize the
#' objective function.
#'
#' @param params A MizerParams object
#' @param species The species for which to match the catch. Optional. By default
#'   all target species are selected. A vector of species names, or a numeric
#'   vector with the species indices, or a logical vector indicating for each
#'   species whether it is to be selected (TRUE) or not.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param lambda The slope of the community spectrum. Default is 2.05.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @param production_lambda A parameter that controls the strength of the penalty
#'  for deviation from the observed production.
#'
#' @return A MizerParams object with the adjusted external mortality, gear
#'   selectivity, catchability and steady-state spectrum for the selected
#'   species.
#' @family match functions
#' @examples
#' params <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
#' plot_catch(params, species = "Hake", catch = celtic_catch)
#' # The function leaves the biomass of the species unchanged
#' all.equal(getBiomass(params), getBiomass(celtic_params), tol = 1e-4)
#' # It also leaves the energy available to an individual for reproduction
#' # and growth unchanged
#' all.equal(getEReproAndGrowth(params), getEReproAndGrowth(celtic_params))
#' # The initial size spectrum is set to the steady state size spectrum
#' params_steady <- steadySingleSpecies(params)
#' all.equal(initialN(params), initialN(params_steady))
#' @export
matchCatch <- function(params, species = NULL, catch, lambda = 2.05,
                       yield_lambda = 1, production_lambda = 1) {
    species <- valid_species_arg(params, species = species,
                                 error_on_empty = TRUE)
    params <- validParams(params)
    if (length(species) > 1) {
        for (s in species) {
            params <- matchCatch(params, species = s, catch = catch,
                                 yield_lambda = yield_lambda,
                                 production_lambda = production_lambda)
        }
        return(params)
    }

    data <- prepare_data(params, species = species, catch,
                         yield_lambda = yield_lambda,
                         production_lambda = production_lambda)
    if (is.null(data)) {
        warning(species, " can not be matched because neither catches nor production are given.")
        return(params)
    }

    sp <- species_params(params)
    gp <- gear_params(params)
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]
    gps <- gp[gp$species == species, ]

    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    if (!"mu_mat" %in% names(sps) || is.na(sps$mu_mat)) {
        # determine external mortality at maturity
        mu_mat <- ext_mort(params)[sp_select, mat_idx]
    } else {
        mu_mat <- sps$mu_mat
    }

    # Initial parameter estimates
    initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50,
                        mu_mat = mu_mat,
                        # we need non-zero catchability to match catch
                        catchability = max(gps$catchability, 1e-8))

    # Set parameter bounds
    # Mortality is bounded by the requirement that the juvenile spectrum of
    # each species must be less steep than the community spectrum.
    # With g(w) = g w^n and mu(w) = mu w^{n-1}, the exponent of the juvenile
    # spectrum is -mu/g-n. The exponent of the community spectrum is -lambda.
    g_mat <- getEReproAndGrowth(params)[sp_select, mat_idx]
    mu_mat_max <- g_mat / w_mat * (lambda - sps$n)
    lower_bounds <- c(l50 = 5, ratio = 0.1, mu_mat = 0,
                      catchability = 1e-8)
    upper_bounds <- c(l50 = Inf, ratio = 0.99, mu_mat = mu_mat_max,
                      catchability = Inf)

    map <- list()
    if (!data$use_counts) {
        map$l50 <- factor(NA)
        map$ratio <- factor(NA)
        lower_bounds <- lower_bounds[-(1:2)]
        upper_bounds <- upper_bounds[-(1:2)]
    }
    if (data$yield_lambda == 0) {
        map$catchability <- factor(NA)
        lower_bounds <- lower_bounds[-4]
        upper_bounds <- upper_bounds[-4]
    }
    # Prepare the objective function.
    obj <- MakeADFun(data = data,
                     parameters = initial_params,
                     map = map,
                     DLL = "mizerEcopath",
                     silent = TRUE)

    # Perform the optimization.
    optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                           lower = lower_bounds, upper = upper_bounds,
                           control = list(trace = 0))

    # Set model to use the optimal parameters
    w_select <- w(params) %in% data$w
    optimal_params <- update_params(params, species, optim_result$par,
                                    data, w_select)

    return(optimal_params)
}
