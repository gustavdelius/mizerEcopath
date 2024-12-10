#' Match the observed catch and yield
#'
#' This function adjusts various model parameters for the selected species so
#' that the model in steady state reproduces the observed catch size
#' distribution and the observed yield.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' The function sets new values for the following parameters:
#' * `l50`: The size at which the gear selectivity is 50%.
#' * `l25`: The size at which the gear selectivity is 25%.
#' * `catchability`: The catchability of the gear.
#' * `mu_mat`: The external mortality at maturity.
#' * `w_mat25`: The size at which 25% of individuals are mature.
#' It uses these to recalculate the corresponding rate arrays in the params
#' object. It sets the initial size spectrum to the steady state size spectrum.
#' The total biomass of each species remains unchanged. Only the selected
#' species are adjusted.
#'
#' The function estimates these parameters by minimizing an objective function.
#' The objective function is the negative log likelihood of the observed catch
#' size distribution given the probabilities predicted by the model plus the sum
#' of squares difference between the log of the observed yield and the log of
#' the predicted yield, multiplied by `yield_lambda`.
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
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @return A MizerParams object with the adjusted gear selectivity and
#'   catchability for the selected species
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
matchCatch <- function(params, species = NULL, catch, yield_lambda = 1) {
    species <- valid_species_arg(params, species = species,
                                 error_on_empty = TRUE)
    if (length(species) > 1) {
        for (s in species) {
            params <- matchCatch(params, species = s, catch = catch,
                                 yield_lambda = yield_lambda)
        }
        return(params)
    }

    data <- prepare_data(params, species = species, catch,
                         yield_lambda = yield_lambda)
    if (is.null(data)) { # There is no catch data for this species
        return(params)
    }

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
                        # we need non-zero catchability to match catch
                        catchability = max(gps$catchability, 1e-8))

    # Prepare the objective function.
    obj <- MakeADFun(data = data,
                     parameters = initial_params,
                     DLL = "mizerEcopath",
                     silent = TRUE)

    # Set parameter bounds
    lower_bounds <- c(l50 = 5, ratio = 0.1, mu_mat = 0, U = 1,
                      catchability = 1e-8)
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
