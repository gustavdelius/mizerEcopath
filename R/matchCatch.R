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
#' * `l50_right`: The size at which the gear selectivity is 50% (right side for double_sigmoid_length selectivity).
#' * `l25_right`: The size at which the gear selectivity is 25% (right side for double_sigmoid_length selectivity).
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
#'   * `species`: The species for which the catch is recorded.
#'   * `gear`: The gear used to collect the catch.
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `catch`: The observed count for each bin.
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
#' all.equal(getBiomass(params, use_cutoff = TRUE),
#'           getBiomass(celtic_params, use_cutoff = TRUE),
#'           tol = 1e-4)
#' # It also leaves the energy available to an individual for reproduction
#' # and growth unchanged
#' all.equal(getEReproAndGrowth(params),
#'           getEReproAndGrowth(celtic_params))
#' # The initial size spectrum is set to the steady state size spectrum
#' params_steady <- steadySingleSpecies(params)
#' all.equal(initialN(params), initialN(params_steady))
#' @export
#'
matchCatch <- function(params, species = NULL, catch, lambda = 2.05,
                       yield_lambda = 1, production_lambda = 1, mu_mat_lim = 5, map = NULL)
{
    species <- valid_species_arg(params, species = species, error_on_empty = TRUE)
    params <- validParams(params)

    if (length(species) > 1) {
        for (s in species) {
            params <- matchCatch(params, species = s, catch = catch, lambda = lambda,
                                 yield_lambda = yield_lambda, production_lambda = production_lambda,
                                 mu_mat_lim = mu_mat_lim, map = map)
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

    if(is.vector(data$counts)){
        data$counts <- as.matrix(data$counts, ncol=1)
        colnames(data$counts) <- gps$gear}

    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    if (!"mu_mat" %in% names(sps) || is.na(sps$mu_mat)) {
        # determine external mortality at maturity
        mu_mat <- ext_mort(params)[sp_select, mat_idx]
    } else {
        mu_mat <- sps$mu_mat
    }

    # Initial parameter estimates
    # differenet parametrization for selectivity (note that C++ fails with NAs)
    initial_params <- list(

        logit_l50 = qlogis((gps$l50 - min(data$l))/(max(data$l) - min(data$l))),

        log_ratio_left = qlogis((gps$l50 - gps$l25)/gps$l50),

        log_l50_right_offset = ifelse( gps$sel_func == 'double_sigmoid_length',
                                       log(pmax(1e-3, gps$l50_right - gps$l50)), 1),

        log_ratio_right = ifelse( gps$sel_func == 'double_sigmoid_length',
                                  log((gps$l25_right - gps$l50_right)/gps$l50_right), 1),

        log_catchability = log(ifelse(gps$catchability <= 0, 1e-8, gps$catchability)),

        mu_mat = mu_mat,

        m = ifelse(any(names(sps)=='m'), sps$m, 1))

    # Set parameter bounds
    # Mortality is bounded by the requirement that the juvenile spectrum of
    # each species must be less steep than the community spectrum.
    # With g(w) = g w^n and mu(w) = mu w^{n-1}, the exponent of the juvenile
    # spectrum is -mu/g-n. The exponent of the community spectrum is -lambda.
    g_mat <- getEReproAndGrowth(params)[sp_select, mat_idx]
    mu_mat_max <- g_mat / w_mat * (lambda - sps$n)

    lower_bounds_list <- list(
        logit_l50 = rep(-10, length(data$sel_func)),
        log_ratio_left = rep(-10, length(data$sel_func)),
        log_l50_right_offset = rep(-10, length(data$sel_func)),
        log_ratio_right = rep(-10, length(data$sel_func)),
        log_catchability = rep(-10, length(data$sel_func)),
        mu_mat = 0.2,
        m = sps$n * 1.01
    )

    upper_bounds_list <- list(
        logit_l50 = rep(10, length(data$sel_func)),
        log_ratio_left = rep(10, length(data$sel_func)),
        log_l50_right_offset = rep(10, length(data$sel_func)),
        log_ratio_right = rep(10, length(data$sel_func)),
        log_catchability = rep(10, length(data$sel_func)),
        mu_mat = min(mu_mat_max,mu_mat_lim),
        m = 3
    )

    if (!is.null(map)) {
        # Traducir los nombres amigables a los nombres internos del optimizador
        name_translation <- c(
            "l50" = "logit_l50",
            "l25" = "log_ratio_left",
            "l50_right" = "log_l50_right_offset",
            "l25_right" = "log_ratio_right",
            "catchability" = "log_catchability",
            "mu_mat" = "mu_mat",
            "m" = "m"
        )

        transformed_map <- list()
        for (orig_name in names(map)) {
            if (orig_name %in% names(name_translation)) {
                transformed_map[[name_translation[[orig_name]]]] <- map[[orig_name]]
            } else {
                transformed_map[[orig_name]] <- map[[orig_name]]
            }
        }
        map <- transformed_map
    } else {
        map <- list()
    }

    if (!data$use_counts) {
        map$logit_l50 <- factor(rep(NA, length(data$sel_func)))
        map$log_ratio_left <- factor(rep(NA, length(data$sel_func)))
        map$log_l50_right_offset <- factor(rep(NA, length(data$sel_func)))
        map$log_ratio_right <- factor(rep(NA, length(data$sel_func)))
    }
    if (data$yield_lambda == 0) {
        map$log_catchability <- factor(rep(NA, length(data$sel_func)))
    }

    lower_bounds <- NULL
    upper_bounds <- NULL

    for (p in names(initial_params)) {
        lb <- lower_bounds_list[[p]]
        ub <- upper_bounds_list[[p]]

        if (!is.null(map[[p]])) {
            keep <- !is.na(as.vector(map[[p]]))
            lb <- lb[keep]
            ub <- ub[keep]
        }
        lower_bounds <- c(lower_bounds, lb)
        upper_bounds <- c(upper_bounds, ub)
    }

    # Prepare the objective function.

    obj <- TMB::MakeADFun(data = data, parameters = initial_params, map = map,
                          DLL = "mizerEcopath", silent = TRUE, debug = FALSE)

    # Perform the optimization.
    optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                           lower = lower_bounds, upper = upper_bounds,
                           control = list(trace = 0))

    # Set model to use the optimal parameters
    w_select <- w(params) %in% data$w

    # Get all parameters, including fixed ones
    pars_all <- obj$env$last.par.best
    if (is.null(pars_all)) pars_all <- obj$env$last.par

    pars_list <- obj$env$parList(pars_all)

    optimal_params <- update_params(params, species, pars_list, data)

    return(optimal_params)
}
