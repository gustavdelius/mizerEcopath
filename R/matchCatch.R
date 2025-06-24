#' Match observed catch, yield and production
#'
#' This function adjusts the gear‐selectivity, catchability and mortality for
#' the selected species so that a steady‑state model reproduces the observed
#' catch size distribution, the observed yield and the observed production, if
#' available.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' The function estimates optimal values for the following parameters:
#' * `catchability`: The catchability of the gear.
#' * `mu_mat`: The external mortality at maturity.
#' * `l25`: The size at which the gear selectivity first reaches 25%.
#' * `l50`: The size at which the gear selectivity first reaches 50%.
#'
#' If a *double‑sigmoid* selectivity function is used where selectivity drops
#' of again at large sizes, as in
#' `mizer::double_sigmoid_length()`, then also the parameters for the right
#' sigmoid are estimated:
#'
#' * `l50_right` – length at 50 % *descending* selectivity
#' * `l25_right` – length at 25 % *descending* selectivity
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
    # Accept count or catch as input
    if (!"count" %in% names(catch) && "catch" %in% names(catch)) {
        catch$count <- catch$catch
    }
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

    use_double_sigmoid <- identical(gps$sel_func, "double_sigmoid_length")

    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    # ── upper bound on mu_mat so juvenile slope < community spectrum ──
    g_mat     <- getEReproAndGrowth(params)[sp_select, mat_idx]
    mu_mat_max <- g_mat / w_mat * (lambda - sps$n)
    if (!"mu_mat" %in% names(sps) || is.na(sps$mu_mat)) {
        # determine external mortality at maturity
        mu_mat <- ext_mort(params)[sp_select, mat_idx]
    } else {
        mu_mat <- sps$mu_mat
    }

    ## ---- New parameterisation: optimise l50, ratio, Δ50, r_right, mu_mat, catchability ----

    # a) Build the initial parameter vector. Δ50 = l50_right - l50  (always > 0)
    initial_params <- list(
        l50          = gps$l50,
        ratio        = gps$l25 / gps$l50,
        r_right      = gps$l25_right / gps$l50_right,
        d50          = gps$l50_right - gps$l50,
        mu_mat       = mu_mat,
        # We need non-zero catchability
        catchability = pmax(gps$catchability, 1e-8)
    )

    # b) One simple, global bounds list
    default_bounds <- list(
        l50          = c(5,  Inf),
        ratio        = c(0.1, 0.8),
        d50          = c(5,  Inf),
        mu_mat       = c(0.2, mu_mat_max),
        catchability = c(1e-8, Inf),
        r_right      = c(1.3, 4)
    )

    # For single‐sigmoid we still *pass* these parameters to TMB,
    # but give them a small, finite value in the dome‐shaped code
    if (!use_double_sigmoid) {
        initial_params["d50"] <- default_bounds$d50[1]       # 10
        initial_params["r_right"] <- default_bounds$r_right[1] # 1.1
    }

    lower_bounds <- sapply(default_bounds, `[`, 1)
    upper_bounds <- sapply(default_bounds, `[`, 2)

    if (getOption("mizerEcopath.debug.matchCatch", FALSE)) {
        message("\nDEBUG: Optim bounds for ", species, ":\n")
        print(data.frame(lower = lower_bounds, upper = upper_bounds))
    }

    # Lock parameters where necessary
    map <- list()

    # Lock right-hand sigmoid parameters if using single sigmoid
    if (!use_double_sigmoid) {
        map$d50 <- factor(NA)
        map$r_right <- factor(NA)
        keep <- !names(lower_bounds) %in% c("d50", "r_right")
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
    }

    # Lock all selectivity parameters if no catch size data
    if (!data$use_counts) {
        map$l50 <- factor(NA)
        map$ratio <- factor(NA)
        map$d50 <- factor(NA)
        map$r_right <- factor(NA)
        keep <- !names(lower_bounds) %in% c("l50", "ratio", "d50", "r_right")
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
    }

    # Lock catchability if yield not to be matched
    if (data$yield_lambda == 0) {
        map$catchability <- factor(NA)
        keep <- names(lower_bounds) != "catchability"
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
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
