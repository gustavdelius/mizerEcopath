# In this script we calibrate the cod model to the observed catch size distribution
# and total yield. We use the L-BFGS-B optimization algorithm to find the parameter
# values that minimize the negative log-likelihood of the observed catch data.
# The objective function also includes a penalty for deviation from the observed
# total yield.

library(mizerEcopath)
library(TMB)
params <- celtic_params
sp <- species_params(params)
gp <- gear_params(params)
catch <- celtic_catch

params <- matchGrowth(params, keep = "biomass")
params <- matchCatch(params, catch = catch)

species <- valid_species_arg(params, 5)

# Model does not fit the observed catch yet:
plot_catch(params, species, catch)


data <- prepare_data(params, species = species, catch,
                     production_lambda = 1)
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

# Initial parameter estimates
initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50,
                    mu_mat = mu_mat,
                    # we need non-zero catchability to match catch
                    catchability = max(gps$catchability, 1e-8))

# Prepare the objective function.
obj <- MakeADFun(data = data,
                 parameters = initial_params,
                 DLL = "mizerEcopath",
                 silent = TRUE)

# Set parameter bounds
lower_bounds <- c(l50 = 5, ratio = 0.1, mu_mat = 0,
                  catchability = 1e-8)
upper_bounds <- c(l50 = Inf, ratio = 0.99, mu_mat = Inf,
                  catchability = Inf)

# Perform the optimization.
optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                       lower = lower_bounds, upper = upper_bounds,
                       control = list(trace = 0))

# Set model to use the optimal parameters
w_select <- w(params) %in% data$w
optimal_params <- update_params(params, species, optim_result$par,
                                data$biomass, w_select)


plot_catch(optimal_params, species, catch)

report <- obj$report()

# Also the yield is approximately matched:
gps$yield_observed
report$model_yield
getYield(optimal_params)[sp_select]
# If you want a better match you can increase the `yield_lambda` parameter

sps$production_observed
report$model_production
getSomaticProduction(optimal_params)[sp_select]

# Check that TMB code and mizer agree on size distribution
plot(data$w, report$N, type = "l", log = "y")
lines(data$w, initialN(optimal_params)[sp_select, w_select], col = "red")

# Biomass is matched perfectly, by design
data$biomass
report$total_biomass
getBiomass(optimal_params, min_w = data$w[1], max_w = max(data$w))[sp_select]
sum(initialN(optimal_params)[sp_select, w_select] * data$w * data$dw)


