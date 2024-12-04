# In this script we calibrate the cod model to the observed catch size distribution
# and total yield. We use the L-BFGS-B optimization algorithm to find the parameter
# values that minimize the negative log-likelihood of the observed catch data.
# The objective function also includes a penalty for deviation from the observed
# total yield.

library(mizerEcopath)
library(TMB)
p <- celtic_params
sp <- species_params(p)
gp <- gear_params(p)
catch <- celtic_catch

# Set parameter bounds
lower_bounds <- c(l50 = 5, ratio = 0.1, M = 0, U = 1, catchability = 0)
upper_bounds <- c(l50 = Inf, ratio = 0.99, M = Inf, U = 20, catchability = Inf)

species <- valid_species_arg(p, 8)

sp_select <- sp$species == species
sps <- sp[sp_select, ]
gps <- gp[gp$species == species, ]

# Model does not fit the observed catch yet:
plot_catch(p, species, catch)
# and it has the wrong total yield:
gps$yield_observed
getYield(p)[sp_select]

# Initial parameter estimates
initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50, M = 2, U = 10,
                    catchability = gps$catchability)

# Prepare the objective function.
data <- prepare_data(p, species, catch, yield_lambda = 1e7)
obj <- MakeADFun(data = data,
                 parameters = initial_params,
                 DLL = "mizerEcopath")

# Perform the optimization. This starts with the initial parameter estimates and
# iteratively updates them to minimize the objective function.
optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                       lower = lower_bounds, upper = upper_bounds,
                       control = list(trace = 0))

optim_result$par
report <- obj$report()

# Set model to use the optimal parameters
optimal_params <- update_params(p, species, optim_result$par)
# and plot the model catch again against the observed catch
plot_catch(optimal_params, species, catch)

# Also the yield is approximately matched:
gps$yield_observed
report$model_yield
getYield(optimal_params)[sp_select]
# If you want a better match you can increase the `yield_lambda` parameter

# Biomass is matched perfectly, by design
sps$biomass_observed
getBiomass(optimal_params)[sp_select]
