library(mizerEcopath)
library(here)

source(here("inst/prepareLogLikelihoodFn.R"))
source(here("inst/objective_function.R"))

# Load catch size distribution
catch <- readRDS(here("../CelticSea/data/ecopath_catch_v_FMP.rds"))
catch$dl <- 1
catch$gear <- as.character(catch$gear)
catch$gear[catch$gear == "commercial"] <- "total"
catch <- catch[catch$species == "Cod" & catch$gear == "total", ]
catch$count <- catch$catch

# Extract a single-species model from an existing model
p <- readRDS(here("../CelticSea/ecopath/Lauria/lauria_HOM_2.rds"))
species <- species_params(p)$species
remove <- species[species != "Cod"]
p <- removeSpecies(p, remove)
sp <- species_params(p)
sp$M <- mean(getExtMort(p) / w(p)^sp$d)
gp <- gear_params(p)
model_lengths <- (w(p) / sp$a)^(1/sp$b)

# Prepare the log-likelihood function
log_likelihood_fn <- prepare_log_likelihood(catch)

# Initial parameter estimates
initial_params <- c(l50 = gp$l50, ratio = gp$l25 / gp$l50, M = 0)

# Set parameter bounds (if necessary)
lower_bounds <- c(l50 = 5, ratio = 0.1, M = 0)
upper_bounds <- c(l50 = Inf, ratio = 0.99, M = Inf)

# Perform the optimization
optim_result <- optim(
    par = initial_params,
    fn = objective_function,
    method = "L-BFGS-B",  # Allows for parameter bounds
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(fnscale = 1, maxit = 1000)
)

# Extract the optimal parameters
optimal_params <- optim_result$par
print(optimal_params)

# Optimal negative log-likelihood value
optimal_neg_log_likelihood <- optim_result$value
print(optimal_neg_log_likelihood)

# Plot observed counts
observed_data <- catch
hist_data <- data.frame(
    length = observed_data$length + observed_data$dl / 2,
    count = observed_data$count
)

plot(hist_data$length, hist_data$count, type = 'h', lwd = 2, col = 'blue',
     xlab = 'Length', ylab = 'Count', main = 'Observed Data and Fitted PDF')

# Generate fitted PDF with optimal parameters
pdf_result <- catch_pdf(optimal_params)
pdf_lengths <- as.numeric(names(pdf_result))
pdf_values <- as.numeric(pdf_result)

# Scale the PDF for visualization
pdf_values_scaled <- pdf_values * max(hist_data$count) / max(pdf_values)

# Add the fitted PDF to the plot
lines(pdf_lengths, pdf_values_scaled, col = 'red', lwd = 2)
legend('topright', legend = c('Observed Counts', 'Fitted PDF'),
       col = c('blue', 'red'), lwd = 2)

