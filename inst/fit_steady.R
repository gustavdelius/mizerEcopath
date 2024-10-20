# In this script we will use `optim()` to fit the parameters of a mizer model to
# given ecopath parameters and given size-distribution data.

library(mizerEcopath)
library(here)

# Extract a single-species model from an existing model
p <- readRDS(here("../CelticSea/ecopath/Lauria/lauria_HOM_2.rds"))
species <- species_params(p)$species
remove <- species[species != "Cod"]
p <- removeSpecies(p, remove)
sp <- species_params(p)
gp <- gear_params(p)
isAllometric(p, tol = 2e-5)

w_mat_idx <- sum(w(p) < sp$w_mat)
E <- getEReproAndGrowth(p)[1, w_mat_idx]


# Load catch size distribution
catch <- readRDS(here("../CelticSea/data/ecopath_catch_v_FMP.rds"))
catch$dl <- 1
catch$gear <- as.character(catch$gear)
catch$gear[catch$gear == "commercial"] <- "total"
catch <- catch[catch$species == "Cod" & catch$gear == "total", ]
catch$count <- catch$catch
lengths <- (min(catch$length) - 5):(max(catch$length) * 1.1)
weights <- sp$a * lengths ^ sp$b

# We keep `E` and `ks` fixed because we can change them later to the values
# that give the desired values for `Q` and `age_mat`. This later change will
# not affect the steady state of the model because the change in growth rate
# that it leads to can be compensated by changing all other rates by the same
# factor. We just have to make sure that the parameters we tune and the data
# we tune to are dimensionless or have dimension involving length or weight only
#
# We will attempt to estimate the parameters `mu/E`, `catchability/E`, `L50`,
# `L25`, `w_repro_max`, and `m` for each species in the model.

# In this first iteration we will take `w_mat` as known and fixed. In future
# we might estimate it as well while penalising deviations from the specified
# value.

# The objective function will determine the steady state with given parameters
# and calculate the distance from the following data:
# P/Q, C/Q, Catch-size-distribution
# It will also penalise deviations from a reproductive efficiency of 0.01.
objective_fn <- function(pars) {
    sp$mu <- pars[1]
    gp$catchability <- pars[2]
    gp$L50 <- pars[3]
    gp$L25 <- pars[4] * L50
    sp$w_repro_max <- pars[6]
    # sp$m <- pars[7]
    species_params(p) <- sp
    gear_params(p) <- gp
    p <- steadySingleSpecies(p)

    catch_dens <- initialN(p) * getFMort(p) * dw(p)
    catch_density <- stats::approx(x = w(p), y = catch_dens,
                                   xout = weights)$y
    catch_prob <- catch_density / sum(catch_density)
    names(catch_prob) <- lengths
    LL <- log_likelihood(catch_prob, catch)
}

log_likelihood <- function(p, df) {
    # Ensure the probability vector 'p' has names corresponding to the integer lengths
    if (is.null(names(p))) {
        stop("The probability vector 'p' must have names corresponding to the lengths.")
    }

    # Get all lengths from the probability distribution
    all_lengths <- names(p)

    # Initialize counts to zero for all lengths
    counts_all <- setNames(rep(0, length(all_lengths)), all_lengths)

    # Update counts with observed counts from the data frame
    observed_lengths <- as.character(df$length)
    counts_all[observed_lengths] <- df$count

    # Total number of observations
    N <- sum(counts_all)

    # Ensure all probabilities are positive
    if (any(p <= 0)) {
        stop("All probabilities must be positive.")
    }

    # Normalize probabilities in case they do not sum to 1
    p_normalized <- p / sum(p)

    # Compute the log-likelihood using the multinomial distribution formula
    LL <- lgamma(N + 1) - sum(lgamma(counts_all + 1)) +
        sum(counts_all * log(p_normalized))

    return(LL)
}

counts <- setNames(rep(0, length(lengths)), lengths)
# Update counts with observed counts from the data frame
observed_lengths <- as.character(catch$length)
counts[observed_lengths] <- catch$catch

plot(lengths, catch_prob, type = "l")
lines(lengths, counts / sum(counts), col = "red")




