
library(mizerEcopath)
library(Matrix)      # sparse matrices, fast %*%
library(dplyr)
library(ggplot2)

# Load age at size data compiled by Jess ----

survey <- readRDS(here::here("inst", "extdata",
                             "Celtic_Sea_Size_at_Age_Data.rds"))
# Aggregate to cm bins
Y_obs <- survey |>
    filter(Scientific_name == "Gadus morhua") |>
    # remove rows with NA in any column
    filter(!is.na(LngtClass) & !is.na(Age) & !is.na(CANoAtLngt)) |>
    # round down to cm
    mutate(LngtClass = floor(LngtClass)) |>
    # aggregate counts
    group_by(Quarter, LngtClass, Age) |>
    summarise(CANoAtLngt = sum(CANoAtLngt, na.rm = TRUE), .groups = "drop") |>
    transmute(qtr = Quarter,
              lenBin = as.integer(LngtClass),
              ageClass = as.integer(Age),
              count = CANoAtLngt)

# survey length‑bin edges [cm]
S_edge <- min(Y_obs$lenBin):(max(Y_obs$lenBin) + 1)
s <- length(S_edge) - 1 # number of survey bins

# highest survey age
K <- max(Y_obs$ageClass)

# list of fractional survey dates
f_h <- c(0.125, 0.625, 0.875)
quarter_h <- c(1, 3, 4) # quarter numbers
H <- length(f_h)


# Evolve a cohort over time in model ----

## Set up example ----
# We will use the Celtic Sea Cod for this example
params <- celtic_params
species <- "Cod"
sps <- species_params(params)[species, ]

# model length‑bin edges [cm]
w <- w(params)
l <- (w/sps$a)^(1/sps$b)
L_edge <- l
m <- length(L_edge) - 1  # number of model bins

# Set initial condition
n_init <- rep(0, length(l))
n_init[2] <- 1  # Point source at small size
initialN(params)[species, ] <- n_init

## Solve the PDE ----

# fine‑age grid (centres) and widths  [years]
dt <- 0.05
t_max <- 6
age_mid   <- seq(dt, t_max,  by = dt)
nsteps <- length(age_mid)
age_width <- rep(dt, nsteps)

n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

## Plot the cohort densities at selected time points ----

# Select time points to plot
plot_steps <- c(21, 41, 61, 81, 101, 121)
plot_times <- (plot_steps - 1) * dt
plot_labels <- paste0("age=", round(plot_times, 2))

# Prepare data frame for ggplot2
N <- length(n_hist[1, ])
plot_df <- data.frame(
    Weight = rep(w, times = length(plot_steps)),
    Length = rep(l, times = length(plot_steps)),
    Abundance = as.vector(t(n_hist[plot_steps, ]) + 1e-10),
    t = rep(plot_times, each = N),
    Time = factor(rep(plot_labels, each = N), levels = plot_labels)
)

# Plot
p <- ggplot(plot_df, aes(x = Length, y = Abundance*exp(t*0.5), color = Time)) +
    geom_line(linewidth = 1) +
    #scale_x_log10() +
    #scale_y_log10() +
    xlim(1e-3, 140) +
    labs(x = "Length (cm)", y = "Abundance density") +
    theme_minimal()
print(p)

## Convert to numbers from densities ----

N_mat <- n_hist[2:(nsteps + 1), 1:m]
N_mat <- sweep(N_mat, 2, dw(params)[1:m], "*")

# Build aggretation matrices ----

## AGE aggregation A_list[[h]] :  (K+1) × length(age_mid) ----
build_A <- function(f) {
    low  <- age_mid - age_width/2
    high <- age_mid + age_width/2

    A_dense <- matrix(0, nrow = K + 1, ncol = length(age_mid))

    for (k in 0:K) {
        a_low  <- k - 1 + f      #   [k-1+f, k+f)
        a_high <- k     + f
        overlap <- pmax(0, pmin(high, a_high) - pmax(low, a_low))
        A_dense[k + 1, ] <- overlap / age_width       # proportion of each cell
    }
    Matrix(A_dense, sparse = TRUE)
}
A_list <- lapply(f_h, build_A)

## LENGTH aggregation B :  s × m   (surveyLen × modelLen) ----
# model length bins are defined by L_edge
low_L  <- L_edge[-length(L_edge)]
high_L <- L_edge[-1]
# survey length bins are defined by S_edge
low_S  <- S_edge[-length(S_edge)]
high_S <- S_edge[-1]

B_dense <- matrix(0, nrow = s, ncol = m)
for (j in 1:s) {
    for (l in 1:m) {
        overlap <- max(0, min(high_L[l], high_S[j]) - max(low_L[l], low_S[j]))
        B_dense[j, l] <- overlap / (high_L[l] - low_L[l])   # fraction of model bin
    }
}
B <- Matrix(B_dense, sparse = TRUE)   # store once

# Model probabilities ----

# Create list `P_list` of matrices of probabilities for each quarter
# P_list[[h]] is a matrix of size (K+1) × s, where each column j
# corresponds to a survey length bin and each row k corresponds to an age class.
# Each entry P_list[[h]][k, j] is the probability p_{k|j}^{(h)} that in
# quarter h a random fis hin length bin j is of age k

# helper to compute  for one quarter
P_quarter <- function(A_h) {
    N_age   <- A_h %*% N_mat           # aggregate ages      (K+1 × m)
    N_len   <- N_age %*% t(B)          # aggregate lengths   (K+1 × s)
    totals  <- colSums(N_len)          # survey-length totals
    P <- N_len / rep(totals, each = K + 1)
    P[, totals == 0] <- 0              # guard 0/0
    P                                  # (K+1) × s
}
# Calculate probabilities for each quarter
P_list <- lapply(A_list, P_quarter)


# Survey data in matrix form ----

# Y_list[[h]] is a matrix of size (K+1) × s, where each column j
# corresponds to a survey length bin and each row k corresponds to an age class.
# Each entry Y_list[[h]][k, j] is the count of age k in length bin j
# for quarter h.

# helper to compute  for one quarter
Y_quarter <- function(h) {
    Y_mat <- matrix(0, nrow = K + 1, ncol = s)  # (K+1) × s
    y_h <- Y_obs[Y_obs$qtr == quarter_h[h], ]
    for (k in 0:K) {
        for (j in 1:s) {
            # find the count for this age and length bin
            count <- y_h$count[y_h$ageClass == k & y_h$lenBin == S_edge[j]]
            if (length(count) > 0) {
                Y_mat[k + 1, j]  <- count
            }

        }
    }
    return(Y_mat)
}
# Calculate counts for each quarter
Y_list <- lapply(seq_along(quarter_h), Y_quarter)


# Multinomial log likelihood ----
# log L = Σ y log p
loglik <- 0
for (h in seq_along(P_list)) {
    P_h <- P_list[[h]] + 1e-15
    Y_h <- Y_list[[h]]
    loglik <- loglik + sum(Y_h * log(P_h), na.rm = TRUE)
}
cat("Log‑likelihood:", loglik, "\n")

# The code above is pure R matrix algebra. For serious fitting:
#   • put the forward model (that builds each N_h) in C++/TMB
#   • keep A_list and B as constants inside the TMB object
#   • calculate log‑lik with the same vectorised trick:  yᵀ log p
#   • gradients come “for free” via automatic differentiation.


# Visualisation of the age-length key ----

# Do it for one quarter
h <- 4
P <- P_h
Y <- Y_h


## Heatmap of differences between observed and model probabilities ----

# Compute differences
diff_matrix <- apply(Y, 2, function(y) y / sum(y)) - P
diff_df <- data.frame(
    Age = rep(0:K, times = s),
    Length = rep(low_S, each = K + 1),
    difference = as.vector(diff_matrix)
)
diff_df$difference[is.nan(diff_df$difference)] <- 0  # replace NaN with 0

ggplot(diff_df, aes(x = Length, y = Age, fill = difference)) +
    geom_tile() +
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
    labs(fill = "Obs - Model") +
    theme_minimal()


## Plot quartiles of age at length ----

ages <- 0:K

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
    # Returns vector of quantiles for weighted data, using only actual observations
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cum_w <- cumsum(w) / sum(w)
    sapply(probs, function(p) {
        idx <- which(cum_w >= p)[1]
        x[idx]
    })
}

result <- data.frame(
    column = integer(0),
    type = character(0),
    median = numeric(0),
    q25 = numeric(0),
    q75 = numeric(0)
)

for (j in 1:s) {
    # Model
    mod_quants <- weighted_quantile(ages, P[,j])
    result <- rbind(result, data.frame(
        column = low_S[j],
        type = "Model",
        median = mod_quants[2],
        q25 = mod_quants[1],
        q75 = mod_quants[3]
    ))

    # Observed
    Yj <- Y[,j]
    if (sum(Yj) > 0) {
        Yj_prob <- Yj / sum(Yj)
        obs_quants <- weighted_quantile(ages, Yj_prob)
    } else {
        obs_quants <- c(NA, NA, NA)
    }
    result <- rbind(result, data.frame(
        column = low_S[j],
        type = "Observed",
        median = obs_quants[2],
        q25 = obs_quants[1],
        q75 = obs_quants[3]
    ))
}

ggplot(result, aes(x = column, y = median, color = type)) +
    geom_line() +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = type), alpha = 0.2, color = NA) +
    labs(x = "Length [cm]", y = "Age (median and IQR)", color = "Type", fill = "Type") +
    theme_minimal()




## Histograms of observed and model probabilities ----
df <- data.frame(
    Age = rep(0:K, times = s),
    Length = rep(low_S, each = K + 1),
    observed = as.vector(apply(Y, 2, function(y) y / sum(y))),
    model = as.vector(P)
)
# Convert to long format
df_long <- pivot_longer(df, cols = c("observed", "model"),
                        names_to = "Type", values_to = "Probability")

# Plot: each facet is one distribution
ggplot(df_long, aes(x = Age, y = Probability, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Length, ncol = 5) +
    theme_minimal()

# The same plot but with lines instead of bars
ggplot(df_long, aes(x = Age, y = Probability, color = Type, linetype = Type)) +
    geom_line() +
    facet_wrap(~ Length, ncol = 5) +
    theme_minimal()

