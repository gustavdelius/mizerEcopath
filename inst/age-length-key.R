
##  0.  LOAD PACKAGES AND DATA ------------------------------------
library(mizerEcopath)
library(Matrix)      # sparse matrices, fast %*%
#library(data.table)  # collapsing raw survey rows
library(dplyr)
library(tidyr)
library(ggplot2)

# Load age at size data compiled by Jess
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

params <- celtic_params
species <- "Cod"
sps <- species_params(params)[species, ]
w <- w(params)
l <- (w / sps$a)^(1/sps$b)

# Set initial condition
n_init <- rep(0, length(l))
n_init[2] <- 1  # Point source at small size
initialN(params)[species, ] <- n_init
# Solve the PDE
dt <- 0.05
nsteps <- 6/dt
n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

# Plot the evolution of n(w, t) at selected time points

# Select time points to plot
plot_steps <- c(21,41,61,81,101)
plot_times <- (plot_steps - 1) * dt
plot_labels <- paste0("age=", round(plot_times, 2))

# Prepare data frame for ggplot2
N <- length(n_hist[1, ])
w <- w(params)
l <- (w/sps$a)^(1/sps$b)
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
    xlim(1e-3, 120) +
    labs(x = "Length (cm)", y = "Abundance density",
         title = "Evolution of n(l, t) from point source at small size") +
    theme_minimal()
print(p)

##  1.  INPUTS THAT COME FROM ELSEWHERE ------------------------
# fine‑age grid (centres) and widths  [years]
age_mid   <- seq(0.05, 6,  by = 0.05)
age_width <- rep(0.05, length(age_mid))

# model length‑bin edges [cm]
L_edge <- l
# survey length‑bin edges [cm]
S_edge <- min(Y_obs$lenBin):(max(Y_obs$lenBin) + 1)

# list of fractional survey dates
f_h <- c(0.125, 0.625, 0.875)
quarter_h <- c(1, 3, 4) # quarter numbers

# list of model snapshots at those dates
# Each N_h is an (age × modelLen) matrix  dim = length(age_mid) × (length(L_edge)-1)
#
N_mat <- n_hist[2:(length(age_mid) + 1), 1:(length(L_edge) - 1)]
# Convert from densities to numbers
N_mat <- N_mat * dw(params)[1:(length(L_edge) - 1)]

K <- max(Y_obs$ageClass) # highest survey age read
H <- length(f_h)
m <- length(L_edge) - 1
s <- length(S_edge) - 1


#  2.  BUILD AGGREGATION MATRICES -----------------------------
## 2.1  AGE aggregation A_list[[h]] :  (K+1) × length(age_mid) ----
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

## 2.2  LENGTH aggregation B :  s × m   (surveyLen × modelLen) ----
low_L  <- L_edge[-length(L_edge)]
high_L <- L_edge[-1]
low_S  <- S_edge[-length(S_edge)]
high_S <- S_edge[-1]

B_dense <- matrix(0, nrow = s, ncol = m)
for (j in 1:s) {
    for (l in 1:m) {
        overlap <- max(0, min(high_L[l], high_S[j]) - max(low_L[l], low_S[j]))
        B_dense[j, l] <- overlap / (high_L[l] - low_L[l])       # fraction of model bin
    }
}
B <- Matrix(B_dense, sparse = TRUE)   # store once

#  3.  PREDICTION → AGE–LENGTH KEY IN SURVEY BINS -------------
# helper to compute p_{k|j}^{(h)} for one quarter
calc_P <- function(A_h) {
    N_age   <- A_h %*% N_mat           # aggregate ages      (K+1 × m)
    N_len   <- N_age %*% t(B)          # aggregate lengths   (K+1 × s)
    totals  <- colSums(N_len)          # survey-length totals
    P <- N_len / rep(totals, each = K + 1)
    P[, totals == 0] <- 0              # guard 0/0
    P                                  # (K+1) × s
}
P_list <- lapply(A_list, calc_P)

Y_quarter <- function(h) {
    Y_mat <- matrix(0, nrow = K + 1, ncol = s)  # (K+1) × s
    y_h <- Y_obs[Y_obs$qtr == quarter_h[h], ]
    for (k in 0:K) {
        for (j in 1:s) {
            # find the count for this age and length bin
            count <- y_h$count[y_h$ageClass == k & y_h$lenBin == S_edge[j]]
            if (length(count) > 0) {
                Y_mat[k+1, j]  <- count
            }

        }
    }
    return(Y_mat)
}
Y_list <- lapply(seq_along(quarter_h), Y_quarter)

#  4.  MULTINOMIAL LOG‑LIKELIHOOD   log L = Σ y log p ---------
loglik <- 0
for (h in seq_along(P_list)) {
    P_h <- P_list[[h]] + 1e-15
    Y_h <- Y_list[[h]]
    loglik <- loglik + sum(Y_h * log(P_h), na.rm = TRUE)
}
cat("Log‑likelihood:", loglik, "\n")


# 5. VISUALISATION OF THE AGE–LENGTH KEY -------------------

P <- P_h
Y <- Y_h

# Assume P and Y are matrices of equal dimension (n, m)
n <- nrow(P); m <- ncol(P)

# Prepare data frame
df <- data.frame(
    Age = rep(0:K, times = m),
    Length = rep(low_S, each = n),
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

ggplot(df_long, aes(x = Age, y = Probability, color = Type, linetype = Type)) +
    geom_line() +
    facet_wrap(~ Length, ncol = 5) +
    theme_minimal()

# Compute differences
diff_matrix <- apply(Y, 2, function(y) y / sum(y)) - P
diff_df <- data.frame(
    Age = rep(0:K, times = m),
    Length = rep(low_S, each = n),
    difference = as.vector(diff_matrix)
)
diff_df$difference[is.nan(diff_df$difference)] <- 0  # replace NaN with 0

ggplot(diff_df, aes(x = Length, y = Age, fill = difference)) +
    geom_tile() +
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
    labs(fill = "Obs - Model") +
    theme_minimal()


# Normalize observed counts per column
Y_norm <- apply(Y, 2, function(y) y / sum(y))  # n_bins x n_columns matrix
Y_norm[is.nan(Y_norm)] <- 0  # replace NaN with 0

# Compute summary statistics for each bin
summary_df <- data.frame(
    bin = 1:ncol(P),

    obs_mean = apply(Y_norm, 2, mean),
    obs_q25  = apply(Y_norm, 2, quantile, 0.25),
    obs_q75  = apply(Y_norm, 2, quantile, 0.75),

    mod_mean = apply(P, 2, mean),
    mod_q25  = apply(P, 2, quantile, 0.25),
    mod_q75  = apply(P, 2, quantile, 0.75)
)

# Prepare for ggplot
plot_df <- summary_df %>%
    pivot_longer(
        cols = -bin,
        names_to = c("type", ".value"),
        names_pattern = "(obs|mod)_(.+)"
    ) %>%
    mutate(type = ifelse(type == "obs", "Observed", "Model"))

ggplot(plot_df, aes(x = bin, y = mean, color = type, fill = type)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, colour = NA) +
    labs(
        x = "Bin",
        y = "Probability",
        color = "Type",
        fill = "Type",
        title = "Mean and Quartiles of Probability Distributions"
    ) +
    theme_minimal()

# Assume P and Y are matrices (n_bins x n_columns)
# vals is a numeric vector of length n_bins: the value represented by each row/bin

vals <- 0:K

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
    # Returns vector of quantiles for weighted data
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cum_w <- cumsum(w) / sum(w)
    sapply(probs, function(p) approx(cum_w, x, xout = p, ties = "ordered")$y)
}

n_col <- ncol(P)
result <- data.frame(
    column = integer(0),
    type = character(0),
    median = numeric(0),
    q25 = numeric(0),
    q75 = numeric(0)
)

for(j in 1:n_col) {
    # Model
    mod_quants <- weighted_quantile(vals, P[,j])
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
        obs_quants <- weighted_quantile(vals, Yj_prob)
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





#  WRAP IN AN OPTIMISER        (e.g. nlminb, TMB, optimx, …) ----
# The code above is pure R matrix algebra. For serious fitting:
#   • put the forward model (that builds each N_h) in C++/TMB
#   • keep A_list and B as constants inside the TMB object
#   • calculate log‑lik with the same vectorised trick:  yᵀ log p
#   • gradients come “for free” via automatic differentiation.
