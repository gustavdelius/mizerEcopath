
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
species <- "Herring"
sps <- species_params(params)[species, ]
w <- w(params)
l <- (w / sps$a)^(1/sps$b)

# Set initial condition
n_init <- rep(0, length(l))
n_init[2] <- 1  # Point source at small size
initialN(params)[species, ] <- n_init
# Solve the PDE
dt <- 0.05
nsteps <- 200
n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

##  1.  INPUTS THAT COME FROM ELSEWHERE ------------------------
# fine‑age grid (centres) and widths  [years]
age_mid   <- seq(0.05, 10,  by = 0.05)
age_width <- rep(0.05, length(age_mid))

# model length‑bin edges [cm]
L_edge <- l
# survey length‑bin edges [cm]
S_edge <- min(Y_obs$lenBin):(max(Y_obs$lenBin) + 1)

# list of fractional survey dates
f_h <- c(0.125, 0.375, 0.625, 0.875)             # one per quarter

# list of model snapshots at those dates
# Each N_h is an (age × modelLen) matrix  dim = length(age_mid) × (length(L_edge)-1)
#
N_mat <- n_hist[2:(length(age_mid) + 1), 1:(length(L_edge) - 1)]
# Convert from densities to numbers
N_mat <- N_mat * w[1:(length(L_edge) - 1)]

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
    P <- N_len / rep(totals, each = K+1)
    P[, totals == 0] <- 0              # guard 0/0
    P                                           # (K+1) × s
}
P_list <- lapply(A_list, calc_P)


#  4.  MULTINOMIAL LOG‑LIKELIHOOD   log L = Σ y log p ---------
loglik <- 0
for (h in seq_along(P_list)) {
    P_h <- P_list[[h]]
    # merge counts for this quarter
    y_h <- Y_obs[Y_obs$qtr == h, ]
    # vector‑index into P_h  (convert ageClass 0..K, lenBin 1..s)
    idx <- (y_h$ageClass + 1) + (y_h$lenBin - 1) * (K + 1)
    p_vec <- as.vector(P_h)[idx] + 1e-15 # guard 0/0
    y_vec <- y_h$count
    loglik <- loglik + sum(y_vec * log(p_vec), na.rm = TRUE)
}
cat("Log‑likelihood:", loglik, "\n")


# 5. VISUALISATION OF THE AGE–LENGTH KEY -------------------

## build one long data-frame of model probabilities
model_long <- lapply(seq_along(P_list), function(h) {
    P <- as.matrix(P_list[[h]])
    colnames(P) <- 1:ncol(P)
    rownames(P) <- 0:K
    as.data.frame(P) |>
        mutate(ageClass = 0:K, qtr = h) |>
        pivot_longer(-c(ageClass, qtr), names_to = "lenBin", values_to = "prob")
}) |> bind_rows() |> mutate(source = "Model") |>
    select(ageClass, lenBin, prob, source, qtr)

# Check that the model probabilities sum to 1 per length bin
model_check <- model_long |>
    group_by(qtr, lenBin) |>
    summarise(prob = sum(prob)) |>
    ungroup()


## add observed proportions (counts normalised to 1 per length bin) -------
obs_long <- Y_obs |>
    group_by(qtr, lenBin) |>
    mutate(prob = count / sum(count)) |>
    ungroup() |>
    mutate(source = "Observed") |>
    select(ageClass, lenBin, prob, source, qtr)

plot_df <- rbind(model_long, obs_long)

ggplot(plot_df,
       aes(x = factor(ageClass), y = prob,
           fill = source, alpha = source)) +
    geom_col(position = "dodge", width = .9) +
    facet_grid(lenBin ~ qtr, labeller = "label_both") +
    scale_alpha_manual(values = c(0.8, 0.4), guide = "none") +
    scale_fill_manual(values = c(Model = "grey20", Observed = "steelblue")) +
    labs(x = "Survey age class", y = "Proportion within length bin",
         fill = "", title = "Age–length keys: model vs observed")

library(ggridges)
ggplot(model_long[model_long$source=="Model", ],
       aes(x = ageClass, y = factor(lenBin), height = prob,
           group = lenBin)) +
    geom_ridgeline(fill = "grey80", colour = "grey30", alpha = .6) +
    geom_ridgeline(data = obs_long,   # overlay observed outline
                   colour = "steelblue", size = .3, fill = NA) +
    labs(y = "Length bin", x = "Age class",
         title = "Model (filled) vs observed (outline) age distributions")


#  WRAP IN AN OPTIMISER        (e.g. nlminb, TMB, optimx, …) ----
# The code above is pure R matrix algebra. For serious fitting:
#   • put the forward model (that builds each N_h) in C++/TMB
#   • keep A_list and B as constants inside the TMB object
#   • calculate log‑lik with the same vectorised trick:  yᵀ log p
#   • gradients come “for free” via automatic differentiation.
