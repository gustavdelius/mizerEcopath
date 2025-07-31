
##  0.  LOAD PACKAGES AND DATA ------------------------------------
# (Matrix for sparse algebra, data.table for counts)
library(Matrix)      # sparse matrices, fast %*%
library(data.table)  # collapsing raw survey rows
library(dplyr)

# Load age at size data compiled by Jess
survey <- readRDS(here::here("inst", "extdata",
                             "Celtic_Sea_Size_at_Age_Data.rds"))
# Aggregate to cm bins
survey <- survey |>
    # remove rows with NA in any column
    filter(!is.na(LngtClass) & !is.na(Age) & !is.na(CANoAtLngt)) |>
    # round down to cm
    mutate(LngtClass = floor(LngtClass)) |>
    # aggregate counts
    group_by(Quarter, LngtClass, Age) |>
    summarise(CANoAtLngt = sum(CANoAtLngt, na.rm = TRUE), .groups = "drop")

params <- celtic_params
species <- "Herring"
sps <- species_params(params)[species, ]
w <- w(params)
l <- sps$a * w^sps$b

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
S_edge <- min(survey$LntClass):(max(survey$LngtClass) + 1)

# list of fractional survey dates   f_h  (e.g. Q1 = 0.125, Q2 = 0.375 …)
f_h <- c(0.125, 0.375, 0.625, 0.875)             # one per quarter

# list of model snapshots at those dates
# Each N_h is an (age × modelLen) matrix  dim = length(age_mid) × (length(L_edge)-1)
#
N_list <- n_hist[2:(length(age_mid) + 1), 1:(length(L_edge) - 1)]

# observed survey counts collapsed to a data.table
#  columns: qtr (1..H), lenBin (1..s), ageClass (0..K), count
Y_obs <- survey |>
    transmute(qtr = Quarter,
              lenBin = as.integer(LngtClass),
              ageClass = as.integer(Age),
              count = CANoAtLngt)

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
predict_ALK <- function(N_h, A_h) {
    #  Step 1: aggregate ages
    N_age <- A_h %*% N_h                     # (K+1) × m
    #  Step 2: aggregate lengths
    N_len <- N_age %*% t(B)                  # (K+1) × s
    #  Step 3: conditional age distribution at each survey length bin
    colSums_N <- colSums(N_len)              # length‑specific totals
    P <- N_len / rep(colSums_N, each = K + 1)  # sweep but fast
    P[, colSums_N == 0] <- 0                 # guard 0/0
    P                                         # matrix (K+1) × s
}

P_list <- mapply(predict_ALK, N_list, A_list, SIMPLIFY = FALSE)


#  4.  MULTINOMIAL LOG‑LIKELIHOOD   log L = Σ y log p ---------
loglik <- 0
for (h in 1:H) {
    P_h <- P_list[[h]]
    # merge counts for this quarter
    y_h <- Y_obs[Y_obs$qtr == h, ]
    # vector‑index into P_h  (convert ageClass 0..K, lenBin 1..s)
    idx <- (y_h$ageClass + 1) + (y_h$lenBin - 1) * (K + 1)
    p_vec <- as.vector(P_h)[idx]     # fast linear access
    y_vec <- y_h$count
    if (any(p_vec == 0 & y_vec > 0))
        return(-Inf)                   # impossible cell
    loglik <- loglik + sum(y_vec * log(p_vec))
}
cat("Log‑likelihood:", loglik, "\n")

#  5.  WRAP IN AN OPTIMISER        (e.g. nlminb, TMB, optimx, …) ----
# The code above is pure R matrix algebra. For serious fitting:
#   • put the forward model (that builds each N_h) in C++/TMB
#   • keep A_list and B as constants inside the TMB object
#   • calculate log‑lik with the same vectorised trick:  yᵀ log p
#   • gradients come “for free” via automatic differentiation.
