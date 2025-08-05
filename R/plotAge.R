plotAge <- function(params, species, age_at_length,
                    dt = 0.05, t_max = 6,
                    f_h = c(0.125, 0.625, 0.875),
                    quarter_h = c(1, 3, 4), h = 3) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)

    # Aggregate to cm bins
    Y_obs <- age_at_length |>
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
    ages <- 0:K

    # Evolve a cohort over time in model ----

    ## Set up example ----
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
    age_mid   <- seq(dt, t_max,  by = dt)
    nsteps <- length(age_mid)
    age_width <- rep(dt, nsteps)

    n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

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
        # Matrix::Matrix(A_dense, sparse = TRUE)
        A <- A_dense
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
    #B <- Matrix::Matrix(B_dense, sparse = TRUE)   # store once
    B <- B_dense

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

    # Panelled plot for all quarters
    all_results <- data.frame()
    for (i in seq_along(quarter_h)) {
        P <- P_list[[i]]
        Y <- Y_list[[i]]
        for (j in 1:s) {
            # Model
            mod_quants <- weighted_quantile(ages, P[,j])
            all_results <- rbind(all_results, data.frame(
                column = low_S[j],
                type = "Model",
                median = mod_quants[2],
                q25 = mod_quants[1],
                q75 = mod_quants[3],
                quarter = as.factor(quarter_h[i])
            ))
            # Observed
            Yj <- Y[,j]
            if (sum(Yj) > 0) {
                Yj_prob <- Yj / sum(Yj)
                obs_quants <- weighted_quantile(ages, Yj_prob)
            } else {
                obs_quants <- c(NA, NA, NA)
            }
            all_results <- rbind(all_results, data.frame(
                column = low_S[j],
                type = "Observed",
                median = obs_quants[2],
                q25 = obs_quants[1],
                q75 = obs_quants[3],
                quarter = as.factor(quarter_h[i])
            ))
        }
    }

    ggplot(all_results, aes(x = column, y = median, color = type)) +
        geom_line() +
        geom_ribbon(aes(ymin = q25, ymax = q75, fill = type), alpha = 0.2, color = NA) +
        labs(x = "Length [cm]", y = "Age (median and IQR)", color = "Type", fill = "Type") +
        theme_minimal() +
        facet_wrap(~ quarter, ncol = 1)

}
