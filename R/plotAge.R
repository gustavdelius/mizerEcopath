#' Plot mean age-at-length from model and survey data
#'
#' This function compares the mean age-at-length predicted by the model to
#' observed survey data for a given species. It aggregates observed
#' age-at-length data into length bins, simulates a cohort using the model, and
#' plots the mean age-at-length with confidence intervals for the observed data.
#'
#' @param params A MizerParams object.
#' @param species Name of the species to plot.
#' @param age_at_length A data frame with columns `Scientific_name`, `Quarter`,
#'   `LngtClass`, `Age`, and `CANoAtLngt` giving observed age-at-length data.
#' @param dt Time step for the model simulation (years). Default is 0.05.
#' @param t_max Maximum age to simulate (years). Default is 6.
#' @param f_h Numeric vector of time of surveys in the different quarters.
#'   Default c(0.125, 0.625, 0.875).
#' @param quarter_h Numeric vector of survey quarters. Default is
#'   c(1, 3, 4).
#' @param plot One of "mean" or "quantiles". If "mean", the mean age-at-length
#'   is plotted. If "quantiles", the 25%, 50% and 75% quartiles are plotted.
#'
#' @return A ggplot2 object showing mean or quantiles of the age-at-length for
#'   model and observed data, with 95% confidence intervals for the observed
#'   means.
#'
#' @import ggplot2
#' @import dplyr
#' @examples
#' # See package vignette for example usage
#' @export
plotAge <- function(params, species, age_at_length,
                    dt = 0.05, t_max = NULL,
                    f_h = c(0.125, 0.625, 0.875),
                    quarter_h = c(1, 3, 4),
                    plot = c("mean", "quantiles")) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)
    plot <- match.arg(plot)

    sci_name <- species_params(params)[species, "SciName"]

    # Aggregate to cm bins, preserving quarter information
    Y_obs <- age_at_length |>
        filter(Scientific_name == sci_name) |>
        # remove rows with NA in any column
        filter(!is.na(LngtClass) & !is.na(Age) & !is.na(CANoAtLngt)) |>
        # round down to cm
        mutate(LngtClass = floor(LngtClass)) |>
        # aggregate counts by quarter
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
    if (is.null(t_max)) {
        t_max <- K + 2  # default maximum time
    }
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

    # Survey data in matrix form by quarter ----

    # helper to compute survey data for one quarter
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

    # Model: aggregate N_len across quarters, weighted by observed counts ----

    # helper to compute model probabilities for one quarter
    P_quarter <- function(A_h) {
        N_age   <- A_h %*% N_mat           # aggregate ages      (K+1 × m)
        N_len   <- N_age %*% t(B)          # aggregate lengths   (K+1 × s)
        totals  <- colSums(N_len)          # survey-length totals
        P <- N_len / rep(totals, each = K + 1)
        P[, totals == 0] <- NA              # guard 0/0
        P                                  # (K+1) × s
    }
    # Calculate probabilities for each quarter
    P_list <- lapply(A_list, P_quarter)

    # Weight model counts by observed counts in each quarter and size bin
    N_len_weighted <- matrix(0, nrow = K + 1, ncol = s)
    for (h in seq_along(P_list)) {
        P_h <- P_list[[h]]
        Y_h <- Y_list[[h]]
        # For each length bin, weight the model probabilities by observed counts
        for (j in 1:s) {
            total_obs_count <- sum(Y_h[, j])
            if (total_obs_count > 0) {
                N_len_weighted[, j] <- N_len_weighted[, j] + P_h[, j] * total_obs_count
            }
        }
    }

    # Convert to probabilities
    totals <- colSums(N_len_weighted)
    #P <- N_len_weighted / rep(totals, each = K + 1)
    P <- sweep(N_len_weighted, 2, totals, "/")
    P[, totals == 0] <- NA

    if (plot == "mean") {
        # Prepare data for plotting ----
        all_results <- data.frame()
        for (j in 1:s) {
            # Model (mean only, no CI)
            mod_mean <- sum(ages * P[,j])
            all_results <- rbind(all_results, data.frame(
                column = low_S[j],
                type = "Model",
                mean = mod_mean,
                lower = NA,
                upper = NA
            ))
            # Observed (mean and multinomial CI) - aggregate across all quarters
            Yj <- rowSums(sapply(Y_list, function(Y_h) Y_h[, j]))
            if (sum(Yj) > 0) {
                obs_stats <- obs_mean_ci(ages, Yj)
            } else {
                obs_stats <- c(NA, NA, NA)
            }
            all_results <- rbind(all_results, data.frame(
                column = low_S[j],
                type = "Observed",
                mean = obs_stats[1],
                lower = obs_stats[2],
                upper = obs_stats[3]
            ))
        }

        ggplot(all_results, aes(x = column, y = mean, color = type)) +
            geom_line() +
            geom_ribbon(data = subset(all_results, type == "Observed"),
                        aes(ymin = lower, ymax = upper, fill = type),
                        alpha = 0.2, color = NA) +
            labs(x = "Length [cm]", y = "Age (mean and 95% CI)", color = "Type", fill = "Type") +
            theme_minimal()
    } else {

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
            Yj <- rowSums(sapply(Y_list, function(Y_h) Y_h[, j]))
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

    }
}


#' Compute mean and confidence interval for mean from multinomial counts
#'
#' Given a vector of values and corresponding counts, computes the mean and an approximate confidence interval for the mean
#' using a normal approximation to the multinomial distribution.
#'
#' @param x Vector of values (e.g., ages).
#' @param counts Vector of counts for each value.
#' @param conf Confidence level for the interval. Default is 0.95.
#'
#' @return A named numeric vector with elements `mean`, `lower`, and `upper`.
#' @keywords internal
obs_mean_ci <- function(x, counts, conf = 0.95) {
    n <- sum(counts)
    if (n == 0) return(c(mean = NA, lower = NA, upper = NA))
    p <- counts / n
    m <- sum(x * p)
    v <- sum(p * (x - m)^2) / n
    se <- sqrt(v)
    z <- qnorm(1 - (1 - conf) / 2)
    c(mean = m, lower = m - z * se, upper = m + z * se)
}


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
