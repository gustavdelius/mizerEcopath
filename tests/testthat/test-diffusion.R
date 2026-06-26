# test-diffusion.R
# Tests that matchCatch()'s TMB objective solves the same diffusion-aware
# steady state as mizer, and that D_ext is optimised correctly.
library(testthat)
library(mizer)

# Single-gear version of the example model, consistent with `celtic_catch`.
make_celtic_single_gear <- function() {
    params <- celtic_params
    gp <- gear_params(params)
    keep <- do.call(rbind, lapply(split(gp, gp$species), function(g) {
        row <- g[which.max(g$yield_observed), , drop = FALSE]
        row$gear <- "total"
        row$yield_observed <- sum(g$yield_observed)
        row
    }))
    rownames(keep) <- paste(keep$species, keep$gear, sep = ", ")
    gear_params(params) <- keep
    initial_effort(params) <- 1
    params
}

# Catch data sampled from the model's predicted catch density (N * F_mort).
model_catch_data <- function(params, species, gear = "total", n_total = 1e5) {
    sp <- species_params(params)
    s <- sp$species == species
    sps <- sp[s, ]
    lengths <- (params@w / sps$a)^(1 / sps$b)
    l_max <- (sps$w_max / sps$a)^(1 / sps$b)
    rate <- params@initial_n[s, ] * getFMort(params)[s, ] * params@dw
    df <- aggregate(rate, list(length = floor(lengths)), sum)
    names(df)[2] <- "rate"
    df <- df[df$rate > 0 & (df$length + 1) <= l_max, ]
    count <- round(df$rate / sum(df$rate) * n_total)
    keep <- count > 0
    data.frame(species = species, gear = gear, length = df$length[keep],
               dl = 1, catch = count[keep])
}

# Build the TMB objective for a single (sigmoid) gear species, at parameter
# values matching the model except for the supplied log_D_ext.
build_obj <- function(params, species, catch, log_D_ext) {
    data <- prepare_data(params, species = species, catch)
    gps <- gear_params(params)[gear_params(params)$species == species, ]
    s <- species_params(params)$species == species
    sps <- species_params(params)[s, ]
    mat_idx <- sum(params@w < sps$w_mat)
    mu_mat <- ext_mort(params)[s, mat_idx]
    ip <- list(
        logit_l50 = qlogis((gps$l50 - min(data$l)) / (max(data$l) - min(data$l))),
        log_ratio_left = qlogis((gps$l50 - gps$l25) / gps$l50),
        log_l50_right_offset = 1,
        log_ratio_right = 1,
        log_catchability = log(gps$catchability),
        mu_mat = mu_mat,
        m = sps$m,
        log_D_ext = log_D_ext
    )
    map <- list(log_l50_right_offset = factor(NA), log_ratio_right = factor(NA))
    obj <- TMB::MakeADFun(data = data, parameters = ip, map = map,
                          DLL = "mizerEcopath", silent = TRUE)
    obj$fn(obj$par)  # populate last.par / report
    obj
}

# Independent R implementation of the diffusion-aware upwind steady-state
# tridiagonal solve (mizer's get_transport_coefs_upwind at dt = 1, Thomas solve),
# given growth g_j, mortality mu_j, diffusion d_j, bin widths dw and egg density N0.
r_steady_diffusion <- function(g, mu, d, dw, N0) {
    K <- length(dw)
    a <- b <- cc <- S <- numeric(K)
    a[1] <- 0; b[1] <- 1; cc[1] <- 0; S[1] <- N0
    for (j in 2:K) {
        inv <- 1 / dw[j]
        a[j] <- -inv * (g[j - 1] + 0.5 * d[j - 1] / dw[j - 1])
        b[j] <- mu[j] + inv * (g[j] + 0.5 * d[j] / dw[j] + 0.5 * d[j] / dw[j - 1])
        cc[j] <- if (j < K) -inv * 0.5 * d[j + 1] / dw[j] else 0
        S[j] <- 0
    }
    cp <- sp <- N <- numeric(K)
    cp[1] <- cc[1] / b[1]; sp[1] <- S[1] / b[1]
    for (j in 2:K) {
        den <- b[j] - a[j] * cp[j - 1]
        cp[j] <- cc[j] / den
        sp[j] <- (S[j] - a[j] * sp[j - 1]) / den
    }
    N[K] <- sp[K]
    for (j in (K - 1):1) N[j] <- sp[j] - cp[j] * N[j + 1]
    N
}

test_that("C++ objective solves the documented diffusion steady-state scheme", {
    p <- make_celtic_single_gear()
    species <- "Cod"
    catch <- model_catch_data(p, species)

    # Test across a range of diffusion strengths, including ~zero diffusion.
    for (log_D in c(-40, log(1), log(5.5), log(50))) {
        obj <- build_obj(p, species, catch, log_D)
        rep <- obj$report()
        N_r <- r_steady_diffusion(rep$growth, rep$mort, rep$d_diff,
                                  obj$env$data$dw, rep$N[1])
        expect_equal(rep$N, N_r, tolerance = 1e-9,
                     info = paste("log_D_ext =", log_D))
    }
})

test_that("with D_ext ~ 0 the scheme reduces to the pure-advection recursion", {
    p <- make_celtic_single_gear()
    species <- "Cod"
    catch <- model_catch_data(p, species)
    obj <- build_obj(p, species, catch, log_D_ext = -40)
    rep <- obj$report()
    g <- rep$growth; mu <- rep$mort; dw <- obj$env$data$dw
    N_adv <- numeric(length(dw))
    N_adv[1] <- rep$N[1]
    for (i in 2:length(dw)) {
        N_adv[i] <- N_adv[i - 1] * g[i - 1] / (g[i] + mu[i] * dw[i])
    }
    expect_equal(rep$N, N_adv, tolerance = 1e-8)
})

test_that("matchCatch optimises D_ext by default and writes it back", {
    p <- make_celtic_single_gear()
    species <- "Cod"
    catch <- model_catch_data(p, species)

    D_ext_before <- diffusion_coefficient(
        p, species_params(p)$species == species,
        species_params(p)$n[species_params(p)$species == species])

    result <- suppressWarnings(matchCatch(p, species = species, catch = catch))
    sp_res <- species_params(result)[species_params(result)$species == species, ]

    # D_ext is a fitted parameter, recorded in species_params, and moved from start
    expect_true("D_ext" %in% names(sp_res))
    expect_false(isTRUE(all.equal(sp_res$D_ext, D_ext_before)))
    # respects the optimiser bounds exp([-20, 15])
    expect_gte(sp_res$D_ext, exp(-20))
    expect_lte(sp_res$D_ext, exp(15))
    # the ext_diffusion slot is kept consistent with the fitted D_ext
    s <- species_params(result)$species == species
    expect_equal(unname(result@ext_diffusion[s, ]),
                 unname(sp_res$D_ext * result@w^(sp_res$n + 1)),
                 tolerance = 1e-8)
})

test_that("map fixes D_ext at its initial value", {
    p <- make_celtic_single_gear()
    species <- "Cod"
    catch <- model_catch_data(p, species)

    D_ext_before <- diffusion_coefficient(
        p, species_params(p)$species == species,
        species_params(p)$n[species_params(p)$species == species])

    result <- suppressWarnings(matchCatch(p, species = species, catch = catch,
                                          map = list(D_ext = factor(NA))))
    sp_res <- species_params(result)[species_params(result)$species == species, ]
    expect_equal(sp_res$D_ext, D_ext_before, tolerance = 1e-6)
})
