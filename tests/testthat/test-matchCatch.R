# test-matchCatch.R
#
# Tests for matchCatch(), which fits gear selectivity, catchability and
# mortality so the model reproduces observed catch size distributions,
# yield and production.
#
# Input validation:
#   - errors when the catch contains a gear not in the model
#   - errors when catch data is missing required columns
#   - warns on non-existent species and errors when none are selected
#
# Core fitting behaviour:
#   - runs without error for a single species
#   - recovers a model-generated catch size distribution
#   - recovers the true parameters from a self-consistent catch
#   - recovers selectivity and catchability from the catch alone
#   - fits only the catch's gears, leaving other gears untouched
#   - updates gear selectivity parameters for the selected species
#   - preserves biomass for the adjusted species
#   - handles multiple species
#   - works with multiple gears, incl. double_sigmoid_length
#
# Optional / free parameters:
#   - respects the yield_lambda parameter
#   - fits m when freed via map (m > n)
#   - moves D_ext from its starting value when free
#   - mu_mat_lim caps the optimised external mortality
#   - map fixes a specified parameter at its initial value
#
# Missing-data and edge cases:
#   - warns and returns params unchanged with no catch and no production
#   - still runs when yield_observed is missing
#   - still runs when production_observed is missing
#   - with empty catch, matches only yield and production (selectivity
#     left unchanged)

library(testthat)
library(mizer)

# Collapse celtic_params to a single gear per species. Used only by the
# self-contained "recovers a model-generated distribution" test below, where the
# catch is generated from this same model.
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

# Expected catch in each 1-cm length bin: the integral of the catch density
# (N * F_mort) over the bin, using the same piecewise-linear interpolation of
# the density that matchCatch()'s likelihood uses internally. Generating the
# synthetic data with the *same* binning as the likelihood is what makes the
# model self-consistent, so that a fit can recover the parameters that
# generated it. A cruder floor()-histogram disagrees by a few percent, and
# because the mortality-diffusion valley is so flat that small discrepancy is
# enough to shift the recovered mortality substantially.
model_catch_rate <- function(params, species, bins) {
    sp <- species_params(params)
    s <- which(sp$species == species)
    sps <- sp[s, ]
    w <- params@w
    dens <- approxfun(w, params@initial_n[s, ] * getFMort(params)[s, ], rule = 2)
    vapply(bins, function(L) {
        lo <- sps$a * L^sps$b
        hi <- sps$a * (L + 1)^sps$b
        grid <- sort(unique(c(w[w > lo & w < hi], lo, hi)))
        v <- dens(grid)
        sum(diff(grid) * (head(v, -1) + tail(v, -1)) / 2)
    }, numeric(1))
}

# Generate catch data from the model itself, in the `celtic_catch` format
# (integer counts over 1-cm length bins summing to ~n_total).
model_catch_data <- function(params, species, gear = "total", n_total = 1e5) {
    sp <- species_params(params)
    s <- sp$species == species
    sps <- sp[s, ]
    w <- params@w
    dens <- params@initial_n[s, ] * getFMort(params)[s, ]
    l_max <- (sps$w_max / sps$a)^(1 / sps$b)
    bins <- seq(floor((min(w[dens > 0]) / sps$a)^(1 / sps$b)), floor(l_max) - 1)
    rate <- model_catch_rate(params, species, bins)
    bins <- bins[rate > 0]
    rate <- rate[rate > 0]
    count <- round(rate / sum(rate) * n_total)
    data.frame(species = species, gear = gear, length = bins[count > 0],
               dl = 1, catch = count[count > 0])
}

# Relative frequency of the model's predicted catch over the given 1-cm bins.
model_catch_hist <- function(params, species, length_bins) {
    rate <- model_catch_rate(params, species, length_bins)
    rate / sum(rate)
}

# Build params with two gears for Hake (sigmoid_length + double_sigmoid_length)
# and matching two-gear catch data.
make_two_gear_hake <- function() {
    params <- celtic_params
    gp <- gear_params(params)
    gp$l50_right <- NA_real_
    gp$l25_right <- NA_real_

    hake_row <- gp["Hake, Gillnet", ]

    gear1 <- hake_row
    gear1$gear <- "small_mesh"
    gear1$yield_observed <- hake_row$yield_observed * 0.4

    gear2 <- hake_row
    gear2$gear <- "large_mesh"
    gear2$sel_func <- "double_sigmoid_length"
    gear2$l50_right <- 110
    gear2$l25_right <- 120
    gear2$yield_observed <- hake_row$yield_observed * 0.6

    gear_params(params) <- rbind(gp[gp$species != "Hake", ], gear1, gear2)

    # Aggregate the per-gear Hake catch into a single size distribution and
    # rescale to integer-like counts (celtic_catch holds small fractional
    # densities), then split it between the two new gears.
    hake_catch <- celtic_catch[celtic_catch$species == "Hake", ]
    hake_catch$catch <- ave(hake_catch$catch, hake_catch$length, FUN = sum)
    hake_catch <- hake_catch[!duplicated(hake_catch$length), ]
    hake_catch$catch <- round(hake_catch$catch / sum(hake_catch$catch) * 1e5)
    small_catch <- hake_catch
    small_catch$gear <- "small_mesh"
    small_catch$catch <- round(hake_catch$catch * 0.4)
    large_catch <- hake_catch
    large_catch$gear <- "large_mesh"
    large_catch$catch <- round(hake_catch$catch * 0.6)

    list(params = params, catch = rbind(small_catch, large_catch))
}

test_that("matchCatch errors when the catch contains a gear not in the model", {
    bad_catch <- celtic_catch
    bad_catch$gear[bad_catch$species == "Hake"] <- "phantom_gear"
    expect_error(
        matchCatch(celtic_params, species = "Hake", catch = bad_catch),
        "are not among the model's gears"
    )
})

test_that("matchCatch fits only the catch's gears, leaving others untouched", {
    gp_before <- gear_params(celtic_params)
    result <- matchCatch(celtic_params, species = "Cod", catch = celtic_catch) |>
        suppressWarnings()
    gp_after <- gear_params(result)

    cod <- gp_before$species == "Cod"
    catch_gears <- unique(celtic_catch$gear[celtic_catch$species == "Cod"])
    fished <- gp_before$gear[cod] %in% catch_gears

    # The survey gear (no catch data) is left exactly as it was.
    untouched <- gp_before[cod, ][!fished, c("l50", "l25", "catchability")]
    after_untouched <- gp_after[cod, ][!fished, c("l50", "l25", "catchability")]
    expect_equal(untouched, after_untouched)

    # The commercial gears (with catch data) have been refitted.
    expect_false(isTRUE(all.equal(
        gp_before[cod, ][fished, c("l50", "l25", "catchability")],
        gp_after[cod, ][fished, c("l50", "l25", "catchability")])))
})

test_that("matchCatch runs without error with valid inputs for a single species", {
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    expect_s4_class(result, "MizerParams")
})

test_that("matchCatch recovers a model-generated catch size distribution", {
    # Generate catch data from a single-gear model, then check that the fitted
    # model reproduces that size distribution closely. (Recovery of the
    # parameters themselves is tested separately below.)
    sg_params <- make_celtic_single_gear()
    species <- "Cod"
    catch <- model_catch_data(sg_params, species)
    result <- matchCatch(sg_params, species = species, catch = catch) |>
        suppressWarnings()

    observed <- catch$catch / sum(catch$catch)
    modelled <- model_catch_hist(result, species, catch$length)
    tv_distance <- 0.5 * sum(abs(observed - modelled))
    expect_lt(tv_distance, 0.05)
})

test_that("matchCatch recovers the true parameters from a self-consistent catch", {
    # With the catch generated using the same bin-integration as the likelihood,
    # and the yield and production supplied, matchCatch() recovers the
    # parameters that generated the data even when started far from them. The
    # data is generated from the single-species steady state, which the fit's
    # solver can reproduce exactly.
    species <- "Haddock"
    sg <- steadySingleSpecies(make_celtic_single_gear(), species = species)
    s  <- species_params(sg)$species == species
    catch <- model_catch_data(sg, species)

    gp_true   <- gear_params(sg)[gear_params(sg)$species == species, ]
    l50_true  <- gp_true$l50
    q_true    <- gp_true$catchability
    mu_true   <- species_params(sg)[s, "mu_mat"]
    D_true    <- species_params(sg)[s, "D_ext"]
    prod_true <- getProduction(sg)[species]

    # Start from deliberately wrong values: selectivity 8 cm too high,
    # mortality 4x too high, diffusion 4x too low, catchability reset to 1.
    pp <- sg
    gp <- gear_params(pp); gi <- which(gp$species == species)
    gp$l50[gi] <- l50_true + 8; gp$l25[gi] <- gp$l50[gi] - 4
    gp$catchability[gi] <- 1
    gear_params(pp) <- gp
    spp <- species_params(pp)
    spp$mu_mat[s] <- mu_true * 4
    spp$D_ext[s]  <- D_true / 4
    spp$production_observed[s] <- prod_true
    pp@ext_diffusion[s, ] <- spp$D_ext[s] * pp@w^(spp$n[s] + 1)
    species_params(pp) <- spp

    fit <- suppressWarnings(matchCatch(
        pp, species = species, catch = catch,
        yield_lambda = 1, production_lambda = 1, map = list(m = factor(NA))))
    gp_fit <- gear_params(fit)[gear_params(fit)$species == species, ]
    sp_fit <- species_params(fit)[s, ]

    expect_equal(gp_fit$l50, l50_true, tolerance = 0.02)         # selectivity
    expect_equal(gp_fit$catchability, q_true, tolerance = 0.03) # from the yield
    expect_equal(sp_fit$mu_mat, mu_true, tolerance = 0.12)      # mortality
    expect_equal(sp_fit$D_ext, D_true, tolerance = 0.12)        # diffusion
})

test_that("matchCatch recovers selectivity and catchability from the catch alone", {
    # Even without a production constraint, the catch shape determines the
    # selectivity and the yield determines the catchability. (Mortality and
    # diffusion are recovered too, but only loosely without production - see the
    # catch_identifiability vignette - so they are not asserted tightly here.)
    species <- "Haddock"
    sg <- steadySingleSpecies(make_celtic_single_gear(), species = species)
    s  <- species_params(sg)$species == species
    catch <- model_catch_data(sg, species)

    gp_true  <- gear_params(sg)[gear_params(sg)$species == species, ]
    l50_true <- gp_true$l50
    q_true   <- gp_true$catchability

    pp <- sg
    gp <- gear_params(pp); gi <- which(gp$species == species)
    gp$l50[gi] <- l50_true + 8; gp$l25[gi] <- gp$l50[gi] - 4
    gp$catchability[gi] <- 1
    gear_params(pp) <- gp
    spp <- species_params(pp)
    spp$mu_mat[s] <- spp$mu_mat[s] * 4
    species_params(pp) <- spp

    fit <- suppressWarnings(matchCatch(
        pp, species = species, catch = catch,
        yield_lambda = 1, production_lambda = 0, map = list(m = factor(NA))))
    gp_fit <- gear_params(fit)[gear_params(fit)$species == species, ]

    expect_equal(gp_fit$l50, l50_true, tolerance = 0.02)
    expect_equal(gp_fit$catchability, q_true, tolerance = 0.03)
})

test_that("matchCatch throws error if catch data is missing required columns", {
    bad_catch_data <- celtic_catch
    bad_catch_data$dl <- NULL

    expect_error(
        matchCatch(celtic_params, species = "Hake", catch = bad_catch_data),
        "must contain columns"
    )
})

test_that("matchCatch can handle multiple species", {
    multi_species <- c("Hake", "Herring")
    result <- matchCatch(celtic_params, species = multi_species, catch = celtic_catch) |>
        suppressWarnings()

    expect_s4_class(result, "MizerParams")
    gp <- gear_params(result)
    expect_true(all(multi_species %in% gp$species))
})

test_that("matchCatch preserves biomass for the adjusted species", {
    biomass_before <- getBiomass(celtic_params, use_cutoff = TRUE)
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    biomass_after <- getBiomass(result, use_cutoff = TRUE)

    expect_equal(biomass_before["Hake"], biomass_after["Hake"], tolerance = 1e-3)
})

test_that("matchCatch updates gear selectivity parameters for the selected species", {
    gp_before <- gear_params(celtic_params)
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    gp_after <- gear_params(result)

    sp_idx <- which(gp_after$species == "Hake")
    expect_false(all(gp_before[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_after[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch respects the yield_lambda parameter", {
    result_default <- matchCatch(celtic_params, species = "Hake",
                                 catch = celtic_catch) |>
        suppressWarnings()
    gp_default <- gear_params(result_default)

    result_mod <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch,
                             yield_lambda = 10) |>
        suppressWarnings()
    gp_mod <- gear_params(result_mod)

    sp_idx <- which(gp_mod$species == "Hake")
    expect_false(all(gp_default[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_mod[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch handles invalid species inputs gracefully", {
    # Muffle the unrelated `erepro` adjustment warning emitted while fitting Cod
    # so that only the invalid-species warning is asserted here.
    expect_warning(
        withCallingHandlers(
            matchCatch(celtic_params, species = c("Cod", "NotASpecies"),
                       catch = celtic_catch),
            warning = function(w) {
                if (!grepl("do not exist", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "The following species do not exist: NotASpecies."
    )

    expect_error(
        matchCatch(celtic_params, species = character(0), catch = celtic_catch),
        "No species have been selected."
    )
})

test_that("matchCatch works with multiple gears including double_sigmoid_length", {
    tg <- make_two_gear_hake()
    result <- suppressWarnings(matchCatch(tg$params, species = "Hake",
                                          catch = tg$catch))
    expect_s4_class(result, "MizerParams")
    gp <- gear_params(result)[gear_params(result)$species == "Hake", ]
    expect_equal(nrow(gp), 2)

    sigmoid_row <- gp[gp$gear == "small_mesh", ]
    double_row  <- gp[gp$gear == "large_mesh", ]

    expect_false(all(sigmoid_row[, c("l50", "l25", "catchability")] ==
                     gear_params(tg$params)[gear_params(tg$params)$gear == "small_mesh",
                                            c("l50", "l25", "catchability")]))
    expect_true(!is.na(double_row$l50_right))
    expect_true(!is.na(double_row$l25_right))
    expect_true(is.na(sigmoid_row$l50_right))
})

test_that("matchCatch with m free produces a valid result", {
    result <- suppressWarnings(matchCatch(celtic_params, species = "Hake",
                                          catch = celtic_catch,
                                          map = list(m = NULL)))
    expect_s4_class(result, "MizerParams")
    m_val <- species_params(result)[species_params(result)$species == "Hake", "m"]
    n_val <- species_params(celtic_params)[species_params(celtic_params)$species == "Hake", "n"]
    expect_gt(m_val, n_val)
})

test_that("matchCatch with D_ext free moves it from its starting value", {
    sp <- species_params(celtic_params)
    s <- sp$species == "Hake"
    D_ext_before <- diffusion_coefficient(celtic_params, s, sp$n[s])
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    D_ext_after <- species_params(result)[s, "D_ext"]
    expect_false(isTRUE(all.equal(D_ext_before, D_ext_after)))
    expect_gte(D_ext_after, exp(-20))
})

test_that("matchCatch warns and returns params unchanged when no catch and no production", {
    params_no_prod <- celtic_params
    sp <- species_params(params_no_prod)
    sp[sp$species == "Hake", "production_observed"] <- NA
    species_params(params_no_prod) <- sp

    empty_catch <- celtic_catch[integer(0), ]

    expect_warning(
        result <- matchCatch(params_no_prod, species = "Hake", catch = empty_catch),
        "can not be matched"
    )
    result@time_modified <- params_no_prod@time_modified
    expect_equal(result, params_no_prod)
})

test_that("matchCatch still runs when yield_observed is missing", {
    params_no_yield <- celtic_params
    gp <- gear_params(params_no_yield)
    gp[gp$species == "Hake", "yield_observed"] <- NA
    gear_params(params_no_yield) <- gp

    result <- suppressWarnings(matchCatch(params_no_yield, species = "Hake",
                                          catch = celtic_catch))
    expect_s4_class(result, "MizerParams")
})

test_that("matchCatch still runs when production_observed is missing", {
    params_no_prod <- celtic_params
    sp <- species_params(params_no_prod)
    sp[sp$species == "Hake", "production_observed"] <- NA
    species_params(params_no_prod) <- sp

    result <- suppressWarnings(matchCatch(params_no_prod, species = "Hake",
                                          catch = celtic_catch))
    expect_s4_class(result, "MizerParams")
})

test_that("matchCatch with empty catch matches only yield and production", {
    empty_catch <- celtic_catch[integer(0), ]
    result <- suppressWarnings(matchCatch(celtic_params, species = "Hake",
                                          catch = empty_catch))
    expect_s4_class(result, "MizerParams")
    # Selectivity parameters (l50, l25) must be unchanged since there is no
    # size-distribution data to fit them against.
    gp_before <- gear_params(celtic_params)[gear_params(celtic_params)$species == "Hake",
                                            c("l50", "l25")]
    gp_after  <- gear_params(result)[gear_params(result)$species == "Hake",
                                     c("l50", "l25")]
    expect_equal(gp_before, gp_after)
})

test_that("mu_mat_lim caps the optimised external mortality", {
    tight_lim <- 0.5
    result <- suppressWarnings(matchCatch(celtic_params, species = "Hake",
                                          catch = celtic_catch,
                                          mu_mat_lim = tight_lim))
    mu_mat <- species_params(result)[species_params(result)$species == "Hake", "mu_mat"]
    expect_lte(mu_mat, tight_lim + 1e-6)
})

test_that("map argument fixes the specified parameter at its initial value", {
    l50_before <- gear_params(celtic_params)[gear_params(celtic_params)$species == "Hake", "l50"]
    result <- suppressWarnings(matchCatch(celtic_params, species = "Hake",
                                          catch = celtic_catch,
                                          map = list(l50 = factor(NA))))
    l50_after <- gear_params(result)[gear_params(result)$species == "Hake", "l50"]
    expect_equal(l50_before, l50_after)
})
