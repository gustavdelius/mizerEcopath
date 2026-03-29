# test-matchCatch.R
library(testthat)
library(mizer)

# Helper: build celtic_params with two gears for Hake
# (sigmoid_length + double_sigmoid_length) and corresponding two-gear catch data.
make_two_gear_hake <- function() {
    params <- celtic_params
    gp <- gear_params(params)
    gp$l50_right <- NA_real_
    gp$l25_right <- NA_real_

    hake_row <- gp[gp$species == "Hake", ]

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

    hake_catch <- celtic_catch[celtic_catch$species == "Hake", ]
    small_catch <- hake_catch
    small_catch$gear <- "small_mesh"
    small_catch$catch <- round(hake_catch$catch * 0.4)
    large_catch <- hake_catch
    large_catch$gear <- "large_mesh"
    large_catch$catch <- round(hake_catch$catch * 0.6)

    list(params = params, catch = rbind(small_catch, large_catch))
}

test_that("matchCatch runs without error with valid inputs for a single species", {
    # Run the function for Hake
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    # Running it again should not make further changes except on time-stamp
    result2 <- matchCatch(result, species = "Hake", catch = celtic_catch)
    result2@time_modified <- result@time_modified
    expect_equal(result, result2)
})

test_that("matchCatch throws error if catch data is missing required columns", {
    # Remove 'dl' column from catch data
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
        # Expect that changes have been applied to both species
        expect_true(all(multi_species %in% gp$species))
})

test_that("matchCatch preserves biomass for the adjusted species", {
    # The documentation suggests biomass remains unchanged for the adjusted species
    biomass_before <- getBiomass(celtic_params, use_cutoff = TRUE)
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    biomass_after <- getBiomass(result, use_cutoff = TRUE)

    # Check that biomass difference is small for Hake
    expect_equal(biomass_before["Hake"], biomass_after["Hake"], tolerance = 1e-3)
})

test_that("matchCatch updates gear selectivity parameters for the selected species", {
    gp_before <- gear_params(celtic_params)

    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch) |>
        suppressWarnings()
    gp_after <- gear_params(result)

    sp_idx <- which(gp_after$species == "Hake")
    # Expect that some key gear parameters have changed
    expect_false(all(gp_before[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_after[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch respects the yield_lambda parameter", {
    # Run with default yield_lambda
    result_default <- matchCatch(celtic_params, species = "Hake",
                                 catch = celtic_catch) |>
        suppressWarnings()
    gp_default <- gear_params(result_default)

    # Run with altered yield_lambda
    result_mod <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch, yield_lambda = 10)
    gp_mod <- gear_params(result_mod)

    sp_idx <- which(gp_mod$species == "Hake")
    # Different yield_lambda should lead to different parameter estimates
    expect_false(all(gp_default[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_mod[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch handles invalid species inputs gracefully", {
    # Species not in params
    expect_warning(
        matchCatch(celtic_params, species = c("Cod", "NotASpecies"),
                   catch = celtic_catch),
        "The following species do not exist: NotASpecies."
    )

    # Empty species vector
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

    # Both gears get updated selectivity parameters
    expect_false(all(sigmoid_row[, c("l50", "l25", "catchability")] ==
                     gear_params(tg$params)[gear_params(tg$params)$gear == "small_mesh",
                                            c("l50", "l25", "catchability")]))
    # double_sigmoid_length gear has non-NA right-side parameters
    expect_true(!is.na(double_row$l50_right))
    expect_true(!is.na(double_row$l25_right))
    # sigmoid gear leaves right-side parameters as NA
    expect_true(is.na(sigmoid_row$l50_right))
})

test_that("matchCatch with m free produces a valid result", {
    result <- suppressWarnings(matchCatch(celtic_params, species = "Hake",
                                          catch = celtic_catch,
                                          map = list(m = NULL)))
    expect_s4_class(result, "MizerParams")
    m_val <- species_params(result)[species_params(result)$species == "Hake", "m"]
    n_val <- species_params(celtic_params)[species_params(celtic_params)$species == "Hake", "n"]
    # m must be > n (lower bound in the optimiser is n * 1.01)
    expect_gt(m_val, n_val)
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
    # Params returned unchanged (modulo timestamp)
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
    # size-distribution data to fit them against
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
