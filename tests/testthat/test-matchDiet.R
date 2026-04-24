test_that("matchDiet throws error for non-MizerParams input", {
    expect_error(matchDiet("not_params"),
                 "params must be a MizerParams object.")
})

test_that("matchDiet(params, diet_matrix = getDietMatrix(params)) leaves params unchanged", {
    dm <- getDietMatrix(celtic_params)
    result <- matchDiet(celtic_params, diet_matrix = dm)
    result@time_modified <- celtic_params@time_modified
    expect_equal(result, celtic_params)

    # Same with North Sea model
    params <- NS_params
    params@ext_encounter[] <- getEncounter(NS_params)
    dm <- getDietMatrix(params)
    result <- matchDiet(params, diet_matrix = dm)
    result@time_modified <- params@time_modified
    expect_equal(result, params)
})

test_that("matchDiet achieves the target diet matrix", {
    sp <- celtic_params@species_params$species
    dm <- getDietMatrix(celtic_params)

    # Reduce predation on species 2 by 20%, leaving other columns unchanged.
    # matchDiet normalises rows internally, so this shifts the proportion of
    # species 2 in each predator's diet downward.
    dm_target <- dm
    dm_target[, sp[2]] <- dm_target[, sp[2]] * 0.8

    result <- matchDiet(celtic_params, diet_matrix = dm_target)

    # Expected absolute flows: normalise rows of dm_target then scale by Q
    Q <- getConsumption(celtic_params)
    expected <- sweep(
        dm_target[sp, sp] / rowSums(dm_target[sp, , drop = FALSE]),
        1, Q, "*"
    )
    expect_equal(getDietMatrix(result)[sp, sp], expected, tolerance = 1e-6)

    # Same as above but with North Sea model
    params <- NS_params
    sp <- params@species_params$species
    params@ext_encounter[] <- getEncounter(NS_params)
    dm <- getDietMatrix(params)
    dm_target <- dm
    dm_target[, sp[2]] <- dm_target[, sp[2]] * 0.9
    result <- matchDiet(params, diet_matrix = dm_target)
    # Expected absolute flows: normalise rows of dm_target then scale by Q
    Q <- getConsumption(params)
    expected <- sweep(
        dm_target[sp, sp] / rowSums(dm_target[sp, , drop = FALSE]),
        1, Q, "*"
    )
    expect_equal(getDietMatrix(result)[sp, sp], expected, tolerance = 1e-6)
})

# checkDietMatrix errors ---------------------------------------------------

test_that("checkDietMatrix errors on non-matrix input", {
    dm <- getDietMatrix(celtic_params)
    expect_error(matchDiet(celtic_params, diet_matrix = as.data.frame(dm)),
                 "`diet_matrix` must be a matrix.")
})

test_that("checkDietMatrix errors on non-numeric matrix", {
    dm <- getDietMatrix(celtic_params)
    dm_char <- matrix(as.character(dm), nrow = nrow(dm), ncol = ncol(dm),
                      dimnames = dimnames(dm))
    expect_error(matchDiet(celtic_params, diet_matrix = dm_char),
                 "`diet_matrix` must be numeric.")
})

test_that("checkDietMatrix errors on NaN entries", {
    dm <- getDietMatrix(celtic_params)
    dm[1, 1] <- NaN
    expect_error(matchDiet(celtic_params, diet_matrix = dm),
                 "`diet_matrix` contains NaNs.")
})

test_that("checkDietMatrix errors on NA entries", {
    dm <- getDietMatrix(celtic_params)
    dm[1, 1] <- NA
    expect_error(matchDiet(celtic_params, diet_matrix = dm),
                 "`diet_matrix` contains NAs.")
})

test_that("checkDietMatrix errors on negative entries", {
    dm <- getDietMatrix(celtic_params)
    dm[1, 1] <- -0.1
    expect_error(matchDiet(celtic_params, diet_matrix = dm),
                 "`diet_matrix` contains negative values.")
})

test_that("checkDietMatrix errors when a predator row sums to zero", {
    dm <- getDietMatrix(celtic_params)
    dm[1, ] <- 0
    expect_error(matchDiet(celtic_params, diet_matrix = dm),
                 "some species do not eat anything.")
})

# convertDietMatrix errors -------------------------------------------------

test_that("convertDietMatrix errors when model species are missing from rows", {
    dm <- getDietMatrix(celtic_params)
    sp <- celtic_params@species_params$species
    # Drop the first species row
    expect_error(matchDiet(celtic_params, diet_matrix = dm[-1, , drop = FALSE]),
                 "diet_matrix does not include all model species as rows.")
})

# Edge cases ---------------------------------------------------------------

test_that("matchDiet works with a single species", {
    sp <- celtic_params@species_params$species
    params1 <- removeSpecies(celtic_params, sp[-1])
    dm <- getDietMatrix(params1)
    result <- matchDiet(params1, diet_matrix = dm)
    result@time_modified <- params1@time_modified
    expect_equal(result, params1)
})

test_that("matchDiet handles a predator whose diet is entirely non-species", {
    dm <- getDietMatrix(celtic_params)
    sp <- celtic_params@species_params$species
    # Zero out all species-species entries for the first predator so it eats
    # only "other" prey; theta should become 0 (not NaN) for that row.
    dm[sp[1], sp] <- 0
    result <- matchDiet(celtic_params, diet_matrix = dm)
    expect_equal(interaction_matrix(result)[sp[1], ], setNames(rep(0, length(sp)), sp))
})

# matchDiet aggregates non-species prey columns into other -----------------

test_that("matchDiet aggregates non-species prey columns into other", {
    dm <- getDietMatrix(celtic_params)
    sp <- celtic_params@species_params$species

    # Split the first non-species column into two artificial prey columns;
    # the result should be identical to passing the original diet matrix.
    non_sp_cols <- setdiff(colnames(dm), sp)
    if (length(non_sp_cols) > 0) {
        col <- non_sp_cols[[1]]
        dm_split <- cbind(
            dm[, sp, drop = FALSE],
            exotic_prey_A = dm[, col] * 0.6,
            exotic_prey_B = dm[, col] * 0.4,
            dm[, non_sp_cols[-1], drop = FALSE]
        )
        result_original <- matchDiet(celtic_params, diet_matrix = dm)
        result_split    <- matchDiet(celtic_params, diet_matrix = dm_split)
        result_split@time_modified <- result_original@time_modified
        expect_equal(result_original, result_split)
    } else {
        skip("No non-species columns in getDietMatrix output to test aggregation")
    }
})

test_that("matchDiet normalises rows so scale of input does not matter", {
    dm <- getDietMatrix(celtic_params)
    # Scale all rows by arbitrary species-specific constants; the result
    # should be identical because rows are normalised internally.
    scale_factors <- seq(0.5, 2, length.out = nrow(dm))
    dm_scaled <- dm * scale_factors  # recycles row-wise (column-major order)
    result_original <- matchDiet(celtic_params, diet_matrix = dm)
    result_scaled   <- matchDiet(celtic_params, diet_matrix = dm_scaled)
    result_scaled@time_modified <- result_original@time_modified
    expect_equal(result_original, result_scaled)
})
