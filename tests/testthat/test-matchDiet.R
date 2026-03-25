test_that("matchDiet throws error for non-MizerParams input", {
    expect_error(matchDiet("not_params"),
                 "params must be a MizerParams object.")
})

test_that("matchDiet(params, diet_matrix = getDietMatrix(params)) leaves params unchanged", {
    dm <- getDietMatrix(celtic_params)
    result <- matchDiet(celtic_params, diet_matrix = dm)
    result@time_modified <- celtic_params@time_modified
    expect_equal(result, celtic_params)
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
})

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
