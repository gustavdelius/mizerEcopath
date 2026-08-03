params <- celtic_params

test_that("plotReproductiveEfficiency returns the correct data", {
    df <- plotReproductiveEfficiency(params, return_data = TRUE)

    expect_named(df, c("Species", "erepro", "dd_survival", "efficiency"))
    expect_equal(nrow(df), nrow(params@species_params))
    expect_equal(as.character(df$Species), params@species_params$species)
    expect_equal(df$erepro, unname(params@species_params$erepro))
    # The realised efficiency is the product of the two factors
    expect_equal(df$efficiency, df$erepro * df$dd_survival)
    # The survival of density dependence is 1 - reproduction_level
    expect_equal(df$dd_survival, unname(1 - getReproductionLevel(params)))
})

test_that("plotReproductiveEfficiency respects the species argument", {
    df <- plotReproductiveEfficiency(params, species = c("Cod", "Hake"),
                                     return_data = TRUE)
    expect_equal(as.character(df$Species), c("Cod", "Hake"))

    all <- plotReproductiveEfficiency(params, return_data = TRUE)
    expect_equal(df$efficiency,
                 all$efficiency[all$Species %in% c("Cod", "Hake")])

    bogus <- function() plotReproductiveEfficiency(params, species = "Yeti")
    expect_warning(expect_error(bogus(), "No species have been selected"),
                   "species do not exist")
})

test_that("plotReproductiveEfficiency copes with zero reproduction", {
    p <- params
    p@species_params$erepro[[1]] <- 0
    expect_equal(getRDI(p)[[1]], 0)

    df <- plotReproductiveEfficiency(p, return_data = TRUE)
    # 0/0 must not become NaN
    expect_true(is.na(df$dd_survival[[1]]))
    expect_true(is.na(df$efficiency[[1]]))
    # The remaining species are unaffected
    unchanged <- plotReproductiveEfficiency(params, return_data = TRUE)
    expect_equal(df$efficiency[-1], unchanged$efficiency[-1])
    # The species with no reproduction is dropped from the plot rather than
    # breaking the log scale
    expect_s3_class(plotReproductiveEfficiency(p), "ggplot")
})

test_that("plotReproductiveEfficiency returns a plot", {
    p <- plotReproductiveEfficiency(params)
    expect_s3_class(p, "ggplot")
    expect_identical(ggplot2::layer_scales(p)$y$name, "Reproductive efficiency")
})
