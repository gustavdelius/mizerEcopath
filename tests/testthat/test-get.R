test_that("getSomaticProduction integrates g N over all sizes by default", {
    expected <- as.vector((getEGrowth(NS_params) * initialN(NS_params)) %*%
                               dw(NS_params))
    names(expected) <- species_params(NS_params)$species

    expect_equal(getSomaticProduction(NS_params), expected)
})

test_that("getSomaticProduction can be restricted to a size range", {
    min_w <- 10
    max_w <- 1000
    sel <- NS_params@w >= min_w & NS_params@w <= max_w
    expected <- as.vector((getEGrowth(NS_params)[, sel, drop = FALSE] *
                                initialN(NS_params)[, sel, drop = FALSE]) %*%
                               dw(NS_params)[sel])
    names(expected) <- species_params(NS_params)$species

    expect_equal(getSomaticProduction(NS_params, min_w = min_w, max_w = max_w),
                 expected)
})
