test_that("getSomaticProduction integrates flux over all sizes by default", {
    flux <- sweep(getFlux(NS_params), 2, dw(NS_params), "*")
    expected <- rowSums(flux)
    names(expected) <- species_params(NS_params)$species

    expect_equal(getSomaticProduction(NS_params), expected)
})

test_that("getSomaticProduction can be restricted to a size range", {
    min_w <- 10
    max_w <- 1000
    sel <- NS_params@w >= min_w & NS_params@w <= max_w
    flux <- sweep(getFlux(NS_params), 2, dw(NS_params), "*")
    expected <- rowSums(flux[, sel, drop = FALSE])
    names(expected) <- species_params(NS_params)$species

    expect_equal(getSomaticProduction(NS_params, min_w = min_w, max_w = max_w),
                 expected)
})
