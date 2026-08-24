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

test_that("the integrals use bin averaging when second_order_w asks for it", {
    p <- NS_params
    second_order_w(p) <- TRUE
    N <- initialN(p)
    ba <- function(K) {
        n <- ncol(K)
        K[, -n] <- 0.5 * (K[, -n] + K[, -1])
        K
    }

    expect_equal(getSomaticProduction(p),
                 setNames(as.vector((N * ba(getEGrowth(p))) %*% dw(p)),
                          species_params(p)$species))
    expect_equal(getMetabolicRespiration(p),
                 setNames(as.vector((N * ba(metab(p))) %*% dw(p)),
                          species_params(p)$species))
    expect_equal(getConsumption(p),
                 drop((N * ba(getEncounter(p) * (1 - getFeedingLevel(p)))) %*%
                          dw(p)))
    expect_equal(getZB(p),
                 setNames(as.vector((N * ba(sweep(getMort(p), 2, w(p), "*"))) %*%
                                        dw(p)),
                          species_params(p)$species))

    # Bin averaging does change the answer on NS_params' coarse grid
    expect_false(isTRUE(all.equal(getSomaticProduction(p),
                                  getSomaticProduction(NS_params))))
})

test_that("the size range is bin-averaged along with the rest of the weight", {
    p <- NS_params
    second_order_w(p) <- TRUE
    # Splitting the size range in two must give back the total, because the
    # two halves of the straddling bin add up to the whole bin.
    w_split <- w(p)[50]
    expect_equal(getSomaticProduction(p, max_w = w_split) +
                     getSomaticProduction(p, min_w = w_split * 1.0001),
                 getSomaticProduction(p))
})

test_that("getDietMatrix bin-averages along the predator size dimension", {
    p <- NS_params
    second_order_w(p) <- TRUE
    diet <- getDiet(p, proportion = FALSE)
    # Trapezoidal average along the predator size dimension (the 2nd one)
    n <- length(w(p))
    K <- diet
    K[, -n, ] <- 0.5 * (diet[, -n, ] + diet[, -1, ])
    expected <- K |>
        sweep(c(1, 2), initialN(p), "*") |>
        sweep(2, dw(p), "*") |>
        aperm(c(1, 3, 2)) |>
        rowSums(dims = 2)

    expect_equal(getDietMatrix(p), expected)
})

test_that("getDietMatrix agrees with getConsumption", {
    # Only on mizer's default first-order path. With bin averaging switched on,
    # mizer's own getDiet() applies the prey-bin quadrature twice and comes out
    # a factor (1 + beta) / 2 above getEncounter() * (1 - f), so this identity
    # cannot hold there however we discretise the predator-size integral.
    # See https://github.com/sizespectrum/mizer/issues/474.
    expect_equal(rowSums(getDietMatrix(NS_params)), getConsumption(NS_params),
                 tolerance = 1e-3)
})
