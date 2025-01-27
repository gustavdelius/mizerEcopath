params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))
params_single <- newSingleSpeciesParams()
params_single <- setPredKernel(params_single,
                               pred_kernel = getPredKernel(params_single))

test_that("getDiffusion with order=1 can reproduce growth", {
    g <- getEReproAndGrowth(params)
    K <- metab(params)
    d1 <- g + K
    expect_equal(getDiffusion(params, order = 1), d1)
})

test_that("Diffusion follows correct power law for juveniles", {
    # TODO
})
