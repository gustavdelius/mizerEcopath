context("getDiffusion uses n argument")

params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))


n_doubled <- initialN(params) * 2

test_that("getDiffusion depends on n", {
  d_default <- getDiffusion(params)
  d_modified <- getDiffusion(params, n = n_doubled)
  expect_false(identical(d_default, d_modified))
})
