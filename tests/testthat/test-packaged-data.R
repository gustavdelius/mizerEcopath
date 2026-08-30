test_that("packaged MizerParams data uses the current S3 representation", {
    models <- list(
        celtic_params,
        ns_3_spp_model_initial,
        ns_3_spp_model_trial,
        ns_3_spp_model_final
    )

    for (params in models) {
        expect_s3_class(params, "MizerParams")
        expect_true(is.list(params))
        expect_silent(validParams(params))
    }
})
