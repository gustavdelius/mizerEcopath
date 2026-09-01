# Upgrade the MizerParams objects shipped in data/ to mizer's S3 representation.
#
# `validParams()` performs the supported conversion from models saved by older
# mizer versions and reconciles their species-parameter provenance. Package data
# must use `usethis::use_data()` rather than `saveParams()` so that the objects
# remain available through R's lazy-data mechanism.

devtools::load_all(".")

celtic_params <- validParams(celtic_params, info_level = 0)
ns_3_spp_model_initial <- validParams(ns_3_spp_model_initial, info_level = 0)
ns_3_spp_model_trial <- validParams(ns_3_spp_model_trial, info_level = 0)
ns_3_spp_model_final <- validParams(ns_3_spp_model_final, info_level = 0)

stopifnot(
    is.list(celtic_params),
    is.list(ns_3_spp_model_initial),
    is.list(ns_3_spp_model_trial),
    is.list(ns_3_spp_model_final)
)

usethis::use_data(celtic_params, overwrite = TRUE, compress = "gzip")
usethis::use_data(ns_3_spp_model_initial, overwrite = TRUE, compress = "gzip")
usethis::use_data(ns_3_spp_model_trial, overwrite = TRUE, compress = "gzip")
usethis::use_data(ns_3_spp_model_final, overwrite = TRUE, compress = "gzip")
