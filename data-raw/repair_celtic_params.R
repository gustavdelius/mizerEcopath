# Make `celtic_params` safe against recalculation.
#
# A MizerParams object holds each species parameter twice: the value the user
# supplied, in `given_species_params`, and the value the model actually uses, in
# `species_params`. Whenever a species parameter changes, mizer rebuilds the
# second table from the first, so any parameter that only ever reached
# `species_params` is silently reset to its default at that point.
#
# `celtic_params` was in that state for six parameters, among them
# `interaction_resource` (0, i.e. the species do not feed on the resource) and
# `pred_kernel_type`. A single `species_params(celtic_params) <- ...` therefore
# turned the resource feeding back on and swapped the predation kernel for the
# lognormal default, changing consumption by up to 66%. That is what made
# `matchConsumption()` and `matchCatch()` produce a different model as soon as
# anything downstream touched a species parameter.
#
# The repair only copies values from `species_params` into
# `given_species_params`; it records what the model already used, so the model
# itself is unchanged. Every other slot is verified to be identical below.

devtools::load_all(".")

# Record every species parameter that a recalculation would otherwise change.
# Recording one parameter can change how another is derived, so we iterate to a
# fixed point.
record_calculated_params <- function(params) {
    recorded <- character(0)
    repeat {
        probe <- params
        suppressMessages(species_params(probe) <- species_params(probe))
        changed <- Filter(
            function(col) {
                col %in% names(probe@species_params) &&
                    !isTRUE(all.equal(params@species_params[[col]],
                                      probe@species_params[[col]]))
            },
            names(params@species_params)
        )
        if (length(changed) == 0) break
        for (col in changed) {
            params@given_species_params[[col]] <- params@species_params[[col]]
        }
        recorded <- union(recorded, changed)
    }
    attr(params, "recorded") <- recorded
    params
}

repaired <- record_calculated_params(celtic_params)
message("Recorded: ", toString(attr(repaired, "recorded")))
attr(repaired, "recorded") <- NULL

# The repair must not change the model, only the record of how it was specified.
for (s in setdiff(slotNames(repaired), c("time_modified", "given_species_params"))) {
    stopifnot(isTRUE(all.equal(slot(celtic_params, s), slot(repaired, s))))
}

# And the model must now survive a recalculation untouched. The one thing a
# recalculation is still allowed to do is drop `beta` and `sigma`: with
# `pred_kernel_type = "power_law"` recorded, mizer no longer needs the two
# lognormal-kernel parameters and removes them. (Before the repair the kernel
# type reverted to the lognormal default and those two values were used, which
# is where most of the change in consumption came from.)
check <- repaired
suppressMessages(species_params(check) <- species_params(check))
for (s in setdiff(slotNames(repaired), c("time_modified", "species_params"))) {
    stopifnot(isTRUE(all.equal(slot(repaired, s), slot(check, s))))
}
stopifnot(setdiff(names(repaired@species_params),
                  names(check@species_params)) == c("beta", "sigma"))
for (col in names(check@species_params)) {
    stopifnot(isTRUE(all.equal(repaired@species_params[[col]],
                               check@species_params[[col]])))
}

celtic_params <- repaired
usethis::use_data(celtic_params, overwrite = TRUE, compress = "gzip")
