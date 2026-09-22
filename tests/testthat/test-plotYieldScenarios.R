test_that("plotYieldScenarios returns yields relative to the baseline", {
    baseline <- project(NS_params, effort = 1, t_max = 2, t_save = 1,
                        progress_bar = FALSE)
    reduced <- project(NS_params, effort = 0.5, t_max = 2, t_save = 1,
                       progress_bar = FALSE)
    increased <- project(NS_params, effort = 1.5, t_max = 2, t_save = 1,
                         progress_bar = FALSE)

    data <- plotYieldScenarios(
        baseline,
        list(Reduced = reduced, Increased = increased),
        return_data = TRUE
    )

    expect_named(data, c("time", "Species", "yield", "baseline_yield",
                         "Scenario", "relative_yield"))
    expect_equal(levels(data$Scenario), c("Reduced", "Increased"))
    expect_equal(
        data$relative_yield,
        100 * data$yield / data$baseline_yield,
        tolerance = 1e-12
    )
    expect_s3_class(
        plotYieldScenarios(baseline, list(Reduced = reduced)),
        "ggplot"
    )
})

test_that("plotYieldScenarios filters species and gears jointly", {
    baseline <- project(NS_params, effort = 1, t_max = 1, t_save = 1,
                        progress_bar = FALSE)
    scenario <- project(NS_params, effort = 0.5, t_max = 1, t_save = 1,
                        progress_bar = FALSE)

    data <- plotYieldScenarios(baseline, list(Reduced = scenario),
                               species = "Cod", gears = "Otter",
                               return_data = TRUE)

    expected <- as.data.frame.table(getYieldGear(scenario),
                                    responseName = "value")
    expected <- expected[expected$sp == "Cod" & expected$gear == "Otter", ]
    expect_equal(data$yield, expected$value)
})

test_that("plotYieldScenarios facets selected species by default", {
    baseline <- project(NS_params, effort = 1, t_max = 1, t_save = 1,
                        progress_bar = FALSE)
    scenario <- project(NS_params, effort = 0.5, t_max = 1, t_save = 1,
                        progress_bar = FALSE)

    data <- plotYieldScenarios(
        baseline, list(Reduced = scenario),
        species = c("Cod", "Herring"), gears = "Otter", return_data = TRUE
    )
    plot <- plotYieldScenarios(
        baseline, list(Reduced = scenario),
        species = c("Cod", "Herring"), gears = "Otter"
    )

    expect_equal(levels(data$Species), c("Cod", "Herring"))
    expect_s3_class(plot$facet, "FacetWrap")

    combined <- plotYieldScenarios(
        baseline, list(Reduced = scenario),
        species = c("Cod", "Herring"), gears = "Otter",
        facet_species = FALSE, return_data = TRUE
    )
    expect_equal(as.character(combined$Species), rep("Combined", nrow(combined)))
})
