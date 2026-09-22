#' Plot scenario yields relative to a baseline
#'
#' Compares selected scenario yields with a baseline simulation at every saved
#' time point. A value of 100 means that a scenario has the same yield as the
#' baseline, while 50 and 200 mean half and twice the baseline yield,
#' respectively.
#'
#' @param baseline A `mizer::MizerSim` object defining the baseline scenario.
#' @param scenarios A non-empty named list of `mizer::MizerSim` objects. The
#'   simulations must have the same saved times as `baseline`.
#' @param species Optional character vector of species to include. By default,
#'   all species are included.
#' @param gears Optional character vector of fishing gears to include. By
#'   default, all gears are included.
#' @param facet_species A boolean indicating whether to facet the plot by
#'   species. Defaults to `TRUE`; selected gears are summed within each species.
#'   Set to `FALSE` to sum the selected species and gears into a single yield.
#' @param return_data A boolean value that determines whether the data frame
#'   underlying the plot is returned instead of the plot itself. Defaults to
#'   `FALSE`.
#'
#' @return A ggplot2 object, or a data frame if `return_data = TRUE`. The data
#'   frame has one row per scenario and saved time point, with columns `time`,
#'   `Species`, `Scenario`, `yield`, `baseline_yield`, and `relative_yield`.
#' @export
#'
#' @examples
#' \dontrun{
#' baseline <- project(NS_params, t_max = 20, effort = 1)
#' reduced <- project(NS_params, t_max = 20, effort = 0.5)
#' increased <- project(NS_params, t_max = 20, effort = 1.5)
#' plotYieldScenarios(baseline,
#'                    list("Reduced effort" = reduced,
#'                         "Increased effort" = increased),
#'                    species = c("Cod", "Herring"), gears = "Otter")
#' }
plotYieldScenarios <- function(baseline, scenarios, species = NULL,
                               gears = NULL, facet_species = TRUE,
                               return_data = FALSE) {
    if (!methods::is(baseline, "MizerSim")) {
        stop("`baseline` must be a MizerSim object.")
    }
    if (!is.list(scenarios) || length(scenarios) == 0L) {
        stop("`scenarios` must be a non-empty named list of MizerSim objects.")
    }
    scenario_names <- names(scenarios)
    if (is.null(scenario_names) || anyNA(scenario_names) ||
            any(scenario_names == "") || anyDuplicated(scenario_names)) {
        stop("`scenarios` must have unique, non-empty names.")
    }
    if (!all(vapply(scenarios, methods::is, logical(1), class2 = "MizerSim"))) {
        stop("Every element of `scenarios` must be a MizerSim object.")
    }
    if (!is.logical(facet_species) || length(facet_species) != 1L ||
            is.na(facet_species)) {
        stop("`facet_species` must be TRUE or FALSE.")
    }

    yield_data <- function(sim) {
        data <- as.data.frame.table(getYieldGear(sim), responseName = "value")
        data$time <- as.numeric(as.character(data$time))
        if (!is.null(species)) {
            missing_species <- setdiff(species, unique(data$sp))
            if (length(missing_species) > 0L) {
                stop("The following species are not in every simulation: ",
                     paste(missing_species, collapse = ", "), ".")
            }
            data <- data[data$sp %in% species, , drop = FALSE]
        }
        if (!is.null(gears)) {
            missing_gears <- setdiff(gears, unique(data$gear))
            if (length(missing_gears) > 0L) {
                stop("The following gears are not in every simulation: ",
                     paste(missing_gears, collapse = ", "), ".")
            }
            data <- data[data$gear %in% gears, , drop = FALSE]
        }
        if (facet_species) {
            data <- stats::aggregate(value ~ time + sp, data = data, FUN = sum)
            names(data)[names(data) == "sp"] <- "Species"
        } else {
            data <- stats::aggregate(value ~ time, data = data, FUN = sum)
            data$Species <- "Combined"
        }
        data
    }

    baseline_yield <- yield_data(baseline)
    names(baseline_yield)[names(baseline_yield) == "value"] <- "baseline_yield"
    result <- lapply(seq_along(scenarios), function(i) {
        scenario_yield <- yield_data(scenarios[[i]])
        names(scenario_yield)[names(scenario_yield) == "value"] <- "yield"
        data <- merge(scenario_yield, baseline_yield, by = c("time", "Species"),
                      all = FALSE, sort = TRUE)
        if (nrow(data) != nrow(baseline_yield) ||
                nrow(data) != nrow(scenario_yield)) {
            stop("Each scenario must have the same saved times as `baseline`.")
        }
        data$Scenario <- scenario_names[[i]]
        data$relative_yield <- 100 * data$yield / data$baseline_yield
        data$relative_yield[data$baseline_yield == 0] <- NA_real_
        data
    })
    result <- do.call(rbind, result)
    rownames(result) <- NULL
    result$Species <- factor(result$Species, levels = unique(result$Species))
    result$Scenario <- factor(result$Scenario, levels = scenario_names)

    if (return_data) return(result)

    plot <- ggplot(result, aes(x = .data$time, y = .data$relative_yield,
                               colour = .data$Scenario,
                               linetype = .data$Scenario)) +
        geom_hline(yintercept = 100, colour = "grey50", linetype = "dashed") +
        geom_line() +
        labs(x = "Time [years]", y = "Yield relative to baseline [%]",
             colour = "Scenario", linetype = "Scenario")
    if (facet_species) plot <- plot + facet_wrap(vars(.data$Species))
    plot
}
