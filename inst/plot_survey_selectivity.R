library(dplyr)
library(ggplot2)

required_columns <- c(
    "common.name", "FishLength_cm", "Snfish", "SnfishQ"
)

if (!exists("species_info_GD_length", inherits = TRUE)) {
    stop("`species_info_GD_length` must be loaded before sourcing this script.")
}

missing_columns <- setdiff(required_columns, names(species_info_GD_length))
if (length(missing_columns) > 0) {
    stop(
        "`species_info_GD_length` is missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

survey_selectivity <- species_info_GD_length |>
    transmute(
        species = common.name,
        length_cm = FishLength_cm,
        selectivity = Snfish / SnfishQ
    ) |>
    filter(is.finite(selectivity)) |>
    arrange(species, length_cm)

ggplot(species_data, aes(x = length_cm, y = selectivity)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap()
    labs(
        title = unique(species_data$species),
        x = "Fish length (cm)",
        y = "Survey selectivity (Snfish / SnfishQ)"
    ) +
    theme_minimal()


survey_selectivity_plots <- split(
    survey_selectivity,
    survey_selectivity$species
) |>
    lapply(function(species_data) {
        ggplot(species_data, aes(x = length_cm, y = selectivity)) +
            geom_line(linewidth = 0.7) +
            geom_point(size = 1.5) +
            labs(
                title = unique(species_data$species),
                x = "Fish length (cm)",
                y = "Survey selectivity (Snfish / SnfishQ)"
            ) +
            theme_minimal()
    })

if (interactive()) {
    invisible(lapply(survey_selectivity_plots, print))
}
