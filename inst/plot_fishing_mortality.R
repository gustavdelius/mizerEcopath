library(dplyr)
library(ggplot2)
library(mizer)

required_columns <- c("common.name", "FishLength_cm", "Fmean.yr")

if (!exists("species_info_GD_length", inherits = TRUE)) {
    stop("`species_info_GD_length` must be loaded before sourcing this script.")
}

if (!exists("params_allometric_final_Richard_and_Jess", inherits = TRUE)) {
    stop(
        "`params_allometric_final_Richard_and_Jess` must be loaded before ",
        "sourcing this script."
    )
}

missing_columns <- setdiff(required_columns, names(species_info_GD_length))
if (length(missing_columns) > 0) {
    stop(
        "`species_info_GD_length` is missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

fishing_mortality <- species_info_GD_length |>
    transmute(
        species = common.name,
        length_cm = FishLength_cm,
        mortality = Fmean.yr
    ) |>
    filter(is.finite(mortality)) |>
    arrange(species, length_cm)

observed_length_ranges <- fishing_mortality |>
    group_by(species) |>
    summarise(
        min_length = min(length_cm),
        max_length = max(length_cm),
        .groups = "drop"
    )

model_length_parameters <- species_params(
    params_allometric_final_Richard_and_Jess
) |>
    as.data.frame() |>
    select(species, a, b)

model_fishing_mortality <- getFMort(
    params_allometric_final_Richard_and_Jess
) |>
    as.data.frame() |>
    transmute(
        model_species = Species,
        w,
        mortality = value
    ) |>
    inner_join(
        model_length_parameters,
        by = c("model_species" = "species")
    ) |>
    transmute(
        species = if_else(
            model_species == "Sole",
            "Dover sole",
            model_species
        ),
        length_cm = (w / a)^(1 / b),
        mortality
    ) |>
    inner_join(observed_length_ranges, by = "species") |>
    filter(dplyr::between(length_cm, min_length, max_length)) |>
    select(species, length_cm, mortality) |>
    arrange(species, length_cm)

fishing_mortality_plots <- split(
    fishing_mortality,
    fishing_mortality$species
) |>
    lapply(function(species_data) {
        current_species <- species_data$species[[1]]
        model_data <- model_fishing_mortality |>
            filter(species == current_species)

        ggplot() +
            geom_line(
                data = species_data,
                aes(x = length_cm, y = mortality, colour = "Observed"),
                linewidth = 0.7
            ) +
            geom_point(
                data = species_data,
                aes(x = length_cm, y = mortality, colour = "Observed"),
                size = 1.5
            ) +
            geom_line(
                data = model_data,
                aes(x = length_cm, y = mortality, colour = "Model"),
                linewidth = 1
            ) +
            scale_colour_manual(
                values = c(Observed = "black", Model = "#D55E00"),
                name = NULL
            ) +
            labs(
                title = current_species,
                x = "Fish length (cm)",
                y = expression("Fishing mortality (" * year^-1 * ")")
            ) +
            theme_minimal()
    })

if (interactive()) {
    invisible(lapply(fishing_mortality_plots, print))
}
