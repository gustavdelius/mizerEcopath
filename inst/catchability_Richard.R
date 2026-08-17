# Set the commercial gear selectivity from Richard's length-resolved fishing
# mortality estimates.
#
# `species_info_GD_length` holds, for each species and 1 cm length class, a mean
# fishing mortality `Fmean.yr` together with the number of fish `n` on which it
# is based. The raw F-at-length is very noisy, especially in the largest length
# classes where `n` is small, so we smooth it before using it.
#
# Because mizer decomposes fishing mortality as
#     F_{g,i}(w) = S_{g,i}(w) * Q_{g,i} * E_g,
# a smoothed F-at-length curve is imposed by writing the normalised curve into
# `selectivity(p)` and its maximum into `catchability(p)` (effort stays 1).

library(mizer)
library(dplyr)
library(ggplot2)

f_data <- species_info_GD_length
gear <- "commercial"

# Names used in the data that differ from the model species names
name_map <- c("Dover sole" = "Sole")

## Smoothing ------------------------------------------------------------------

source(file.path("inst", "smooth_at_length.R"))

# F is a positive rate, so smooth it on the log scale.
smooth_F <- function(length_cm, F_yr, n, ...) {
    smooth_at_length(length_cm, F_yr, n, link = log, inv_link = exp, ...)
}

## Fishing mortality on the model's size grid ---------------------------------

species <- species_params(p)$species
sel <- selectivity(p)
catcha <- catchability(p)

# Lengths (cm) corresponding to the weight grid, per species
w_grid <- w(p)
sp_a <- species_params(p)$a
sp_b <- species_params(p)$b

F_curves <- list()  # for plotting

for (i in seq_along(species)) {
    data_name <- names(name_map)[match(species[i], name_map)]
    if (is.na(data_name)) data_name <- species[i]

    d <- f_data |>
        filter(common.name == data_name, Fmean.yr > 0) |>
        arrange(FishLength_cm)
    if (nrow(d) < 5) {
        warning("Not enough F-at-length data for ", species[i], ", skipping.")
        next
    }

    F_fun <- smooth_F(d$FishLength_cm, d$Fmean.yr, d$n)

    l_grid <- (w_grid / sp_a[i])^(1 / sp_b[i])
    F_grid <- F_fun(l_grid)

    # Split into a catchability (the maximum F) and a selectivity in [0, 1]
    F_max <- max(F_grid)
    sel[gear, species[i], ] <- F_grid / F_max
    catcha[gear, species[i]] <- F_max

    F_curves[[species[i]]] <- data.frame(
        species = species[i], length_cm = l_grid, F_yr = F_grid)
}

p <- setFishing(p, selectivity = sel, catchability = catcha)

## Check ----------------------------------------------------------------------

smoothed <- bind_rows(F_curves)
raw <- f_data |>
    mutate(species = ifelse(common.name %in% names(name_map),
                            name_map[common.name], common.name)) |>
    filter(species %in% smoothed$species)

ggplot(smoothed, aes(length_cm, F_yr)) +
    geom_point(data = raw, aes(FishLength_cm, Fmean.yr, size = n),
               inherit.aes = FALSE, alpha = 0.4) +
    geom_line(colour = "blue", linewidth = 0.8) +
    facet_wrap(~species, scales = "free") +
    # The imposed fall-off below the data spans many decades; clip the axis so
    # that the fitted region stays readable.
    scale_y_log10() +
    coord_cartesian(ylim = c(min(raw$Fmean.yr) / 10, NA)) +
    scale_size_area(max_size = 3) +
    labs(x = "Length (cm)", y = "Fishing mortality (1/year)",
         title = "Smoothed F at length (line) over the raw estimates (points)")
