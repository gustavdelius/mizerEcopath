# Set the survey gear selectivity from Richard's length-resolved survey data.
#
# `species_info_GD_length` holds, for each species and 1 cm length class, the
# observed number of fish `Snfish` and the catchability-corrected number
# `SnfishQ`. Their ratio is the survey selectivity at length: the fraction of
# the fish present that the survey gear retains. It is already a genuine
# selectivity in [0, 1], so unlike the commercial fishing mortality it needs no
# rescaling and the survey `catchability` is left alone.
#
# The ratio is much less noisy than `Fmean.yr` — it traces a clean curve over
# the well-sampled lengths — but the sparse large-length classes still scatter,
# so it goes through the same smoother.

library(mizer)
library(dplyr)
library(ggplot2)

f_data <- species_info_GD_length
gear <- "survey"

# Names used in the data that differ from the model species names
name_map <- c("Dover sole" = "Sole")

## Smoothing ------------------------------------------------------------------

source(file.path("inst", "smooth_at_length.R"))

# Selectivity is a proportion, so smooth it on the logit scale. That keeps the
# smoothed values inside (0, 1) without any clamping, and makes the imposed
# fall-off below the smallest observed length act on the odds, which for the
# small selectivities down there is the same as acting on the selectivity.
smooth_selectivity <- function(length_cm, sel, n, ...) {
    smooth_at_length(length_cm, sel, n, link = qlogis, inv_link = plogis, ...)
}

## Selectivity on the model's size grid ---------------------------------------

species <- species_params(p)$species
sel <- selectivity(p)

w_grid <- w(p)
sp_a <- species_params(p)$a
sp_b <- species_params(p)$b

sel_curves <- list()  # for plotting

for (i in seq_along(species)) {
    data_name <- names(name_map)[match(species[i], name_map)]
    if (is.na(data_name)) data_name <- species[i]

    d <- f_data |>
        filter(common.name == data_name) |>
        mutate(sel_obs = Snfish / SnfishQ) |>
        filter(is.finite(sel_obs), sel_obs > 0, sel_obs < 1) |>
        arrange(FishLength_cm)
    if (nrow(d) < 5) {
        warning("Not enough survey selectivity data for ", species[i],
                ", skipping.")
        next
    }

    sel_fun <- smooth_selectivity(d$FishLength_cm, d$sel_obs, d$n)

    l_grid <- (w_grid / sp_a[i])^(1 / sp_b[i])
    sel[gear, species[i], ] <- sel_fun(l_grid)

    sel_curves[[species[i]]] <- data.frame(
        species = species[i], length_cm = l_grid, sel = sel_fun(l_grid))
}

# Only the selectivity is replaced; the survey catchability stays at its
# existing tiny value so that the survey gear still imposes no mortality.
p <- setFishing(p, selectivity = sel)

## Check ----------------------------------------------------------------------

smoothed <- bind_rows(sel_curves)
raw <- f_data |>
    mutate(species = ifelse(common.name %in% names(name_map),
                            name_map[common.name], common.name),
           sel_obs = Snfish / SnfishQ) |>
    filter(species %in% smoothed$species)

ggplot(smoothed, aes(length_cm, sel)) +
    geom_point(data = raw, aes(FishLength_cm, sel_obs, size = n),
               inherit.aes = FALSE, alpha = 0.35) +
    geom_line(colour = "blue", linewidth = 0.8) +
    facet_wrap(~species, scales = "free_x") +
    scale_size_area(max_size = 2.5) +
    coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    labs(x = "Length (cm)", y = "Survey selectivity",
         title = "Smoothed survey selectivity (line) over Snfish / SnfishQ (points)")
