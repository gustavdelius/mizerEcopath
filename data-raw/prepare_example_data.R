# ================================================================
# Purpose: Convert curated example .rds files into .rda datasets
#          for use in the mizerEcopath package.
# Author: James Rimmer
# ================================================================

# Load packages
library(usethis)

# Load objects from data-raw
species_params_example <- readRDS("data-raw/example_species_params.rds")
ecopath_diet_example <- readRDS("data-raw/ecopath_diet_full.rds")
diet_matrix_example <- readRDS("data-raw/example_diet_matrix.rds")

# Save to data/ for internal use by the package
usethis::use_data(
    species_params_example,
    ecopath_diet_example,
    diet_matrix_example,
    overwrite = TRUE
)
