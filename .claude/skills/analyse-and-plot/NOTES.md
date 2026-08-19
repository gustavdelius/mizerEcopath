# Project notes

- `inst/species_info_GD_length.Rdata` contains 571 survey-at-length records for
  11 species. The plotting fields are `common.name`, `FishLength_cm`, `Snfish`,
  and `SnfishQ`; the current data have no zero values in `Snfish`.
- `inst/plot_survey_selectivity.R` calculates the requested survey-selectivity
  measure as `SnfishQ / Snfish` and creates a named ggplot for each species.
- `Fmean.yr` provides fishing mortality at length for every record in that data
  set. `inst/plot_fishing_mortality.R` creates a named ggplot for each species.
- The final allometric model is saved as
  `inst/params_allometric_final_Richard_and_Jess.rds`. It contains nine of the
  11 observed species: model species `Sole` maps to observed `Dover sole`, while
  Angler and Red gurnard have no model counterpart.
