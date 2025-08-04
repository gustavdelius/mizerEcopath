#load in libraries
library(here)
library(dplyr)
library(Hmisc)

# Load age at size data compiled by Jess ----
survey <- readRDS(here::here("inst", "extdata",
                             "Celtic_Sea_Size_at_Age_Data.rds"))
# Aggregate to cm bins
Y_obs <- survey |>
    # remove rows with NA in any column
    filter(!is.na(LngtClass) & !is.na(Age) & !is.na(CANoAtLngt)) |>
    # round down to cm
    mutate(LngtClass = floor(LngtClass)) |>
    # aggregate counts
    group_by(Scientific_name, LngtClass, Age) |>
    summarise(CANoAtLngt = sum(CANoAtLngt, na.rm = TRUE), .groups = "drop")
#add columns ensuring values are integers
Y_obs$lenBin <- as.integer(Y_obs$LngtClass)
Y_obs$ageClass<-as.integer(Y_obs$Age)
Y_obs$count<-as.integer(Y_obs$CANoAtLngt)

#determine quartiles and median for each length and species, weight by counts
Quartile_data <- Y_obs %>%
    group_by(lenBin, Scientific_name) %>%
    summarise(
        Q1 = wtd.quantile(ageClass, weights = count, probs = 0.25, na.rm = TRUE),
        Median = wtd.quantile(ageClass, weights = count, probs = 0.5, na.rm = TRUE),
        Q3 = wtd.quantile(ageClass, weights = count, probs = 0.75, na.rm = TRUE),
        .groups = 'drop'
    )

#plot the data
ggplot(Quartile_data, aes(x = lenBin, y = Median, color = Scientific_name)) +
    geom_line() +
    geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = Scientific_name), alpha = 0.2, color = NA) +
    labs(x = "Length [cm]", y = "Age (median and IQR)", color = "Scientific_name", fill = "Scientific_name") +
    theme_minimal()
