# =============================================================================
# Installation Script
# This script checks for and installs any missing R packages required for the
# ctDNA data analysis and visualization in the CAcTUS trial project.
# =============================================================================

# List of required packages
packages <- c(
  "tidyverse", # Includes ggplot2, dplyr, lubridate
  "haven",     # For .dta files
  "survival",  # For survival analysis
  "survminer",  # For survival visualization
  "epiR",      # For epidemiological analysis
  "lubridate", # For date manipulation
  "ggplot2", # For data visualization
  "dplyr", # For data manipulation
  "EnvStats", # For statistical analysis
  "nlme" # For linear mixed-effects models
)

# Function to check, install, and load packages
install_if_missing <- function(p) {
  # 'character.only = TRUE' tells R that 'p' is a variable containing the name
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}

# Run the function over the list
invisible(sapply(packages, install_if_missing))

message("--- All libraries are installed and loaded! ---")