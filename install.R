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
  "epiR"      # For epidemiological analysis
  "lubridate", # For date manipulation
  "ggplot2", # For data visualization
  "dplyr" # For data manipulation
)

# Function to install missing packages
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

# Run the installation
lapply(packages, install_if_missing)

message("Installation complete!")