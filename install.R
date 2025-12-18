# List of required packages
packages <- c(
  "tidyverse", # Includes ggplot2, dplyr, lubridate
  "haven",     # For .dta files
  "survival",  # For survival analysis
  "survminer"  # For survival visualization
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