# Cactus Cancer Trial

This repository contains the analyses included in the **CACTUS cancer trial**.

> **Status:** Submitted for publication in *Nature Communications* (2025).  
> **Authors:** Dr. Hitesh Mistry & Dr. Juan Delgado — Systems Forecasting

---

## 📋 Overview
The scripts in this repository facilitate the reproduction of survival analyses and data visualizations presented in the manuscript. The analysis utilizes the following R libraries:
* `tidyverse` & `dplyr` (Data manipulation)
* `haven` (Importing .dta files)
* `lubridate` (Date operations)
* `survival` & `survminer` (Survival analysis and plotting)
* `epiR` (epidemiological analysis)

### Files
The primary files containing most of the analysis are: 
* ctDNA_data_check.R
* primary_secondary_plots_analyses.R

There are a few other side analyses contained in other files:
* time_series_analysis.R
* ctDNA_data_assessment.R
* Cox_PH.R
These analyses might be of use to some users.

### 🚀 Usage
Once the installation is complete, you can run the primary analysis scripts located in the root directory.
Ensure that any required .dta files are placed in the appropriate data directory before running the scripts.

### Usage
To request the data, please refer to Dr Rebecca Lee (rebecca.lee-3@manchester.ac.uk).

---

## 🛠 Installation

To set up the environment on **Ubuntu**, follow these steps:

### 1. Install System Dependencies
R packages on Linux require several system-level libraries to be installed via the terminal first:

```bash
sudo apt update
sudo apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    cmake libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```

### 2. Run the Installation Script

After the system dependencies are installed, run the install.R file to install the necessary R packages:
```bash
Rscript install.R
```

## Contact
For inquiries regarding the forecasting models or data access, please contact the authors at Systems Forecasting (hitesh@systemsforecasting.co.uk and juan@systemsforecasting.co.uk).