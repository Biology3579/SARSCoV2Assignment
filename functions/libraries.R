## Script name: libraries.R
##
## Purpose of script: 
##    A function to load all necessary libraries for the SARS-CoV-2 analysis, to declutter main analysis file.
##    An outline of the use of each library is also provided.
##    A summary of key libraries is also provided in each functions script (but here they are all combined).
##
## Author: Biology3579
##
## Date Created: 2025-03-18 
##
## ---------------------------

# IMPORTANT: The packages should be installed - with the correct version history -
# using the renv::restore function in the main script. 
# Run that line (in the main script) before using this command to ensure that the correct package version is installed!

# Function to load the required libraries for the SARS-CoV-2 analysis ----
load_libraries <- function() {
  
  # Libraries for core data wrangling and processing ----
  library(tidyverse)     # includes dplyr, tidyr, ggplot2, readr for data manipulation and plotting
  library(dplyr)         # data wrangling (mutate, filter, group_by, summarise, rename, arrange)
  library(tidyr)         # data tidying (pivoting, handling missing values, reshaping data)
  library(lubridate)     # date handling (as.Date, floor_date, time-series processing)
  library(here)          # file path management (used for sourcing functions and saving data)
  library(readr)         # reading and writing csv and rds files (for loading datasets)
  library(janitor)       # cleaning column names, handling missing values (used for data cleaning)
  
  # Libraries for data visualization ----
  library(ggplot2)       # core plotting package (used for all plots)
  library(viridis)       # colorblind-friendly colour palettes (used as a reference for consistent variant and region colouring)
  library(scales)        # formatting axes (date labels, continuous scale adjustments)
  library(patchwork)     # arranging multiple plots together (for combining plots)
  library(ggpubr)        # enhanced plotting utilities (e.g., adding p-values, statistical summaries)
  
  # Libraries for statistical modeling and epidemic analysis ----
  library(nls.multstart) # logistic growth model fitting (used for variant growth estimation)
  library(nlstools)      # non-linear regression tools (for logistic model diagnostics)
  library(EpiEstim)      # estimating time-varying reproduction number rt
  library(incidence)     # incidence curve analysis 
  
  # Libraries for time-series analysis and rolling averages ----
  library(zoo)           # rolling averages and time-series transformations (used in smoothing data)
  
  # Libraries for report formatting and reproducibility ----
  library(knitr)         # dynamic reporting (used for formatting reports in rmarkdown/quarto)
  library(tinytex)       # latex-related tasks (for compiling pdf reports)
  
  # Libraries for environment and package management ----
  library(renv)         # managing package environments for reproducibility
}


