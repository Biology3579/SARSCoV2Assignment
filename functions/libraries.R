## Script name: libraries.R
##
## Purpose of script: 
##    A function to load all necessary libraries for the SARS-CoV-2 analysis, to declutter main analysis file.
##    An outline of the use of each library is also provided
##
## Author: Biology3579
##
## Date Created: 2025-03-18 
##
## ---------------------------


# Function to load all required libraries for the SARS-CoV-2 analysis
load_libraries <- function() {
  
  # core data wrangling and processing
  library(tidyverse)     # includes dplyr, tidyr, ggplot2, readr for data manipulation and visualization
  library(dplyr)         # data wrangling (mutate, filter, group_by, summarise, rename, arrange)
  library(tidyr)         # data tidying (pivoting, handling missing values, reshaping data)
  library(lubridate)     # date handling (as.Date, floor_date, time-series processing)
  library(here)          # file path management (used for sourcing functions and saving data)
  library(readr)         # reading and writing csv files (for loading datasets)
  library(janitor)       # cleaning column names, handling missing values (used for data cleaning)
  
  # data visualization
  library(ggplot2)       # core plotting package (used for all visualizations)
  library(viridis)       # colorblind-friendly color palettes (for consistent variant and region coloring)
  library(scales)        # formatting axes (date labels, continuous scale adjustments)
  library(patchwork)     # arranging multiple plots together (for combining visualizations)
  library(ggpubr)        # enhanced plotting utilities (e.g., adding p-values, statistical summaries)
  library(svglite)       # exporting plots as svg (for high-quality vector graphics)
  
  # statistical modeling and epidemic analysis
  library(nls.multstart) # logistic growth model fitting (used for variant growth estimation)
  library(nlstools)      # non-linear regression tools (for logistic model diagnostics)
  library(EpiEstim)      # estimating time-varying reproduction number rt
  library(incidence)     # incidence curve analysis (used in epidemiological practicals)
  
  # time-series analysis and rolling averages
  library(zoo)           # rolling averages and time-series transformations (used in smoothing data)
  library(slider)        # sliding window calculations (alternative to zoo for rolling averages)
  
  # report formatting and reproducibility
  library(knitr)         # dynamic reporting (used for formatting reports in rmarkdown/quarto)
  library(kableExtra)    # enhancing tables in reports (used in markdown-based documents)
  library(tinytex)       # latex-related tasks (for compiling pdf reports)
  
  # environment and package management
  library(renv)         # managing package environments for reproducibility
}


