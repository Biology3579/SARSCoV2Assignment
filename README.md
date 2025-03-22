# SARS-CoV-2 Variant Analysis
This repository contains the scripts and files for analyzing the evolution and spread of SARS-CoV-2 variants. The project focuses on quantifying variant fitness, growth dynamics, and epidemiological trends using genomic and epidemiological data sources.

All the necessary code for the analysis can be found in the [SARSCoV2Analysis.Rmd]([https://github.com/Biology3579/ReproducibleScienceAssignment/blob/main/PenguinAnalysis.Rmd](https://github.com/Biology3579/SARSCoV2Assignment/blob/main/SARSCoV2Analysis.Rmd)) file. This R Markdown document includes data preprocessing, statistical modeling, visualization, and epidemiological interpretation.

## Finished Reports
The repository includes two completed reports summarizing the analysis, results, and conclusions:

PDF Report: A full report containing analysis details, figures, and interpretations.
HTML Report: An interactive HTML version of the report for dynamic viewing.

## Running the Analysis
To reproduce the results, follow these steps to set up your local environment.

### Prerequisites
1. Install Git on your system.
2. Install R and RStudio.

### Installation
**1. Clone the repository:**

First, clone the repository to your local machine: 
1. Open Git Bash or the terminal application
2. Navigate to where you want to clone the repo to _(e.g. Documents)_:
```bash
cd ~/Documents
```
3. Clone the repo to the specificed location:
```bash
git clone [https://github.com/Biology3579/SARSCoV2Assignment.git](https://github.com/Biology3579/SARSCoV2Assignment.git) 
```
This will create a local copy of the repository on your machine.

4. After cloning, navigate to the project directory:
```bash
cd SARSCoV2Assignment
```

**2. Open PenguinAnalysis.Rmd**

Find the `SARSCoV2Assignment.Rmd` file in the main repo.
This is the main file for the analysis.

**3. Restore the R environment**

To ensure the correct R packages and dependencies are installed, run the following command in R:

```r
renv::restore()
```
This will:

- Install the exact package versions specified in the renv.lock file.
- Set up the project environment to match the original computational setup.
*Note: this requires yoi to have renv previously installed.*
If you don't have renv already installed, install it by running the following command in your R console:
```r
install.packages("renv")
```

## Packages
_Packages Used:_
The following R packages are required for the analysis:

*Data Handling & Processing:*
here – Manages project directory paths.
readr – Reads and writes CSV files.
janitor – Cleans and standardizes data.
tidyverse – A collection of data science packages for efficient manipulation and visualization.
dplyr – Data manipulation (filtering, grouping, summarizing).
lubridate – Handles date and time data.
Statistical Modeling & Epidemiology
nls.multstart – Fits logistic growth models for variant analysis.
nlstools – Provides tools for analyzing non-linear regression models.
EpiEstim – Estimates the time-varying reproduction number ($R_t$).
incidence – Analyzes incidence curves in epidemiological studies.

*Visualization*
ggplot2 – Generates all plots and visualizations.
scales – Adjusts axis scales and formats labels.
patchwork – Arranges multiple plots together.
viridis – Provides colorblind-friendly color palettes.
ggpubr – Adds statistical annotations to plots.
Report Formatting
knitr – Converts R Markdown into reports.
tinytex – Generates PDFs from R Markdown.

*Reproducibility*
renv – Manages the R environment, ensuring reproducibility of package versions.

These packages are managed through `renv`, and the necessary versions are specified in the `renv.lock` file.

##Project Structure
```bash
SARSCoV2Assignment/
│-- data/                     # Contains raw and processed datasets
│-- functions/                 # Custom R functions for cleaning, modeling, and plotting
│-- reports/                   # Contains the final HTML and PDF reports
│-- SARSCoV2Analysis.Rmd        # Main analysis script
│-- renv.lock                   # Package dependencies for reproducibility
│-- README.md                   # Project documentation (this file)
```

