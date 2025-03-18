## Script name: cleaning_and_curating.R
##
## Purpose of script: 
##    A collection of functions for cleaning and curating the ...
## Author: Biology3579
##
## Date Created: 2025-03-18 
##
## ---------------------------
#
# Load necessary libraries
library(dplyr)
library(janitor)


# Cleaning ----
# This set of functions helps to ...

# Function to standardize column names and rename key columns
clean_col_names <- function(raw_data) {
  raw_data %>%
    clean_names() %>%  # Converts column names to snake_case
    rename(
      collection_date = date) # Rename date for clarity
}

# A function to remove any empty columns or rows
remove_empty_entries <- function(raw_data) {
  raw_data %>%
    remove_empty(c("rows", "cols"))
}

# A function to ensure all the variables are in the correct format 
clean_data_format <- function(raw_data) {
  raw_data %>%
    mutate(
      collection_date = as.Date(collection_date),
      lineage = as.character(lineage),
      count = as.numeric(count)
    )
}

# Function to identify and rename the variants of focus as the major lineages
classify_major_lineages <- function(raw_data, major_lineages) {
  raw_data %>%
    mutate(lineage = ifelse(lineage %in% major_lineages, lineage, "Other"))
}

# Processing cleaned data ----
# Here I curate the clean penguins data set specifically for my analysis:
#   - I rename columns to make variable names more intuitive and analysis-friendly.
#   - I filter out rows with missing values (NAs)
#   - I ensure that the class type of key variables (e.g., numeric or factor) is correct 
#     to avoid errors during analysis.
#   - I remove variables that are not relevant to my specific analysis (e.g., administrative 
#     or study-specific details) to reduce noise and streamline the dataset.
#   - I create a new variable, `bill_morphology`, which represents the ratio of bill length
#     to bill depth. This new variable is a critical metric for understanding differences 
#     in bill shape across penguin species and serves as a key indicator in subsequent 
#     statistical analyses.
#This step is performed separately from the cleaning process to ensure that the 
#original cleaned dataset remains intact and is not overwritten. 
#This allows for easy reuse of the cleaned data for further analyses.


# Function to aggregate lineage counts over time with 10-day binning.
# # The function enables to sum all the different strains at a given date, which you can then track over time.
curate_lineage_trend_data <- function(data) {
  data %>%
    group_by(collection_date, lineage) %>%  # Aggregate by date and lineage
    summarise(
      lineage_count = sum(count, na.rm = TRUE),  # Total count per lineage per day/week
      .groups = "drop"
    ) %>%
    group_by(collection_date) %>%
    mutate(
      total_count = sum(lineage_count),  # Compute total samples per date
      lineage_frequency = lineage_count / total_count  # Compute frequency per lineage
    ) %>%
    ungroup()
}


# Function to drop specific columns based on a provided list of column names
drop_cols <- function(clean_data, columns_names) {
  clean_data %>%
    select(-all_of(columns_names))
}

# Function to rename the columns for ease of use
rename_columns <- function(clean_data) {
  clean_data %>%
    rename(
      bill_length_mm = culmen_length_mm,  # Rename culmen_length_mm to bill_length_mm
      bill_depth_mm = culmen_depth_mm)   # Rename culmen_depth_mm to bill_depth_mm
}

# Function to filter out rows with missing values (NAs)
remove_NA <- function(raw_data) {
  raw_data %>%
    na.omit()
}

# Function for converting variables to their appropriate types
convert_variables <- function(clean_data) {
  clean_data %>%
    mutate(
      bill_length_mm = as.numeric(bill_length_mm),
      bill_depth_mm = as.numeric(bill_depth_mm),
      body_mass_g = as.numeric(body_mass_g),
      species = as.factor(species),
      sex = as.factor(tolower(sex))) # Convert to lowercase and make it a factor
}

# Function to calculate and add the bill morphology (length-to-depth ratio) to the dataset
add_bill_morphology <- function(clean_data) {
  clean_data %>%
    mutate(bill_morphology = bill_length_mm / bill_depth_mm)  # Add the bill morphology (length-to-depth ratio)
}

# A unified curating function combining all of the curating functions for simplicity of use
curating_penguins_clean <- function(clean_data){
  clean_data %>%
    drop_cols(c(
      "study_name",
      "region",
      "stage", 
      "flipper_length_mm",
      "sample_number",
      "clutch_completion",
      "date_egg",
      "delta_15_n_o_oo",
      "delta_13_c_o_oo",
      "comments")) %>%
    rename_columns() %>%
    remove_NA() %>%
    convert_variables()  %>%  
    add_bill_morphology()  
}
