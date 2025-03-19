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


#Function to ...
process_sanger_data <- function(raw_data, major_lineages) {
  raw_data %>%
    rename(collection_date = date, major_lineage = lineage) %>%  # Standardize column names
    mutate(
      collection_date = as.Date(collection_date),  # Ensure correct date format
      major_lineage = as.character(major_lineage),  # Ensure lineage is a character variable
      count = as.numeric(count)  # Ensure count is numeric
    ) %>%
    classify_major_lineages(major_lineages) %>%  # Identify and classify major lineages
    group_by(collection_date, major_lineage) %>%
    summarise(lineage_count = sum(count, na.rm = TRUE), .groups = "drop") %>%  # Aggregate lineage counts
    group_by(collection_date) %>%
    mutate(
      total_count = sum(lineage_count),  # Compute total samples per date
      lineage_frequency = lineage_count / total_count  # Compute frequency per lineage
    ) %>%
    ungroup()
}


#Function to ...
process_onscis_data <- function(raw_data, bin_size = 10) {
  raw_data %>%
    mutate(
      collection_date = as.Date(collection_date),  
      major_lineage = as.character(major_lineage)  
    ) %>%
    group_by(collection_date, major_lineage) %>%
    summarise(count = n(), .groups = "drop") %>%
    ungroup() %>%
    
    #Apply consistent binning aligned to 1970-01-01 (to match `aggregate()` method)
    mutate(
      collection_date_bin = as.Date(
        floor(as.numeric(collection_date) / bin_size) * bin_size,
        origin = "1970-01-01" )) %>%
    
    #Aggregate after binning
    group_by(collection_date_bin, major_lineage) %>%
    summarise(
      lineage_count = sum(count, na.rm = TRUE),
      .groups = "drop") %>%
    
    #Compute total count per bin for frequency calculation
    group_by(collection_date_bin) %>%
    mutate(
      total_count = sum(lineage_count, na.rm = TRUE),
      lineage_frequency = lineage_count / total_count) %>%
    ungroup()
}

#Fucntion
process_delta_data <- function(raw_data) {
  raw_data %>%
    remove_empty_entries() %>%  # Remove empty rows/columns
    filter(phecname != "") %>%  # Remove rows with missing region names
    mutate(
      date = as.Date(date),  # Convert to Date format
      Delta = as.numeric(Delta)  # Convert logical to numeric (1 = TRUE, 0 = FALSE)
    ) %>%
    group_by(phecname, date) %>%  # Group by region and date
    summarise(
      total_samples = n(),  # Compute total samples per region per date
      delta_count = sum(Delta),  # Count Delta variant occurrences
      delta_frequency = delta_count / total_samples,  # Compute Delta frequency
      .groups = "drop"
    )
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


# Fucntions to ...

extract_growth_phase_dates <- function(data, variant, min_growth_points = 4, fixation_threshold = 0.98) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(major_lineage == variant) %>%
    arrange(collection_date) %>%
    mutate(
      is_nonzero = lineage_frequency > 0,  # Create logical column for nonzero frequencies
      consecutive_nonzero = is_nonzero & lag(is_nonzero, default = FALSE) & lead(is_nonzero, default = FALSE) # Check for consecutive nonzero
    )
  
  # Identify the first date where `min_growth_points` consecutive nonzero values exist
  start_date <- variant_data %>%
    filter(consecutive_nonzero) %>%
    slice_head(n = 1) %>%
    pull(collection_date)
  
  # Identify the first date when frequency reaches the fixation threshold
  end_date <- variant_data %>%
    filter(lineage_frequency >= fixation_threshold) %>%
    arrange(collection_date) %>%
    slice_head(n = 1) %>%
    pull(collection_date)
  
  # If no fixation found, use the last available date
  if (length(end_date) == 0) {
    end_date <- max(variant_data$collection_date, na.rm = TRUE)
  }
  
  # Return a tibble with extracted dates
  return(tibble(variant = variant, start_date = start_date, end_date = end_date))
}


fit_variant_growth_1 <- function(data, variant, start_date, end_date) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(major_lineage == variant)
  
  # Subset data to only include the user-defined growth phase
  growth_phase_data <- variant_data %>%
    filter(collection_date >= start_date & collection_date <= end_date) %>%
    mutate(days_since_start = as.numeric(collection_date - start_date))
  
  # Fit the logistic growth model using nls()
  nls_fit <- nls(
    lineage_frequency ~ logistic_growth(days_since_start, s, f0),
    data = growth_phase_data,
    start = list(s = 0.1, f0 = min(growth_phase_data$lineage_frequency, na.rm = TRUE)) # Initial guesses
  )
  
  # Extract the fitted growth rate (s)
  growth_rate <- coef(nls_fit)["s"]
  
  # Return results as a tibble
  return(list(
    fitted_model = nls_fit, 
    growth_data = growth_phase_data, 
    growth_rate = growth_rate
  ))
}

#
fit_variant_growth <- function(data, variant, start_date, end_date) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(major_lineage == variant)
  
  # Subset data to only include the user-defined growth phase
  growth_phase_data <- variant_data %>%
    filter(collection_date >= start_date & collection_date <= end_date) %>%
    mutate(days_since_start = as.numeric(collection_date - start_date))
  
  # Fit the logistic growth model using nls()
  nls_fit <- nls(
    lineage_frequency ~ logistic_growth(days_since_start, s, f0),
    data = growth_phase_data,
    start = list(s = 0.1, f0 = min(growth_phase_data$lineage_frequency, na.rm = TRUE)) # Initial guesses
  )
  
  # Extract the fitted growth rate (s)
  growth_rate <- coef(nls_fit)["s"]
  
  # Return results as a tibble
  return(list(
    fitted_model = nls_fit, 
    growth_data = growth_phase_data, 
    growth_rate = growth_rate
  ))
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
