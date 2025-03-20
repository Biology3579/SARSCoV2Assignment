## Script name: cleaning_and_processing_data.R
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



# Functions to process the sanger raw data, making it ... ----
#     - First I ...
#     - Then ...
# First, to simplify the long list of variant names assigned by the Pango nomenclature, we group them into broader ‘major variants.’ For this practical, we focus on variants that caused significant waves in the UK since late 2020. 
# These include: Alpha (B.1.1.7), Delta (B.1.617.2), and various Omicron subvariants, including BA.1, BA.2, BA.4, BA.5, and XBB. 
# We put variants from other variant into the 'Other' category.

# Function to clean Sanger data
clean_sanger_data <- function(raw_data, variants) {
  raw_data %>%
    rename(collection_date = date, variant = lineage) %>%  # Rename lineage to variant
    mutate(
      collection_date = as.Date(collection_date),  # Convert date column to Date format
      variant = as.character(variant),  # Ensure variant is a character variable
      count = as.numeric(count),  # Ensure count is numeric
      variant = ifelse(variant %in% variants, variant, "Other")  # Classify major variants
    )
}

# Function to calculate daily counts and frequencies of major variants over time
counts_and_frequencies <- function(cleaned_data) {
  cleaned_data %>%
    group_by(collection_date, variant) %>%
    summarise(variant_count = sum(count, na.rm = TRUE), .groups = "drop") %>%  # Aggregate variant counts
    group_by(collection_date) %>%
    mutate(
      total_count = sum(variant_count),  # Compute total samples per date
      variant_frequency = variant_count / total_count  # Compute frequency per variant
    ) %>%
    ungroup()
}


#Functions to ... ----

# Function to clean ONSCIS data
clean_onscis_data <- function(raw_data) {
  raw_data %>%
    mutate(
      collection_date = as.Date(collection_date),  # Convert to Date format
      variant = as.character(major_lineage)  # Rename major_lineage to variant
    ) %>%
    group_by(collection_date, variant) %>%
    summarise(count = n(), .groups = "drop") %>%  # Aggregate counts
    ungroup()
}

# Function to bin data and calculate frequencies
bin_and_calculate_frequencies <- function(cleaned_data, bin_size = 10) {
  cleaned_data %>%
    mutate(
      collection_date = as.Date(
        floor(as.numeric(collection_date) / bin_size) * bin_size,
        origin = "1970-01-01"  # Aligns binning with the epoch
      )
    ) %>%
    
    # Aggregate after binning
    group_by(collection_date, variant) %>%
    summarise(variant_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    
    # Compute total count per bin for frequency calculation
    group_by(collection_date) %>%
    mutate(
      total_count = sum(variant_count, na.rm = TRUE),
      variant_frequency = variant_count / total_count
    ) %>%
    ungroup()
}


#Functions to ... ----


# Function to clean Delta data
clean_delta_data <- function(raw_data) {
  raw_data %>%
    filter(!is.na(phecname) & phecname != "") %>%  # Remove missing region names
    mutate(
      date = as.Date(date),  # Convert to Date format
      Delta = as.numeric(Delta)  # Convert logical to numeric (1 = TRUE, 0 = FALSE)
    )
}

# Function to calculate counts and frequencies of Delta variant
counts_and_frequencies_delta <- function(cleaned_data) {
  cleaned_data %>%
    group_by(phecname, date) %>%  # Group by region and date
    summarise(
      total_samples = n(),  # Compute total samples per region per date
      delta_count = sum(Delta, na.rm = TRUE),  # Count Delta variant occurrences (handle NAs)
      delta_frequency = ifelse(total_samples > 0, delta_count / total_samples, NA),  # Compute frequency safely
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
