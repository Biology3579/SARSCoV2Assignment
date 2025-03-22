## Script name: cleaning_and_processing_data.R
##
## Purpose of script: 
##    A collection of functions for cleaning, curating, and processing 
##    SARS-CoV-2 sequencing and epidemiological data for downstream analysis.
##
## Author: Biology3579
##
## Date Created: 2025-03-18 
##
## ---------------------------

# Load required libraries ----
library(dplyr)       # Data manipulation (mutate, filter, group_by, summarise, rename)
library(tidyr)       # Data tidying (handling missing values, reshaping data)
library(lubridate)   # Date handling (as.Date, floor_date)

# ---------------------------
# Cleaning Sanger Data ----
# The Sanger dataset contains sequencing data for SARS-CoV-2 variants.
# This function standardizes and filters the raw Sanger data for downstream analysis.
#
# Steps:
#   1. Rename columns for clarity.
#   2. Convert relevant columns to appropriate data types (date, character, numeric).
#   3. Classify lineages into broader variant categories (e.g., Alpha, Delta, Omicron subvariants).
#   4. Ensure consistency in variant naming and handle unknown variants by grouping them as "Other".

clean_sanger_data <- function(raw_data, variants) {
  raw_data %>%
    rename(collection_date = date, variant = lineage) %>%  # Rename for clarity
    mutate(
      collection_date = as.Date(collection_date),  # Convert date column to Date format
      variant = as.character(variant),  # Ensure variant is a character variable
      count = as.numeric(count),  # Ensure count is numeric
      variant = ifelse(variant %in% variants, variant, "Other")  # Categorize variants
    )
}

# ---------------------------
# Calculating Variant Frequencies ----
# This function calculates the daily frequency of each variant over time.
# 
# Steps:
#   1. Group data by date and variant.
#   2. Compute the total count of each variant per day.
#   3. Normalize by the total sequencing count to obtain variant frequencies.

counts_and_frequencies <- function(cleaned_data) {
  cleaned_data %>%
    group_by(collection_date, variant) %>%
    summarise(variant_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    
    # Calculate total number of sequences per date 
    group_by(collection_date) %>%
    mutate(
      total_count = sum(variant_count, na.rm = TRUE),  # total for that date
      variant_frequency = variant_count / total_count   # relative frequency
    ) %>%
    ungroup()
}


# ---------------------------
# Cleaning ONS-CIS Data ----
# The ONS-CIS dataset contains epidemiological data on SARS-CoV-2 variant prevalence.
# This function prepares the dataset for frequency analysis.
#
# Steps:
#   1. Convert `collection_date` to Date format.
#   2. Ensure `variant` is stored as a character variable.
#   3. Aggregate count data by date and variant.

clean_onscis_data <- function(raw_data) {
  raw_data %>%
    mutate(
      collection_date = as.Date(collection_date),  # Convert to Date format
      variant = as.character(major_lineage)  # Rename major_lineage to variant
    ) %>%
    group_by(collection_date, variant) %>%
    summarise(count = n(), .groups = "drop") %>%  # Aggregate counts per variant
    ungroup()
}

# ---------------------------
# Binning Data and Calculating Frequencies ----
# Some analyses require binning data into time intervals (e.g., 10-day bins).
# This function groups data into time bins and recalculates variant frequencies.
#
# Steps:
#   1. Convert dates into time bins based on `bin_size` (default: 10 days).
#   2. Sum counts per variant in each bin.
#   3. Compute relative frequencies within each bin.

bin_and_calculate_frequencies <- function(cleaned_data, bin_size = 10) {
  cleaned_data %>%
    mutate(
      collection_date = as.Date(
        floor(as.numeric(collection_date) / bin_size) * bin_size,
        origin = "1970-01-01"  # Align binning with time origin
      )
    ) %>%
    group_by(collection_date, variant) %>%
    summarise(variant_count = sum(count, na.rm = TRUE), .groups = "drop") %>%  # Aggregate counts per bin
    group_by(collection_date) %>%
    mutate(
      total_count = sum(variant_count, na.rm = TRUE),
      variant_frequency = variant_count / total_count  # Compute frequencies
    ) %>%
    ungroup()
}

# ---------------------------
# Cleaning and Processing Delta Data ----
# This function prepares Delta variant data from raw input files.
#
# Steps:
#   1. Remove rows with missing regional information (`phecname`).
#   2. Convert `date` to Date format.
#   3. Convert logical Delta presence into numeric (1 = Delta detected, 0 = not detected).

clean_delta_data <- function(raw_data) {
  raw_data %>%
    filter(!is.na(phecname) & phecname != "") %>%  # Remove missing region names
    mutate(
      date = as.Date(date),  # Convert to Date format
      Delta = as.numeric(Delta)  # Convert logical to numeric
    )
}

# ---------------------------
# Calculating Delta Variant Frequencies ----
# This function computes the proportion of Delta cases within each region over time.
#
# Steps:
#   1. Aggregate total cases and Delta cases per region and date.
#   2. Compute the proportion of Delta cases per region.

counts_and_frequencies_delta <- function(cleaned_data) {
  cleaned_data %>%
    group_by(phecname, date) %>%  # Group by region and date
    summarise(
      total_samples = n(),  # Total samples per region per date
      delta_count = sum(Delta, na.rm = TRUE),  # Count Delta occurrences
      delta_frequency = ifelse(total_samples > 0, delta_count / total_samples, NA),  # Compute frequency
      .groups = "drop"
    )
}

