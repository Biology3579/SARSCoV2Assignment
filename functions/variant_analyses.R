## Script name: variant_analyses.R
##
## Purpose of script: 
##    A collection of functions for analyzing the spread and growth dynamics 
##    of SARS-CoV-2 variants, including:
##    - Extracting growth phase periods
##    - Fitting logistic growth models
##    - Estimating daily Delta cases
##    - Calculating the reproduction number (Rt)
##
## Author: Biology3579
##
## Date Created: 2025-03-18
##
## ---------------------------


# Load required libraries ----
library(dplyr)        # Data wrangling and manipulation
library(tidyr)        # Data tidying functions
library(ggplot2)      # Data visualization (if needed for plotting)
library(nlstools)     # Non-linear regression tools (for logistic growth model)
library(lubridate)    # Date handling and manipulation
library(purrr)        # Functional programming (for group operations)
library(EpiEstim)     # Estimating R_t (Reproduction number)

# ---------------------------
# Extracting Growth Phase Dates ----
# This function determines the start and end of a variant's growth phase.
#
# Steps:
#   1. Identify the first date where the variant is consistently detected.
#   2. Find the peak frequency date, marking the end of the growth phase.

extract_growth_phase_dates <- function(data, variant) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(variant == !!variant) %>%
    arrange(collection_date) %>%
    mutate(
      is_nonzero = variant_frequency > 0,  # Logical column for nonzero frequencies
      consecutive_nonzero = is_nonzero & lag(is_nonzero, default = FALSE) & lead(is_nonzero, default = FALSE) # Check for consecutive nonzero
    )
  
  # Identify the first date where the variant is consistently detected
  start_date <- variant_data %>%
    filter(consecutive_nonzero) %>%
    slice_head(n = 1) %>%
    pull(collection_date)
  
  # Find the maximum observed frequency for this variant
  max_frequency <- max(variant_data$variant_frequency, na.rm = TRUE)
  
  # Identify the first date when the frequency reaches this maximum value
  end_date <- variant_data %>%
    filter(variant_frequency == max_frequency) %>%
    arrange(collection_date) %>%
    slice_head(n = 1) %>%
    pull(collection_date)
  
  # Return a tibble with extracted dates
  return(tibble(variant = variant, start_date = start_date, end_date = end_date))
}

# ---------------------------
# Fitting Logistic Growth Models ----
# This function fits a logistic growth model to variant frequency data.
#
# Steps:
#   1. Sort the dataset by the time variable.
#   2. Define the time variable relative to the first observation.
#   3. Fit a logistic growth model using non-linear least squares (nls).
#   4. Extract growth rate (`s`) and initial frequency (`f0`).
#   5. Generate smooth logistic predictions over the given timeframe.

fit_logistic_growth_general <- function(data, time_col, frequency_col, group_col) {
  
  # Sort data by the specified time column
  data <- data %>% arrange(.data[[time_col]])
  
  # Define time since first observation
  data <- data %>%
    mutate(time = as.numeric(.data[[time_col]] - min(.data[[time_col]], na.rm = TRUE)))
  
  # Fit logistic model using nls()
  nls_fit <- tryCatch(
    nls(
      formula = as.formula(paste0(frequency_col, " ~ logistic_growth(time, s, f0)")),
      data = data,
      start = list(s = 0.1, f0 = min(data[[frequency_col]][data[[frequency_col]] > 0], na.rm = TRUE)),
      control = nls.control(maxiter = 100)
    ),
    error = function(e) NULL # Handle errors gracefully
  )
  
  # If the model fails, return NULL
  if (is.null(nls_fit)) return(NULL)
  
  # Extract growth rate (`s`) and initial frequency (`f0`)
  coef_values <- coef(nls_fit)
  growth_rate <- coef_values["s"]
  f0 <- coef_values["f0"]
  
  # Generate smooth logistic predictions
  smooth_dates <- seq(min(data[[time_col]], na.rm = TRUE), max(data[[time_col]], na.rm = TRUE), by = "1 day")
  smooth_predictions <- tibble(
    !!time_col := as.Date(smooth_dates),
    predicted_frequency = logistic_growth(as.numeric(smooth_dates - min(data[[time_col]], na.rm = TRUE)), growth_rate, f0),
    !!group_col := unique(data[[group_col]]),
    s = growth_rate,  # Store growth rate
    f0 = f0           # Store initial frequency
  )
  
  return(smooth_predictions)
}

# ---------------------------
# Estimating Daily Delta Cases ----
# This function estimates the number of Delta variant cases in England.
#
# Steps:
#   1. Align daily case data with weekly Sanger sequencing data.
#   2. Assign each daily case count to the closest past weekly Delta proportion.
#   3. Multiply the 7-day averaged case counts by the Delta frequency to estimate cases.

estimate_delta_cases <- function(daily_cases, sanger_data, start_date = "2020-09-05") {
  
  # Filter daily data to start from first available week in Sanger dataset
  daily_cases_filtered <- daily_cases %>%
    filter(date >= start_date) %>%
    mutate(date = as.Date(date))  # Ensure correct date format
  
  # Extract only Delta variant proportions from Sanger dataset
  sanger_delta <- sanger_data %>%
    filter(variant == "B.1.617.2") %>%  # Keep only Delta variant
    select(collection_date, variant_frequency) %>%
    rename(week = collection_date)  # Rename for clarity
  
  # Assign each daily case count to the closest past weekly proportion
  daily_cases_with_delta <- daily_cases_filtered %>%
    mutate(
      week = sanger_delta$week[findInterval(date, sanger_delta$week)]  # Match each date to the closest previous week
    ) %>%
    left_join(sanger_delta, by = "week") %>%  # Merge weekly Delta proportions with daily case counts
    mutate(
      estimated_delta_cases = round(cases_sevendayaveraged * variant_frequency, 0)  # Estimate Delta cases (rounded)
    )
  
  # Convert `week` back to Date format for consistency
  daily_cases_with_delta$week <- as.Date(daily_cases_with_delta$week)
  
  return(daily_cases_with_delta)
}

# ---------------------------
# Estimating the Time-Varying Reproduction Number (Rt) ----
# This function estimates Rt for the Delta variant using EpiEstim.
#
# Steps:
#   1. Filter Delta case estimates for the fixed Rt estimation period.
#   2. Define the serial interval parameters for Delta.
#   3. Use EpiEstim to compute Rt over the selected time window.
#   4. Return the estimated reproduction number.

generate_rt_estimates <- function(delta_estimates) {
  
  # Use the specified start and end dates
  start_date <- as.Date("2021-04-23")
  end_date <- as.Date("2021-11-01")
  
  # Set standard serial interval parameters for Delta variant
  mean_si <- 4.1  # Mean serial interval (days)
  std_si <- 2.8   # Standard deviation of the serial interval
  
  # Step 1: Filter estimated Delta cases within the specified date range
  delta_incidence_data <- delta_estimates %>%
    filter(date >= start_date & date <= end_date) %>%
    select(date, estimated_delta_cases) %>%
    rename(dates = date, I = estimated_delta_cases)  # Rename for compatibility with EpiEstim
  
  # Step 2: Define the serial interval parameters
  serial_interval <- list(mean_si = mean_si, std_si = std_si)
  
  # Step 3: Run EpiEstim to estimate Rt
  rt_results <- estimate_R(
    incid = delta_incidence_data,
    method = "parametric_si",
    config = make_config(serial_interval)
  )
  
  # Return only the R_t estimates
  return(rt_results$R)
}
