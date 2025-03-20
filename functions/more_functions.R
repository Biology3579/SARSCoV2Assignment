

# Fucntions to ...

extract_growth_phase_dates <- function(data, variant) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(variant == !!variant) %>%
    arrange(collection_date) %>%
    mutate(
      is_nonzero = variant_frequency > 0,  # Logical column for nonzero frequencies
      consecutive_nonzero = is_nonzero & lag(is_nonzero, default = FALSE) & lead(is_nonzero, default = FALSE) # Check for consecutive nonzero
    )
  
  # Identify the first date where `min_growth_points` consecutive nonzero values exist
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
    ))
  
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





fit_variant_growth_1 <- function(data, variant, start_date, end_date) {
  
  # Filter data for the specific variant
  variant_data <- data %>%
    filter(variant == variant)
  
  # Subset data to only include the user-defined growth phase
  growth_phase_data <- variant_data %>%
    filter(collection_date >= start_date & collection_date <= end_date) %>%
    mutate(days_since_start = as.numeric(collection_date - start_date))
  
  # Fit the logistic growth model using nls()
  nls_fit <- nls(
    variant_frequency ~ logistic_growth(days_since_start, s, f0),
    data = growth_phase_data,
    start = list(s = 0.1, f0 = min(growth_phase_data$variant_frequency, na.rm = TRUE)) # Initial guesses
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
    filter(variant == variant)
  
  # Subset data to only include the user-defined growth phase
  growth_phase_data <- variant_data %>%
    filter(collection_date >= start_date & collection_date <= end_date) %>%
    mutate(days_since_start = as.numeric(collection_date - start_date))
  
  # Fit the logistic growth model using nls()
  nls_fit <- nls(
    variant_frequency ~ logistic_growth(days_since_start, s, f0),
    data = growth_phase_data,
    start = list(s = 0.1, f0 = min(growth_phase_data$variant_frequency, na.rm = TRUE)) # Initial guesses
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
