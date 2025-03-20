## Script name: plotting.R
##
## Purpose of script: 
##    A collection of functions for plotting and saving figures for COVID-19 
##    variant trends, decluttering the main analysis file.
##
## Author: [Your Name]
##
## Date Created: 2025-03-18 
##
## ---------------------------

# Required libraries
library(ggplot2)
library(dplyr)
library(svglite)
library(patchwork)
library(tidyr)  # Needed for unnest()

# Custom colour palettes for coloured plots ----
## These palettes are adjusted from the Viridis palette from viridisLite, which are designed 
## to be colour-blind friendly and are perceptually uniform in both colour and black and white 
## I have adjusted slighted for better contrast 

# Colour palette for major SARS-CoV-2 variants (High Contrast)
variant_colours <- c(
  "B.1.1.7"   = "#440154",  # Deepest Purple (Darkest)
  "B.1.617.2" = "#433E85",  # Dark Indigo-Blue
  "BA.1"      = "#31688E",  # Blue-Teal
  "BA.2"      = "#26828E",  # Teal-Green
  "BA.2.75"   = "#1F9E89",  # Emerald Green
  "BA.4"      = "#35B779",  # Bright Green
  "BA.5"      = "#6DCD59",  # Lime-Green
  "BA.5.3"    = "#B8DE29",  # Yellow-Green
  "Other"     = "#E7E419",  # Golden Yellow
  "XBB"       = "#FDE725"   # Brightest Yellow (Lightest)
)


# Colour Palette for England Regions (High Contrast)
region_colours <- c(
  "East Midlands"           = "#440154",  # Deep Purple
  "East of England"         = "#404788",  # Dark Blue
  "London"                 = "#287D8E",  # Cyan-Teal
  "North East"             = "#1F9E89",  # Emerald Green
  "North West"             = "#35B779",  # Bright Green
  "South East"             = "#6DCD59",  # Lime Green
  "South West"             = "#B5DE2B",  # Yellow-Green
  "West Midlands"          = "#E7E419",  # Bright Yellow
  "Yorkshire and Humber"   = "#FDE725"   # Lightest Yellow
)



# Plotting functions for stacked area plots (Q1.3) ----
## These functions...

# Function to plot stacked area plot of total variant counts
plot_total_variant_counts <- function(data) {
  ggplot(data, aes(x = collection_date, y = total_count, fill = variant)) +
    geom_area(position = "stack") +  # Ensure proper stacking of counts
    scale_fill_manual(values = variant_colours) +  # Apply colourblind-friendly palette
    labs(
      title = "Total Counts of Major variants Over Time",
      x = "Collection Date",
      y = "Total Count",
      fill = "variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
}


# Function to plot stacked area plot of variant frequencies
plot_variant_frequencies <- function(data) {
  ggplot(data, aes(x = collection_date, y = variant_frequency, fill = variant)) +
    geom_area(position = "fill") +  # Stacked area plot (normalized)
    scale_fill_manual(values = variant_colours) +  # Apply colourblind-friendly palette
    labs(
      title = "Frequency of Major Variants Over Time",
      x = "Collection Date",
      y = "Proportion",
      fill = "variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
}

# Plotting functions for (Q2.1) ----

plot_ba2_frequency_comparison <- function(ba2_sanger, ba2_onscis) {
  ggplot() +
    # Sanger Data (Weekly)
    geom_line(data = ba2_sanger, 
              aes(x = collection_date, y = variant_frequency, color = source, linetype = source), 
              linewidth = 1.2) +
    
    # ONS-CIS Data (10-day Binned)
    geom_line(data = ba2_onscis, 
              aes(x = collection_date, y = variant_frequency, color = source, linetype = source), 
              linewidth = 1.2) + 
    
    scale_color_manual(values = c("Sanger (Weekly)" = "#31688E",  # Dark Purple
                                  "ONS-CIS (10-day Binned)" = "#1F9E89")) +  # Teal Green
    scale_linetype_manual(values = c("Sanger (Weekly)" = "solid", 
                                     "ONS-CIS (10-day Binned)" = "solid")) +
    
    # Labels & Theme
    labs(
      title = "BA.2 Frequency Trajectory (sanger vs ONS-CIS)",
      x = "Collection Date",
      y = "Frequency",
      color = "Dataset Source",
      linetype = "Dataset Source"
    ) +
    
    # Adjust x-axis to focus on the main timeframe of interest
    scale_x_date(limits = c(as.Date("2021-10-30"), as.Date("2022-10-01")), date_labels = "%Y-%m") +
    
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Plotting functions for ... (Q3) ----

# 
plot_variant_frequency_trajectories <- function(variants) {
  # Filter data for selected variants
  filtered_data <- sanger_analysis_data %>%
    filter(variant %in% variants)
  
  # Generate plot
  ggplot(filtered_data, aes(x = collection_date, y = variant_frequency, color = variant)) +
    geom_line(linewidth = 1.2) +  # Line for trends
    geom_point(size = 2, alpha = 0.7) +  # Points for visibility
    scale_color_manual(values = variant_colours) +
    labs(
      title = paste("Frequency Trajectories of", paste(variants, collapse = ", ")),
      x = "Collection Date",
      y = "Frequency",
      color = "Variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal x-axis labels
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}




# ..
plot_logistic_growth_facet <- function(data, growth_phases, variants) {
  # Create lists to store observed and predicted data
  observed_list <- list()
  predicted_list <- list()
  growth_rates <- list()  # Store growth rates for annotation
  
  for (variant in variants) {
    # Extract growth phase dates for the given variant
    phase_dates <- growth_phases %>%
      filter(variant == !!variant)  # ✅ Correct filtering
    
    # Subset the data to the variant's growth phase
    variant_growth_phase <- data %>%
      filter(variant == !!variant, 
             collection_date >= phase_dates$start_date, 
             collection_date <= phase_dates$end_date) %>%
      mutate(variant = variant)  # Add variant column for faceting
    
    # Apply the generalized logistic growth function
    logistic_results <- fit_logistic_growth_general(
      variant_growth_phase, 
      time_col = "collection_date", 
      frequency_col = "variant_frequency", 
      group_col = "variant"
    )
    
    # Skip variant if model failed
    if (is.null(logistic_results)) next
    
    # Extract growth rate (`s`) and initial frequency (`f0`)
    growth_rate <- unique(logistic_results$s)
    
    # Store growth rate for annotation
    growth_rates[[variant]] <- data.frame(
      variant = variant,
      s = round(growth_rate, 4),  # Round for readability
      x_pos = max(variant_growth_phase$collection_date),  # Position at rightmost date
      y_pos = min(variant_growth_phase$variant_frequency, na.rm = TRUE) * 0.5  # Position near bottom
    )
    
    # Store observed and predicted data
    observed_list[[variant]] <- variant_growth_phase
    predicted_list[[variant]] <- logistic_results
  }
  
  # Combine all observed and predicted data
  observed_data <- bind_rows(observed_list)
  predicted_data <- bind_rows(predicted_list)
  growth_rate_data <- bind_rows(growth_rates)  # Convert growth rate list to dataframe
  
  # Create the faceted plot
  ggplot(observed_data, aes(x = collection_date, y = variant_frequency, color = variant)) +
    geom_point(size = 2, alpha = 0.7) +  # Observed data
    geom_line(data = predicted_data, aes(x = collection_date, y = predicted_frequency, color = variant), size = 1) +  
    geom_text(data = growth_rate_data, aes(x = x_pos, y = y_pos, label = paste0("s = ", s)), 
              color = "black", size = 3.5, fontface = "bold", hjust = 1, vjust = 0) +  # ✅ Small text, bottom-right
    scale_color_manual(values = variant_colours) +  # Use predefined colour scheme
    scale_x_date(date_labels = "%Y-%m") +  # ✅ Keep dates in YYYY-MM format
    labs(
      title = "Logistic Growth Fit for SARS-CoV-2 Variants",
      x = "Collection Date",
      y = "Frequency",
      color = "Variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # ✅ Keep readable x-axis labels
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    facet_grid(. ~ variant, scales = "free_x")  # ✅ Side-by-side facets with independent x-axis ranges
}


# Plotting functions for ... (Q4) ----

#Function to ...
plot_delta_frequencies <- function(data) {
  ggplot(data, aes(x = week, y = delta_frequency, color = phecname)) +  # Change 'date' -> 'week'
    geom_line(size = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = region_colours) +
    labs(
      title = "Delta Variant Frequency Over Time (Weekly Smoothed)",
      x = "Collection Date",
      y = "Frequency of Delta",
      color = "Region"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_logistic_growth_region <- function(observed_data, predicted_data) {
  
  # Ensure required columns exist
  if (!all(c("week", "delta_frequency", "phecname") %in% colnames(observed_data))) {
    stop("Error: observed_data must contain 'week', 'delta_frequency', and 'phecname' columns.")
  }
  
  if (!all(c("week", "predicted_frequency", "phecname", "s", "f0") %in% colnames(predicted_data))) {
    stop("Error: predicted_data must contain 'week', 'predicted_frequency', 'phecname', 's', and 'f0' columns.")
  }
  
  # Extract growth rate (`s`) and initial frequency (`f0`) for annotation
  growth_info <- predicted_data %>%
    group_by(phecname) %>%
    slice_max(week, n = 1) %>%  # Get the latest date per region
    mutate(
      label = paste0("s = ", round(s, 3), "\nf0 = ", round(f0, 3)),  # Format text for both `s` and `f0`
      x_pos = week,
      y_pos = 0.05  # Fixed position near the bottom
    )
  
  # Generate plot
  ggplot(observed_data, aes(x = week, y = delta_frequency, color = phecname)) +
    geom_point(size = 1.8, alpha = 0.7) +  # Observed frequencies
    geom_line(data = predicted_data, aes(x = week, y = predicted_frequency, group = phecname, color = phecname),
              size = 1) +  # Logistic growth curves
    geom_text(data = growth_info,  
              aes(x = x_pos, y = y_pos, label = label),  
              color = "black", size = 3, fontface = "bold", hjust = 1, vjust = 0) +  # ✅ Small text, stacked s & f0
    scale_color_viridis_d() +  # Colorblind-friendly palette
    labs(
      title = "Logistic Growth Fit for Delta Variant Frequency by Region",
      x = "Collection Date (Week)",
      y = "Frequency of Delta",
      color = "Region"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~phecname, scales = "free_x")  # ✅ Facet by region, separate x-axis ranges
}


# Function to save plots as SVG ----
save_plot_svg <- function(data, filename, size, scaling, plot_function, ...) {
  size_inches <- size / 2.54  # Convert size from cm to inches
  svglite(filename, width = size_inches, height = size_inches, scaling = scaling)
  
  plot <- plot_function(data, ...)  # Call the respective plot function
  print(plot)  # Print the plot to save it
  
  dev.off()  # Close the device
}
