## Script name: plotting.R
##
## Purpose of script: 
##    This script contains a collection of functions for plotting and saving figures
##    related to COVID-19 variant trends. It is designed to declutter the main analysis
##    file by keeping all visualisation-related functions separate.
##
## Author: Biology3579
##
## Date Created: 2025-03-18 
##
## ---------------------------

# Load required libraries ----
library(ggplot2)      # Data visualization (ggplot-based plotting)
library(dplyr)        # Data wrangling (mutate, filter, group_by, summarise, arrange)
library(tidyr)        # Data tidying (pivoting and restructuring data)
library(lubridate)    # Date handling and manipulation
library(purrr)        # Functional programming (map functions for iteration)
library(viridis)      # Color palettes for colorblind-friendly visuals
library(scales)       # Formatting axes (date labels, continuous scale adjustments)
library(svglite)      # Exporting plots in SVG format

# Define colour palettes ----
## Custom colour palettes adapted from the Viridis palette (viridisLite)
## These palettes are colour-blind friendly and maintain perceptual uniformity
## across both coloured and black-and-white displays.

# Colour palette for major SARS-CoV-2 variants (High Contrast)
variant_colours <- c(
  "B.1.1.7"   = "#440154",  # Deep Purple (Darkest)
  "B.1.617.2" = "#433E85",  # Dark Indigo-Blue
  "BA.1"      = "#31688E",  # Blue-Teal
  "BA.2"      = "#26828E",  # Teal-Green
  "BA.2.75"   = "#1F9E89",  # Emerald Green
  "BA.4"      = "#35B779",  # Bright Green
  "BA.5"      = "#6DCD59",  # Lime-Green
  "BA.5.3"    = "#B8DE29",  # Yellow-Green
  "Other"     = "#E7E419",  # Golden Yellow
  "XBB"       = "#FDE725"   # Bright Yellow (Lightest)
)

# Colour palette for regions in England (High Contrast)
region_colours <- c(
  "East Midlands"          = "#440154",  # Deep Purple
  "East of England"        = "#404788",  # Dark Blue
  "London"                 = "#287D8E",  # Cyan-Teal
  "North East"             = "#1F9E89",  # Emerald Green
  "North West"             = "#35B779",  # Bright Green
  "South East"             = "#6DCD59",  # Lime Green
  "South West"             = "#B5DE2B",  # Yellow-Green
  "West Midlands"          = "#E7E419",  # Bright Yellow
  "Yorkshire and Humber"   = "#FDE725"   # Lightest Yellow
)

# Plotting functions ----
## Functions to generate key plots related to COVID-19 variants.

# Function to plot a stacked area plot of total weekly variant counts over time ----
# This function plots the total number of sars-cov-2 variant cases over time, 
# using a stacked area chart to show relative contributions of different variants.
plot_total_variant_counts <- function(data) {
  ggplot(data, aes(x = collection_date, y = variant_count, fill = variant)) + 
    geom_area(position = "stack") +  # Stacked area plot to show variant counts over time
    scale_fill_manual(values = variant_colours) +  # Use predefined color scheme for variants
    scale_y_continuous(labels = scales::label_number(accuracy = 1)) +  # Format y-axis to display full numbers
    scale_x_date(
      date_breaks = "2 month",  # Add tick marks every 2 months
      date_minor_breaks = "1 month",  # Add minor tick marks every month
      date_labels = "%Y-%m"  # Format date labels as YYYY-MM
    ) +
    labs(
      title = "Total Counts of Major Variants Over Time",  # Plot title
      x = "Collection Date",  # X-axis label
      y = "Total Count",  # Y-axis label
      fill = "Variant"  # Legend title
    ) +
    theme_minimal() +  # Apply a clean minimalistic theme
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # Diagonal text for better readability
      panel.grid.minor.x = element_line(color = "lightgray", size = 0.3)  # Make minor grid lines visible
    )
}


# Function to plot a stacked area plot of variant frequencies over time ----
# This function plots the proportional representation of each sars-cov-2 variant over time.
# instead of showing absolute counts, it normalizes the data so that the total proportion is always 1.
plot_variant_frequencies <- function(data) {
  ggplot(data, aes(x = collection_date, y = variant_frequency, fill = variant)) +
    geom_area(position = "fill") +  # Normalised stacked area plot to show variant proportions
    scale_fill_manual(values = variant_colours) +  # Use predefined color scheme for variants
    scale_x_date(
      date_breaks = "2 month",  # Add tick marks every 2 months
      date_minor_breaks = "1 month",  # Add minor tick marks every month
      date_labels = "%Y-%m"  # Format date labels as YYYY-MM
    ) +
    labs(
      title = "Frequency of Major Variants Over Time",  # Plot title
      x = "Collection Date",  # X-axis label
      y = "Proportion",  # Y-axis label (since proportions are used instead of counts)
      fill = "Variant"  # Legend title
    ) +
    theme_minimal() +  # Apply a clean minimalistic theme
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # Diagonal text for better readability
      panel.grid.minor.x = element_line(color = "lightgray", size = 0.3)  # Make minor grid lines visible
    )
}

# Function to compare ba.2 frequency between Sanger and ONS-CIS datasets ----
# This function creates a line plot comparing the frequency trajectories of ba.2 variant from two sources:
# - Sanger dataset (weekly sequencing data)
# - ONS-CIS dataset (10-day binned genomic survey data)
plot_ba2_frequency_comparison <- function(ba2_sanger, ba2_onscis) {
  ggplot() +
    # Plot ba.2 frequency trajectory from the Sanger dataset
    geom_line(data = ba2_sanger, aes(x = collection_date, y = variant_frequency, color = source), linewidth = 1.2) +
    # Plot ba.2 frequency trajectory from the ONS-CIS dataset
    geom_line(data = ba2_onscis, aes(x = collection_date, y = variant_frequency, color = source), linewidth = 1.2) +
    # Assign colors to differentiate the datasets
    scale_color_manual(values = c("Sanger (Weekly)" = "#31688E", "ONS-CIS (10-day Binned)" = "#1F9E89")) +
    # Set x-axis limits to focus on the period where ba.2 was dominant
    scale_x_date(
      limits = c(as.Date("2021-10-30"), as.Date("2022-10-01")),
      date_breaks = "2 month",        # Add tick marks every 2 months
      date_minor_breaks = "1 month",  # Add minor tick marks every month
      date_labels = "%Y-%m"           # Format date labels as YYYY-MM
    ) +
    labs(
      title = "BA.2 Frequency Trajectory (Sanger vs ONS-CIS)",  # plot title
      x = "Collection Date",  # x-axis label
      y = "Frequency",  # y-axis label
      color = "Dataset Source"  # legend title
    ) +
    theme_minimal() +  # Apply a clean minimalistic theme 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Format x-axis and position legend
}

# Function to plot frequency trajectories of selected variants ----
# This function filters the dataset to only include the specified variants 
# and generates a line plot showing how their frequencies changed over time.
plot_variant_frequency_trajectories <- function(variants) {
  
  # Filter the dataset to only include the selected variants
  filtered_data <- sanger_analysis_data %>%
    filter(variant %in% variants)
  
  # Generate plot
  ggplot(filtered_data, aes(x = collection_date, y = variant_frequency, color = variant)) +
    geom_line(linewidth = 1.2) +  # Plot frequency trends as lines
    geom_point(size = 2, alpha = 0.7) +  # Add points for better visibility of trends
    scale_color_manual(values = variant_colours) +  # Use predefined color scheme for variants
    scale_x_date(
      date_breaks = "2 month",  # Add tick marks every 2 months
      date_minor_breaks = "1 month",  # Add minor tick marks every month
      date_labels = "%Y-%m"  # Format date labels as YYYY-MM
    ) +
    labs(
      title = paste("Frequency Trajectories of", paste(variants, collapse = ", ")),  # Dynamic title based on selected variants
      x = "Collection Date",  # X-axis label
      y = "Frequency",  # Y-axis label
      color = "Variant"  # Legend title
    ) +
    theme_minimal() +  # apply a clean minimalistic theme
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5),  # Make x-axis labels diagonal for readability
          legend.position = "top",  # Move legend to the top for better visibility
          legend.title = element_text(face = "bold"),  # Format legend title for emphasis
          plot.title = element_text(hjust = 0.5, face = "bold"))  # Center and bold the title
}

# Function to plot logistic growth curves for multiple variants
# This function visualizes logistic growth curves for selected variants using a faceted layout.
# it fits a logistic model to the observed variant frequencies and overlays the predicted growth trajectory.
plot_logistic_growth <- function(data, growth_phases, variants) {
  
  # Create lists to store observed and predicted data for each variant
  observed_list <- list()
  predicted_list <- list()
  growth_rates <- list()  # Store growth rates for annotation purposes
  
  for (variant in variants) {
    
    # Extract growth phase start and end dates for the selected variant
    phase_dates <- growth_phases %>%
      filter(variant == !!variant)  # Filter for the specific variant
    
    # Subset the data to include only the relevant growth phase
    variant_growth_phase <- data %>%
      filter(variant == !!variant, 
             collection_date >= phase_dates$start_date, 
             collection_date <= phase_dates$end_date) %>%
      mutate(variant = variant)  # Ensure variant column is retained for faceting
    
    # Fit a logistic growth model to the observed frequency data
    logistic_results <- fit_logistic_growth_general(
      variant_growth_phase, 
      time_col = "collection_date", 
      frequency_col = "variant_frequency", 
      group_col = "variant"
    )
    
    # Extract the estimated growth rate (`s`)
    growth_rate <- unique(logistic_results$s)
    
    # Store the growth rate for annotation on the plot
    growth_rates[[variant]] <- data.frame(
      variant = variant,
      s = round(growth_rate, 4),  # round for readability
      x_pos = max(variant_growth_phase$collection_date),  # position at the rightmost date
      y_pos = min(variant_growth_phase$variant_frequency, na.rm = TRUE) * 0.5  # position near the bottom
    )
    
    # Store the observed and predicted data
    observed_list[[variant]] <- variant_growth_phase
    predicted_list[[variant]] <- logistic_results
  }
  
  # Combine observed and predicted data for all variants
  observed_data <- bind_rows(observed_list)
  predicted_data <- bind_rows(predicted_list)
  growth_rate_data <- bind_rows(growth_rates)  # Convert the growth rate list to a dataframe
  
  # Generate the logistic growth plot with facets for each variant
  ggplot(observed_data, aes(x = collection_date, y = variant_frequency, color = variant)) +
    geom_point(size = 2, alpha = 0.7) +  # Plot observed frequencies as points
    geom_line(data = predicted_data, aes(x = collection_date, y = predicted_frequency, color = variant), size = 1) +  
    geom_text(data = growth_rate_data, aes(x = x_pos, y = y_pos, label = paste0("s = ", s)), 
              color = "black", size = 3.5, fontface = "bold", hjust = 1, vjust = 0) +  # Add growth rate annotation
    scale_color_manual(values = variant_colours) +  # Apply custom color palette
    scale_x_date(date_labels = "%Y-%m") +  # Format x-axis labels as year-month
    labs(
      title = "Logistic Growth Fit for SARS-CoV-2 Variants",
      x = "Collection Date",
      y = "Frequency",
      color = "Variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Ensure readable x-axis labels
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    facet_grid(. ~ variant, scales = "free_x")  # Create side-by-side facets for different variants
}

# Function to plot delta variant frequency across different regions over time
# This function generates a line plot showing how the frequency of delta varied across different regions in England.
plot_delta_frequencies <- function(data) {
  ggplot(data, aes(x = week, y = delta_frequency, color = phecname)) +  # X-axis represents the week, y-axis represents frequency
    geom_line(size = 1.2) +  # Plot delta frequency as lines
    geom_point(size = 2, alpha = 0.7) +  # Add points for better visibility
    scale_color_manual(values = region_colours) +  # Apply predefined color scheme for regions
    labs(
      title = "Delta Variant Frequency Over Time (Weekly Smoothed)",  # plot title
      x = "Collection Date",  # x-axis label
      y = "Frequency of Delta",  # y-axis label
      color = "Region"  # legend title
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Format x-axis labels for readability
}
# Function to plot the logistic growth of the Delta variant by region ----
# This function plots the logistic growth of the Delta variant's frequency over time for different regions.

plot_logistic_growth_region <- function(observed_data, predicted_data) {
  
  # Extract growth rate (s) and initial frequency (f0) to annotate each region
  growth_info <- predicted_data %>%
    group_by(phecname) %>%   # Group data by region (phecname)
    slice_max(week, n = 1) %>%  # Select the latest time point (most recent week)
    mutate(
      label = paste0("s = ", round(s, 3), "\nf0 = ", round(f0, 3)),  # Format text annotation for s and f0
      x_pos = week,  # Position x-axis label at the most recent week
      y_pos = 0.05   # Fixed y-axis position for annotation (near the bottom of plot)
    )
  
  # Generate the logistic growth plot
  ggplot(observed_data, aes(x = week, y = delta_frequency, color = phecname)) +
    geom_point(size = 1.8, alpha = 0.7) +  # Scatter plot: observed frequencies of Delta variant
    geom_line(data = predicted_data, aes(x = week, y = predicted_frequency, group = phecname, color = phecname),
              size = 1) +  # Line plot: predicted logistic growth curves
    geom_text(data = growth_info,  
              aes(x = x_pos, y = y_pos, label = label),  
              color = "black", size = 3, fontface = "bold", hjust = 1, vjust = 0) +  # Text annotation for s and f0
    scale_color_viridis_d() +  # Use a colorblind-friendly palette
    labs(
      title = "Logistic Growth Fit for Delta Variant Frequency by Region",
      x = "Collection Date",
      y = "Frequency of Delta",
      color = "Region"
    ) +
    theme_minimal() +  # Apply a clean, minimalistic theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
    facet_wrap(~phecname, scales = "free_x")  # Create separate plots per region with independent x-axis scales
}


# Function to plot estimated vs sequenced Delta cases ----
# This function compares two different ways of tracking Delta variant cases over time
plot_estimated_vs_sequenced_delta_cases <- function(delta_estimates, sanger_data) {
  
  # Determine the first and last date where estimated Delta cases are non-zero
  min_date <- min(delta_estimates$date[delta_estimates$estimated_delta_cases > 0], na.rm = TRUE)
  max_date <- max(delta_estimates$date[delta_estimates$estimated_delta_cases > 0], na.rm = TRUE)
  
  # Prepare dataset for plotting: Estimated Daily Delta Cases
  plot_data <- delta_estimates %>%
    filter(date >= min_date & date <= max_date) %>%  # Filter only relevant date range
    select(date, estimated_delta_cases) %>%  # Keep only the date and case estimate columns
    rename(value = estimated_delta_cases) %>%  # Rename column for consistency
    mutate(source = "Estimated Daily Delta Cases")  # Add a new column to label data source
  
  # Prepare dataset for plotting: Sanger Weekly Delta Sequences
  sanger_weekly <- sanger_data %>%
    filter(variant == "B.1.617.2", collection_date >= min_date & collection_date <= max_date) %>%  # Filter Delta variant only
    select(collection_date, variant_count) %>%  # Keep only date and variant count
    rename(value = variant_count, date = collection_date) %>%  # Rename columns for consistency
    mutate(source = "Sanger Weekly Delta Sequences")  # Add a column to label the data source
  
  # Combine both datasets for visualization
  combined_plot_data <- bind_rows(plot_data, sanger_weekly)  # Merge estimated and sequenced data into a single dataframe
  
  # Generate comparative plot
  ggplot(combined_plot_data, aes(x = date, y = value, color = source)) +
    geom_line(size = 1.1) +  # Line plot to show trend over time
    geom_point(size = 1.5, alpha = 0.6) +  # Add points for better visualization
    scale_color_manual(values = c("Estimated Daily Delta Cases" = "#1F9E89",  
                                  "Sanger Weekly Delta Sequences" = "#433E85")) +  # Define custom colors for clarity
    labs(
      title = "Estimated Daily Delta Cases vs. Weekly Sanger Sequences",
      x = "Collection Date",
      y = "Number of Cases",
      color = "Source"
    ) +
    scale_x_date(limits = c(min_date, max_date), date_labels = "%Y-%m") +  # Restrict x-axis to the range of Delta cases
    theme_minimal() +  # Apply a minimalistic theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
          legend.position = "top",  # Move legend to the top for clarity
          legend.title = element_text(face = "bold"),  # Make legend title bold
          plot.title = element_text(hjust = 0.5, face = "bold"))  # Center and bold the plot title
}


# Function to plot estimated reproduction number ($R_t$) over time
# This function plots the estimated time-varying reproduction number ($R_t$)
# for delta variant cases, including a shaded confidence interval.
plot_rt_estimates <- function(rt_estimates) {
  
  # Rename columns for readability
  rt_estimates <- rt_estimates %>%
    rename(
      t_mid = t_start,
      Mean_R = `Mean(R)`,
      Std_R = `Std(R)`,
      Lower_CI = `Quantile.0.025(R)`,
      Upper_CI = `Quantile.0.975(R)`
    ) %>%
    mutate(t_mid = as.Date((t_mid + t_end) / 2, origin = "1970-01-01"))  # ensure t_mid is Date type
  
  # Create the plot
  ggplot(rt_estimates, aes(x = t_mid, y = Mean_R)) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "#433E85", alpha = 0.5) +  # confidence interval
    geom_line(color = "#433E85", size = 1) +  # Rt estimate line
    geom_hline(yintercept = 1, linetype = "dashed") +  # dashed line at Rt = 1
    scale_x_date(date_labels = "%b", date_breaks = "1 month") +  # format x-axis by month
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      title = expression("Time-varying reproduction number" ~ (R[t]) ~ "for Delta"),
      x = "Collection Date",
      y = expression("Reproduction number" ~ (R[t]))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}



  