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

# Custom color palette for major SARS-CoV-2 variants ----
lineage_colors <- c(
  "B.1.1.7" = "#56B4E9",  # Sky Blue
  "B.1.617.2" = "#E69F00", # Orange
  "BA.1" = "#009E73",  # Green
  "BA.2" = "#D55E00",  # Vermilion (Red-Orange)
  "BA.2.75" = "#CC79A7",  # Pink-Purple
  "BA.4" = "#0072B2",  # Deep Blue
  "BA.5" = "#F0E442",  # Yellow
  "BA.5.3" = "#882255", # Dark Magenta
  "XBB" = "#999999",  # Medium Gray
  "Other" = "#666666"  # Dark Gray
)


# Function to plot stacked area plot of total variant counts ----
plot_total_lineage_counts <- function(data) {
  ggplot(data, aes(x = collection_date, y = total_count, fill = major_lineage)) +
    geom_area(position = "stack") +  # Ensure proper stacking of counts
    scale_fill_manual(values = lineage_colors) +  # Apply colorblind-friendly palette
    labs(
      title = "Total Counts of Major Lineages Over Time",
      x = "Collection Date",
      y = "Total Count",
      fill = "Lineage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# Function to plot stacked area plot of variant frequencies ----
plot_lineage_frequencies <- function(data) {
  ggplot(data, aes(x = collection_date, y = lineage_frequency, fill = major_lineage)) +
    geom_area(position = "fill") +  # Stacked area plot (normalized)
    scale_fill_manual(values = lineage_colors) +  # Apply custom colors
    labs(
      title = "Frequency of Major Variants Over Time",
      x = "Collection Date",
      y = "Proportion",
      fill = "Lineage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#Function
plot_ba2_frequency_comparison <- function(ba2_cog, ba2_onscis) {
  ggplot() +
    # COG-UK Data (Weekly)
    geom_line(data = ba2_cog, 
              aes(x = collection_date, y = lineage_frequency, color = source), 
              linewidth = 1) +
    
    # ONS-CIS Data (10-day Binned)
    geom_line(data = ba2_onscis, 
              aes(x = collection_date_bin, y = lineage_frequency, color = source), 
              linewidth = 1) + 
    
    # Labels & Theme
    scale_color_manual(values = c("COG-UK (Weekly)" = "#D55E00", 
                                  "ONS-CIS (10-day Binned)" = "#CC6677")) + 
    labs(
      title = "BA.2 Frequency Trajectory (COG-UK vs ONS-CIS)",
      x = "Collection Date",
      y = "Frequency",
      color = "Dataset Source"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ..
plot_logistic_growth_curve <- function(growth_results, variant, color) {
  
  # Extract growth phase data and fitted model
  data <- growth_results$growth_data
  fitted_model <- growth_results$fitted_model
  growth_rate <- growth_results$growth_rate
  
  # Generate smooth sequence of dates
  smooth_dates <- seq(min(data$collection_date), max(data$collection_date), by = "1 day")
  
  # Compute predicted frequencies
  smooth_predictions <- data.frame(
    collection_date = smooth_dates,
    predicted_frequency = logistic_growth(
      as.numeric(smooth_dates - min(data$collection_date)),
      coef(fitted_model)["s"], 
      coef(fitted_model)["f0"]
    )
  )
  
  # Create the plot
  ggplot(data, aes(x = collection_date)) +
    geom_point(aes(y = lineage_frequency), color = "black", size = 2, alpha = 0.7) + # Actual data points
    geom_line(data = smooth_predictions, aes(x = collection_date, y = predicted_frequency), 
              color = color, size = 1) +  # Fitted logistic curve
    annotate(
      "text", 
      x = min(data$collection_date) + 30, 
      y = 0.8, 
      label = paste0("s = ", round(growth_rate, 4)), 
      color = color, 
      size = 5
    ) +
    labs(
      title = paste("Logistic Growth Fit for", variant),
      x = "Collection Date",
      y = "Frequency"
    ) +
    theme_minimal()
}


#Function to ...
plot_delta_frequencies <- function(data) {
  ggplot(data, aes(x = week, y = delta_frequency, color = phecname)) +  # Change 'date' -> 'week'
    geom_line(size = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_viridis_d() + # High contrast, colorblind-friendly
    labs(
      title = "Delta Variant Frequency Over Time (Weekly Smoothed)",
      x = "Collection Date",
      y = "Frequency of Delta",
      color = "Region"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



# Function to save plots as SVG ----
save_plot_svg <- function(data, filename, size, scaling, plot_function, ...) {
  size_inches <- size / 2.54  # Convert size from cm to inches
  svglite(filename, width = size_inches, height = size_inches, scaling = scaling)
  
  plot <- plot_function(data, ...)  # Call the respective plot function
  print(plot)  # Print the plot to save it
  
  dev.off()  # Close the device
}
