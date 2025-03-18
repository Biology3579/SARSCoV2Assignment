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
  ggplot(data, aes(x = collection_date, y = total_count, fill = lineage)) +
    geom_area(position = "stack") +  # Ensure proper stacking of counts
    scale_fill_manual(values = lineage_colors) +  # Apply colorblind-friendly palette
    labs(
      title = "Total Counts of Major Lineages Over Time",
      x = "Collection Date",
      y = "Total Count",
      fill = "Variant"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# Function to plot stacked area plot of variant frequencies ----
plot_lineage_frequencies <- function(data) {
  ggplot(data, aes(x = collection_date, y = lineage_frequency, fill = lineage)) +
    geom_area(position = "fill") +  # Stacked area plot (normalized)
    scale_fill_manual(values = lineage_colors) +  # Apply custom colors
    labs(
      title = "Frequency of Major Variants Over Time",
      x = "Collection Date",
      y = "Proportion",
      fill = "Variant"
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
