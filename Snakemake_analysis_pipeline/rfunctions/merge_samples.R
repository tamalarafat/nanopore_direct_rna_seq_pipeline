# Load necessary libraries
library(stringr)

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Recursive merging function for a list of data frames
merge_all_samples <- function(file_paths, merge_by) {
  # Load the first file as the starting data frame
  merged_data <- loadRData(file_paths[1])
  
  # Loop through remaining files and merge iteratively
  for (i in 2:length(file_paths)) {
    next_data <- loadRData(file_paths[i])
    merged_data <- merge(merged_data, next_data, by = merge_by, all = FALSE)  # Only common sites
  }
  
  return(merged_data)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- unlist(strsplit(args[1], " "))  # Convert input list to a vector of paths
merge_columns <- unlist(strsplit(args[2], ","))  # Columns to merge on
output_file <- args[3]

# Run merging function
merged_data <- merge_all_samples(file_paths, merge_columns)

# Save the final merged data
save(merged_data, file = output_file)
