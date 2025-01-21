# Load necessary libraries
library(stringr)
library(yaml)

# Functions to perform the task

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Find average of columns in the table
getAvrg <- function(df_object, 
                    id_columns, 
                    id_assign){
  
  df_object[id_assign] = rowMeans(df_object[, id_columns])
  
  return(df_object)
}

# Function to perform post processing on the merged modification table
post_processing_modified_sites <- function(file_paths, config_file_path, retain_Ncov_min, retain_min_diff_Modrate) {
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the conditions specified in configuration file - samples belonging to each condition
  conditions <- config$conditions
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Load the merged table
  df = loadRData(file_paths)
  
  # Sorting the data based on ID and start position
  df = df[order(df$ref_seq_id, df$start_pos), ]
  
  # Changing column name; Keeping only one column storing the end position
  colnames(df)[grep(pattern = "end_pos", colnames(df))][1] = "end_pos"
  
  # Removing column with repeated information on end position
  df = df[, !colnames(df) %in% colnames(df)[grep(pattern = "end_pos", colnames(df))][-1]]
  
  # Changing column name; Keeping only one column storing the base code
  colnames(df)[grep(pattern = "modified_base", colnames(df))][1] = "modified_base"
  
  # Removing column with repeated information on end position
  df = df[, !colnames(df) %in% colnames(df)[grep(pattern = "modified_base", colnames(df))][-1]]
  
  # Post filtering steps
  # Filtering step 4 - keeping only the sites that has a coverage of 30 or more in any of the samples (Nvalid_cov >= 30)
  df = df[rowSums(!df[, grep(pattern = "N_valid_cov", colnames(df))] < retain_Ncov_min)  != 0, ]
  
  # Filtering step 5 - removing sites that have 0% modification rate
  df = df[rowSums(df[, grep(pattern = "percent_modified", colnames(df))]) != 0, ]
  
  # Arrange the columns
  temp_cols = c("ref_seq_id", "start_pos", "end_pos", "strand", "modified_base", "N_valid_cov", "N_mod", "N_canonical", "percent_modified")
  
  # Keeping only those columns needed for the downstream analysis
  df = df[, colnames(df) %in% colnames(df)[grep(pattern = paste(temp_cols, collapse = "|"), colnames(df))]]
  
  # Ordering the columns to have sorted information
  temp_col_index = as.numeric(unlist(sapply(temp_cols, function (x){grep(x, colnames(df))})))
  
  # Order the table
  df = df[, temp_col_index]
  
  # Define columns for averaging with prefixes and base names
  temp_avg_cols <- c("avg_nmod_" = "N_mod_", "avg_ncov_" = "N_valid_cov_", "avg_mod_rate_" = "percent_modified_")
  
  # Calculate average of the columns - average of total coverage, number of modifications, and modification rate
  for (i in seq_along(temp_avg_cols)) {
    # Get prefix and base column name
    prefix <- names(temp_avg_cols)[i]
    base_col <- temp_avg_cols[i]
    
    # Apply averaging based on conditions
    for (condition in names(conditions)) {
      id_columns <- paste0(base_col, conditions[[condition]])  # Create column names for averaging
      assign_name <- paste0(prefix, condition)  # New column name with prefix and condition
      
      # Calculate and add the average to df
      df <- getAvrg(df, id_columns, assign_name)
    }
  }
  
  # # Calculate differential modification rate between the conditions
  # # Loop through pairs of groups in the effect order
  # for (i in 1:(length(effect_order) - 1)) {
  #   group1 <- paste0("avg_mod_rate_", effect_order[i])
  #   group2 <- paste0("avg_mod_rate_", effect_order[i + 1])
  #   
  #   # Define the new column name for the differential calculation
  #   diff_col_name <- paste0("diff_mod_rate_", effect_order[i], "_to_", effect_order[i + 1])
  #   
  #   # Calculate and store the differential modification rate
  #   df[[diff_col_name]] <- df[[group1]] - df[[group2]]
  # }
  
  # Assume the first condition in effect_order is the reference (e.g., "control")
  reference_group <- effect_order[1]
  
  # Loop through each condition after the reference
  for (i in 2:length(effect_order)) {
    # Define the groups to compare
    group1 <- paste0("avg_mod_rate_", reference_group)
    group2 <- paste0("avg_mod_rate_", effect_order[i])
    
    # Define the new column name for the differential calculation
    diff_col_name <- paste0("diff_mod_rate_", reference_group, "_to_", effect_order[i])
    
    # Calculate and store the differential modification rate
    df[[diff_col_name]] <- df[[group1]] - df[[group2]]
  }
  
  
  # Ensure 'retain_min_diff_Modrate' is provided
  if (!missing(retain_min_diff_Modrate)) {
    
    # Identify all columns with "diff_mod_rate" in the name
    temp_diff_modr_col <- grep("diff_mod_rate", colnames(df), value = TRUE)
    
    # Proceed only if there is at least one 'diff_mod_rate' column
    if (length(temp_diff_modr_col) > 0) {
      # Create a logical matrix indicating where abs(diff_mod_rate) >= threshold for each column
      diff_mod_rate_pass <- abs(df[, temp_diff_modr_col, drop = FALSE]) >= retain_min_diff_Modrate
      
      # Determine rows to keep based on the number of columns
      if (length(temp_diff_modr_col) > 1) {
        # Multiple columns: Retain rows where all columns meet the threshold
        rows_to_keep <- rowSums(diff_mod_rate_pass) == length(temp_diff_modr_col)
      } else {
        # Single column: Retain rows where the single column meets the threshold
        rows_to_keep <- diff_mod_rate_pass[, 1]
      }
      
      # Apply the filtering
      df <- df[rows_to_keep, ]
    }
  }
  
  # Set row names to default
  rownames(df) <- NULL
  
  return(df)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
retain_Ncov_min <- as.numeric(args[3])
retain_min_diff_Modrate <- as.numeric(args[4])
output_file <- args[5]

# Run post processing function
processed_mod_table <- post_processing_modified_sites(file_paths, config_file_path, retain_Ncov_min, retain_min_diff_Modrate)

# Save the final merged data
save(processed_mod_table, file = output_file)