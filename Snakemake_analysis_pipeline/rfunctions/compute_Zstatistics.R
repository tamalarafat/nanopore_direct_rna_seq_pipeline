# Load necessary libraries
library(stringr)
library(yaml)

# Functions to perform the task

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create a function to perform Z-test
z_test_mod_sites <- function(df_object,
                             nmod_ctrl_id,
                             nmod_ki_id,
                             ncov_ctrl_id,
                             ncov_ki_id)
{
  
  
  nmod_ctrl_id = colnames(df_object)[grep(pattern = nmod_ctrl_id, colnames(df_object))]
  
  nmod_ki_id = colnames(df_object)[grep(pattern = nmod_ki_id, colnames(df_object))]
  
  ncov_ctrl_id = colnames(df_object)[grep(pattern = ncov_ctrl_id, colnames(df_object))]
  
  ncov_ki_id = colnames(df_object)[grep(pattern = ncov_ki_id, colnames(df_object))]
  
  m_list = list()
  
  for (i in c(1:nrow(df_object))){
    
    res <- prop.test(x = c(df_object[i, nmod_ctrl_id], df_object[i, nmod_ki_id]), n = c(df_object[i, ncov_ctrl_id], df_object[i, ncov_ki_id]))
    
    m_list[[i]] = data.frame(pVal_ztest = res$p.value, Xsquared = res$statistic, prop1 = res$estimate[[1]], prop2 = res$estimate[[2]]) 
    
    names(m_list)[i] <- i
  }
  
  temp_df = do.call(rbind, m_list)
  
  return(temp_df)
}

# Function to perform post processing on the merged modification table
compute_Zstatistics <- function(file_paths, config_file_path, filter_by_pval) {
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Assume the first condition in effect_order is the reference (e.g., "control")
  reference_group <- effect_order[1]
  
  # Loop through each condition after the reference
  for (i in 2:length(effect_order)) {
    # Define the groups to compare
    temp_nmod_ref <- paste0("avg_nmod_", reference_group)
    temp_nmod_eff <- paste0("avg_nmod_", effect_order[i])
    temp_ncov_ref <- paste0("avg_ncov_", reference_group)
    temp_ncov_eff <- paste0("avg_ncov_", effect_order[i])
    
    # Perform Z-test on the cleaned table
    temp_zstat_df = z_test_mod_sites(df_object = temp_df, 
                                     nmod_ctrl_id = temp_nmod_ref, 
                                     nmod_ki_id = temp_nmod_eff, 
                                     ncov_ctrl_id = temp_ncov_ref, 
                                     ncov_ki_id = temp_ncov_eff)
    
    # Subset the data to keep only the test stat
    temp_zstat_df = temp_zstat_df[, c(1, 2)]
    
    # Define the new column name for the differential calculation
    diff_col_name <- paste0("_", reference_group, "_to_", effect_order[i])
    
    # Assign column names
    colnames(temp_zstat_df) <- paste0(colnames(temp_zstat_df), diff_col_name)
    
    # Combine data
    temp_df = cbind(temp_df, temp_zstat_df)
  }
  
  # Ensure 'filter_by_pval' is provided
  if (!missing(filter_by_pval)) {
    
    if (filter_by_pval == "NULL"){
      return (temp_df)
    } 
    
    else {
      # Convert to numeric
      filter_by_pval = as.numeric(filter_by_pval)
      # Identify all columns with "diff_mod_rate" in the name
      temp_pval_col <- grep("pVal_ztest", colnames(temp_df), value = TRUE)
      
      # Proceed only if there is at least one 'diff_mod_rate' column
      if (length(temp_pval_col) > 0) {
        # Create a logical matrix indicating where abs(diff_mod_rate) >= threshold for each column
        temp_pval_pass <- abs(temp_df[, temp_pval_col, drop = FALSE]) < filter_by_pval
        
        # Determine rows to keep based on the number of columns
        if (length(temp_pval_col) > 1) {
          # Multiple columns: Retain rows where all columns meet the threshold
          rows_to_keep <- rowSums(temp_pval_pass) == length(temp_pval_col)
        } else {
          # Single column: Retain rows where the single column meets the threshold
          rows_to_keep <- temp_pval_pass[, 1]
        }
        
        # Apply the filtering
        temp_df <- temp_df[rows_to_keep, ]
      }
    }
  }
  
  # Set row names to default
  rownames(temp_df) <- NULL
  
  return(temp_df)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
filter_by_pval <- args[3]
output_file <- args[4]

# Run Z-statistics function
significant_mod_sites <- compute_Zstatistics(file_paths, config_file_path, filter_by_pval)

# Save the table with z-statistics data
save(significant_mod_sites, file = output_file)