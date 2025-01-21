# Load necessary libraries
library(stringr)

# Function - loads and store the RData file to a variable
loadRData <- function(fileName){
  
  #loads an RData file, and returns it
  load(fileName)
  
  get(ls()[ls() != "fileName"])
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
retain_strand <- args[2]
retain_Ndiff_percent <- args[3]
retain_Ndel_percent <- args[4]
append_data <- args[5]
store_dir <- args[6]
store_folder <- args[7]
output_file_name <- args[8]

# Function - to filter out sites
filter_modification_sites <- function(input_file,
                                      retain_strand = NULL,
                                      retain_Ndiff_percent = NULL,
                                      retain_Ndel_percent = NULL,
                                      append_data = FALSE,
                                      store_dir, 
                                      store_folder,
                                      output_file_name){
  
  # Load the RData object
  df = loadRData(input_file)
  
  # Filtering step - by strand
  if (!missing(retain_strand)){
    df = df[df$strand == retain_strand, ]
  }
  
  # Filtering step - by percentage of missing base or number of total delete
  if (!missing(retain_Ndiff_percent)){
    df$total_diff = df$N_valid_cov + df$N_diff
    
    df$Ndiff_percent = (df$N_diff/df$total_diff) * 100
    
    df = df[df$Ndiff_percent <= retain_Ndiff_percent, ]
    
    df = df[, !colnames(df) %in% "total_diff"]
    
    if (append_data == FALSE){
      df = df[, !colnames(df) %in% "Ndiff_percent"]
    }
  }

  # Filtering step - by percentage of missing base or number of total delete
  if (!missing(retain_Ndel_percent)){
    df$total_del = df$N_valid_cov + df$N_delete
    
    df$Ndel_percent = (df$N_delete/df$total_del) * 100
    
    df = df[df$Ndel_percent <= retain_Ndel_percent, ]
    
    df = df[, !colnames(df) %in% "total_del"]
    
    if (append_data == FALSE){
      df = df[, !colnames(df) %in% "Ndel_percent"]
    }
  }
  
  # Store the filtered data
  # Prepare output filename
  # output_file_name <- str_replace(basename(input_file), "\\..*", "")
  save(df, file = output_file_name)
}

# Run the function with provided arguments
filter_modification_sites(input_file, retain_strand, retain_Ndiff_percent, retain_Ndel_percent, append_data, store_dir, store_folder, output_file_name)
