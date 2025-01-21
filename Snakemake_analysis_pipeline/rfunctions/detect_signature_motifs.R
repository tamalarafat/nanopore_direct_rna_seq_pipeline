# Load necessary libraries
library(stringr)
library(yaml)
library(Biostrings)
library(biomaRt)
library(seqinr)
library(dplyr)  # For mutating columns
library(readr)  # for type_convert()

# Functions to perform the task

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to generate fasta subset and generate gene transcript table
detect_signature_motifs <- function(file_paths, config_file_path){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Get motif search parameters
  analyte_type <- config$analyte_type
  
  # Get the motif pattern from the cofig file
  signature_motif_pattern = config$signature_motif_pattern
  
  # Motif column name
  temp_motif_col = paste0("motif_", analyte_type)
  
  # Loop through each signature motif and assign "YES"/"NO" based on matching rows
  for (i in seq_along(signature_motif_pattern)) {
    # Get the signature name and motif pattern
    signature_name <- names(signature_motif_pattern)[i]
    signature_pattern <- paste0("^", signature_motif_pattern[[i]], "$")
    
    # Create a logical vector indicating matches
    matching_rows <- grepl(signature_pattern, temp_df[[temp_motif_col]])
    
    # Directly assign "YES" and "NO" using ifelse (vectorized)
    temp_df[[signature_name]] <- ifelse(matching_rows, "YES", "NO")
  }
  
  return(temp_df)
}


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
output_file <- args[3]

# Process the input file to detect signature motifs
mod_sites_signature_motifs <- detect_signature_motifs(file_paths, config_file_path)

# Save the modified data frame with signature motifs
save(mod_sites_signature_motifs, file = output_file)

# Load conditions and effect direction from config.yaml
config <- yaml::read_yaml(config_file_path)

# Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
signature_motif_pattern <- config$signature_motif_pattern

# Save separate tables for each signature motif with "YES" in the respective column
for (i in seq_along(signature_motif_pattern)) {
  # Get the signature name
  signature_name <- names(signature_motif_pattern)[i]
  
  # Filter rows where the signature motif is marked as "YES"
  temp_df <- mod_sites_signature_motifs[mod_sites_signature_motifs[[signature_name]] == "YES", ]
  
  # Save the filtered data frame for the current signature
  save(temp_df, file = file.path(dirname(output_file), paste0(signature_name, "_significant_mod_sites.RData")))
}