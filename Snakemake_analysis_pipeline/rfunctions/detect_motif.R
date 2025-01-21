# Load necessary libraries
library(stringr)
library(yaml)
library(Biostrings)
library(biomaRt)
library(seqinr)

# Functions to perform the task

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to search motif in the fasta file
find_motif <- function(file_paths, config_file_path, fasta_path){
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get motif search parameters
  analyte_type <- config$analyte_type
  
  # Get the search filter for biomart
  motif_base <- config$motif_base
  motif_window_left <- as.numeric(config$motif_window_left) # can be of any length
  motif_window_right <- as.numeric(config$motif_window_right) # can be of any length
  motif_length <- motif_window_left + motif_window_right
  
  # Read the fasta file
  fasta_sub <- read.fasta(fasta_path)
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Create a column containing 5-mers start position
  temp_df$start_kmer = temp_df$end_pos - motif_window_left
  
  # Create a column containing 5-mers end position
  temp_df$end_kmer = temp_df$end_pos + motif_window_right
  
  # Create a column to store the 5-mers or DRACH motif
  temp_df$kmers = NA
  
  # For each transcript id find the 5-mers 7 DRACH motif
  for (i in c(1:nrow(temp_df))){
    
    if ((temp_df$start_kmer[i] >= 1) & (toupper(as.character(fasta_sub[temp_df$ref_seq_id[i]][[1]])[temp_df$start_kmer[i] + motif_window_left]) == motif_base)) {
      temp_df$kmers[i] = toupper(str_c(as.character(fasta_sub[temp_df$ref_seq_id[i]][[1]])[seq(temp_df$start_kmer[i], temp_df$start_kmer[i] + motif_length)], collapse = ""))
    }
    else {
      temp_df$kmers[i] = NA
    }
  }
  
  if (analyte_type == "RNA"){
    
    # Change the base T to U
    temp_df$kmers = gsub("T", "U", temp_df$kmers)
    
    # Rename a column by its name
    colnames(temp_df)[colnames(temp_df) == "kmers"] <- "motif_RNA"
  }
  
  else {
    # Rename a column by its name
    colnames(temp_df)[colnames(temp_df) == "kmers"] <- "motif_DNA"
  }
  
  return(temp_df)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
fasta_path <- args[3]
output_file <- args[4]

# Get fasta subset and gene transcript table
mod_sites_motifs = find_motif(file_paths, config_file_path, fasta_path)

# Save the gene transcript table
save(mod_sites_motifs, file = output_file)