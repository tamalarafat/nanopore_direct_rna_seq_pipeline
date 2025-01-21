# Load necessary libraries
library(stringr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
store_dir <- args[2]
store_folder <- args[3]
sample_name <- args[4]  # Capture the sample name

# Define the bed_conversion function
bed_conversion <- function(input_file, store_dir, store_folder, sample_name) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(file.path(store_dir, store_folder))) {
    dir.create(file.path(store_dir, store_folder), recursive = TRUE, mode = "0777")
  }
  
  # Column names for the .bed file
  col_names <- c(
    "ref_seq_id", "start_pos", "end_pos", "modified_base_code", "score",
    "strand", "start_pos_rep", "end_pos_rep", "color", "N_valid_cov",
    "percent_modified", "N_mod", "N_canonical", "N_other_modified",
    "N_delete", "N_fail", "N_diff", "N_nocall"
  )
  
  # Read and process the .bed file
  bed_read <- read.table(input_file)
  colnames(bed_read) <- col_names
  bed_read <- bed_read[, !colnames(bed_read) %in% c("score", "start_pos_rep", "end_pos_rep", "color")]
  
  # Modify column names by appending sample name
  colnames(bed_read) <- paste0(colnames(bed_read), "_", sample_name)
  
  # Define a vector for efficient replacement
  replacements <- c("ref_seq_id", "start_pos", "strand")
  for (i in c(1:length(replacements))) {
    colnames(bed_read)[grep(pattern = replacements[i], x = colnames(bed_read))] = replacements[i]
  }
  
  bed_read$ref_seq_id <- bed_read$ref_seq_id <- sub("\\..*$", "", bed_read$ref_seq_id)
  
  # Prepare output filename
  output_file_name <- str_replace(basename(input_file), "\\..*", "")
  save(bed_read, file = file.path(store_dir, store_folder, paste0(output_file_name, ".RData")))
}

# Run the function with provided arguments
bed_conversion(input_file, store_dir, store_folder, sample_name)
