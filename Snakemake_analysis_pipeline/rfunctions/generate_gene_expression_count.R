# Function to install Bioconductor packages without updating existing packages
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Install required Bioconductor packages
bioc_packages <- c("enrichplot", "org.Hs.eg.db", "clusterProfiler", "tximport")
lapply(bioc_packages, install_bioc_if_missing)

# Load necessary libraries
library(stringr)
library(yaml)
library(readr)  # for type_convert()
library(tximport)


# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to prepare and process nanocount output table
prepare_nanocount_data <- function(file_paths, sample_names){
  
  # Let's create a directory to store the GO annotation files
  if (!dir.exists("temporary_data")){
    dir.create("temporary_data", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = "temporary_data/"
  
  # Loop through remaining files and merge iteratively
  for (i in 1:length(file_paths)) {
    
    file_name = basename(file_paths[i])
    
    # Find matches using grep
    matching_sample <- sample_names[sapply(sample_names, function(x) grepl(x, file_name))]
    
    # read the file, in tsv format
    temp = read_tsv(file_paths[i], show_col_types = FALSE)
    
    # rearrange the columns to a desired format
    temp = temp[, c("transcript_name", "transcript_length", "transcript_length", "est_count", "tpm")]
    
    # assign column names to the count table
    colnames(temp) = c("target_id", "length", "eff_length", "est_counts", "tpm")
    
    temp$target_id <- sub(pattern = "\\..*", replacement = "", temp$target_id)
    
    # Sort the data frames using numeric extraction from target_id
    temp <- temp[order(as.numeric(sub("ENST", "", temp$target_id))), ]
    
    # write the table in tsv format
    write.table(temp, file = paste0(temp_dir, matching_sample, ".tsv"), row.names = FALSE, sep = "\t")
  }
}

# Function to generate gene count table
generate_gene_count_table <- function(input_file, config_file_path, gene_transcript_symbols){
  
  # List the files
  input_files = list.files(path = input_file, pattern = ".tsv", full.names = TRUE)
  names(input_files) <- sub(pattern = "\\..*", "", basename(input_files))
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  print(config)
  
  # Print and verify
  print(config[["assign_gene_ID_type"]])
  
  # Get the name of the gene ID conversion type from the configuration file
  conversion_ID_type <- config[["assign_gene_ID_type"]]
  
  # Load symbols and ids table
  temp_symbol_df = loadRData(gene_transcript_symbols)
  
  # Column to retreive from gene id symbol data
  temp_col_names <- c("ensembl_transcript_id", conversion_ID_type)
  
  tx2gene = temp_symbol_df[, temp_col_names]
  
  # Assign colnames to the id transform
  colnames(tx2gene) = c("TXNAME", "GENEID")
  
  tx2gene = tx2gene[!(is.na(tx2gene$GENEID) | tx2gene$GENEID == ""), ]
  
  tx2gene <- tx2gene[match(unique(tx2gene$TXNAME), tx2gene$TXNAME), ]
  
  # reading read count files to generate the gene-count table 
  txi <- tximport(input_files, type = "none", txOut = TRUE, txIdCol = "target_id", abundanceCol = "tpm",
                  countsCol = "est_counts", lengthCol = "length", importer = read_tsv)
  
  # Summarize the transcript abundance file
  txi.sum <- summarizeToGene(txi, tx2gene)
  
  return(txi.sum)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- unlist(strsplit(args[1], " "))  # Convert input list to a vector of paths
sample_names <- unlist(strsplit(args[2], " "))
config_file_path <- args[3]
gene_transcript_symbols <- args[4]
output_file <- args[5]

temp_paths <- "temporary_data" 

prepare_nanocount_data(file_paths, sample_names)

summarized_data <- generate_gene_count_table(temp_paths, config_file_path, gene_transcript_symbols)

save(summarized_data, file = output_file)

unlink(temp_paths, recursive = TRUE)
