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

# Function to annotate the mod sites
annotate_modification_sites <- function(file_paths, 
                                        config_file_path, 
                                        gene_transcript_table){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
  temp_seq_id <- config$biomart_search_filter
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Load symbols and ids table
  temp_symbol_df = loadRData(gene_transcript_table)
  
  # Add Ensembl gene, entrez, and hgnc symbol to the modification table
  temp_mod_ids = temp_df$ref_seq_id
  
  temp_df$hgnc_symbol = temp_symbol_df$hgnc_symbol[match(temp_df$ref_seq_id, temp_symbol_df[[temp_seq_id]])]
  
  temp_df$entrezgene_id = temp_symbol_df$entrezgene_id[match(temp_df$ref_seq_id, temp_symbol_df[[temp_seq_id]])]
  
  temp_df$ensembl_gene_id = temp_symbol_df$ensembl_gene_id[match(temp_df$ref_seq_id, temp_symbol_df[[temp_seq_id]])]
  
  # Set the first condition in effect_order is the reference (e.g., "control")
  reference_group <- effect_order[1]
  
  temp_diff_mod_cols = colnames(temp_df)[grep("diff_mod_rate_", colnames(temp_df))]
  
  temp_group_list = list()
  
  for (i in 2:length(effect_order)){
    
    ## inside the loop
    temp_treat_group <- effect_order[i]
    
    # Define the new column name for the differential calculation
    diff_col_name <- paste0("diff_mod_rate_", reference_group, "_to_", temp_treat_group)
    
    # Genes with higher modification rate in the reference condition
    temp_ref_genes <- unique(str_sort(temp_df$entrezgene_id[temp_df[[diff_col_name]] > 0], numeric = TRUE))
    
    # Genes with higher modification rate in the treatment condition
    temp_trt_genes <- unique(str_sort(temp_df$entrezgene_id[temp_df[[diff_col_name]] < 0], numeric = TRUE))
    
    temp_list <- list(reference_group = temp_ref_genes, temp_treat_group = temp_trt_genes)
    
    names(temp_list) <- c(reference_group, temp_treat_group)
    
    temp_ind = i - 1
    
    temp_group_list[[temp_ind]] <- temp_list
    
    names(temp_group_list)[temp_ind] <- paste0("GO_comparison_", reference_group, "_to_", temp_treat_group)
  }
  
  temp_group_list[["mod_data"]] <- temp_df
  
  return(temp_group_list)
}


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
gene_transcript_table <- args[3] 
output_file <- args[4]

# Get fasta subset and gene transcript table
annotated_mod_sites <- annotate_modification_sites(file_paths, config_file_path, gene_transcript_table)

annotated_mod_table <- annotated_mod_sites[[2]]

# Save the modified data frame with signature motifs
save(annotated_mod_table, file = output_file)
