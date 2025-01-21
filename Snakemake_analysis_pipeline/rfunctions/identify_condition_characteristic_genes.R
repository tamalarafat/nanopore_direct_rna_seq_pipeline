# Load necessary libraries
library(stringr)
library(readr)  # for type_convert()
library(tximport)
library(DESeq2)

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to identify characteristic DEGs for each condition from the differentially expressed genes identified using DESEQ2
identify_characteristic_degs <- function(input_file,
                                         config_file_path){
  
  
  res <- loadRData(file_paths)
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the filtering threshold information: by P-value
  pval_threshold <- config$p_val_threshold_deg
  
  # Get the filtering threshold information: by P-value
  fc_threshold <- config$log2FC_threshold_deg
  
  # Get the filtering threshold information: by adjusted P-value
  adj_pval_threshold <- config$adusted.pval_threshold_deg
  
  # Criteria 1: Identify top differentially expressed genes based on log2FC
  temp_fc = which(abs(res$log2FoldChange) > fc_threshold)
  
  # Criteria 2: Identify top differentially expressed genes based on pval
  if (pval_threshold != "FALSE"){
    pval_threshold <- as.numeric(pval_threshold)
    temp_pval = which(res$pvalue < pval_threshold)
    
    # Combine the filtering criteria
    temp_ind = union(temp_pval, temp_fc)
  } else {
    temp_ind = temp_fc
  }
  
  # Create a dataframe with the genes that passed the filtering criteria
  temp_df = data.frame(ID = rownames(res)[temp_ind], log2FC = res$log2FoldChange[temp_ind], pvalue = res$pvalue[temp_ind], padj = res$padj[temp_ind], lfcSE = res$lfcSE[temp_ind], stat = res$stat[temp_ind])
  
  # Assign rownames to gene IDs
  rownames(temp_df) <- temp_df$ID
  
  # Assign group or condition tags to the foldchange direction
  # Get group/condition information
  temp_string <- mcols(res)[2, 2]
  
  # Extract the strings before and after "vs"
  temp_matches <- regmatches(temp_string, regexec("([^ ]+) vs ([^ ]+)", temp_string))
  
  # Numerator or the reference group
  temp_numerator = temp_matches[[1]][2]
  
  # Denominator or the effect group where we want to observe differences
  temp_denominator = temp_matches[[1]][3]
  
  # Column name
  temp_col <- paste0("observed_expression_in_", temp_numerator)
  
  # Add group information to the data table
  temp_df[[temp_col]] = ifelse(temp_df$log2FC > 0, "UP", "DOWN")
  
  # Add group information to the data table
  temp_df$condition = ifelse(temp_df$log2FC > 0, temp_numerator, temp_denominator)
  
  # Rearrange the data table based on FC
  temp_df = temp_df[order(temp_df$log2FC, decreasing = TRUE), ]
  
  # Criteria 2: Identify top differentially expressed genes based on pval
  if (adj_pval_threshold != "FALSE"){
    adj_pval_threshold <- as.numeric(adj_pval_threshold)
    temp_adj_pval = which(temp_df$padj < adj_pval_threshold)
    temp_df = temp_df[temp_adj_pval, ]
  }
  
  # If to return by groups
  temp_list = list()
  
  # Add the numerator group
  temp_list[[temp_numerator]] = temp_df[temp_df$condition == temp_numerator, ]
  
  # Add the denominator group
  temp_list[[temp_denominator]] = temp_df[temp_df$condition == temp_denominator, ]
  
  # Add the denominator group
  temp_list[["condition_characteristic_deg"]] = temp_df
  
  return(temp_list)
} 

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2]
output_deg_table <- args[3]
temp_dir <- paste0(dirname(output_deg_table), "/")

characteristic_degs <- identify_characteristic_degs(file_paths, config_file_path)

characteristic_genes <- characteristic_degs[[3]]
write.table(characteristic_genes, file = output_deg_table, row.names = FALSE, sep = "\t")

write.table(characteristic_degs[[1]], file = paste(temp_dir, names(characteristic_degs)[1], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")

write.table(characteristic_degs[[2]], file = paste(temp_dir, names(characteristic_degs)[2], "_condition_characteristic_genes.tsv"), row.names = FALSE, sep = "\t")
