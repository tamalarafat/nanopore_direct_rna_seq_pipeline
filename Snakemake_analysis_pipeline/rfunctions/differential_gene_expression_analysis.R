# Set CRAN mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Check for BiocManager and install if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Function to install CRAN packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = FALSE)  # Avoid installing dependencies unless necessary
  }
}

# Function to install Bioconductor packages without updating existing packages
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Install required CRAN packages
cran_packages <- c("ggplot2", "ggrepel", "ggseqlogo", "reshape2", "magrittr")
lapply(cran_packages, install_if_missing)

# Install required Bioconductor packages
bioc_packages <- c("enrichplot", "org.Hs.eg.db", "clusterProfiler", "tximport", "DESeq2")
lapply(bioc_packages, install_bioc_if_missing)

# Load necessary libraries
library(stringr)
library(readr)  # for type_convert()
library(tximport)
library(DESeq2)
library(magrittr)


# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

perform_deg_deseq <- function(input_file,
                              config_file_path){
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the conditions specified in configuration file - samples belonging to each condition
  conditions <- config$conditions
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Reference condition to which expression will be compared
  reference_group <- effect_order[1]
  
  # Get the cell line information
  cell_line <- config$cell_line
  
  # Controlling for cell line
  control_for_cell_line <- config$control_for_cell_line_effect
  
  # Gene filtering threshold - filter out genes if at least this amount of expression is detected in each of the samples
  gene_filter_threshold <- as.numeric(config$gene_exp_threshold_filter)
  
  # Load the gene count table
  temp_gct <- loadRData(input_file)
  
  # Get the sample names
  sample_names = colnames(temp_gct[["counts"]])
  
  # Get the cell line names
  sample_cell_line <- sapply(sample_names, function(sample) {
    name <- names(cell_line)[sapply(cell_line, function(x) sample %in% x)]
    if (length(name) > 0) name else NA  # Return NA if not found
  })
  
  # Get the condition names
  sample_condition <- sapply(sample_names, function(sample) {
    name <- names(conditions)[sapply(conditions, function(x) sample %in% x)]
    if (length(name) > 0) name else NA  # Return NA if not found
  })
  
  # Generate a sample table to track the coldata
  coldata <- data.frame(names = factor(sample_names), 
                        cell = factor(sample_cell_line),
                        condition = factor(sample_condition))
  
  # Assign rownames to the coldata
  rownames(coldata) <- coldata$names
  
  # Assign the reference group: the group against which comparison will be performed
  if (levels(coldata$condition)[2] != reference_group){
    coldata$condition %<>% relevel(levels(coldata$condition)[2])
  }
  
  
  # Initialize the DESeq object - controlling for cell line effect
  if (control_for_cell_line == TRUE){
    dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                    colData = coldata, 
                                    design = ~ cell + condition)
  }
  # Initialize the DESeq object - without controlling for cell line effect
  else {
    dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                    colData = coldata, 
                                    design = ~ condition)
  }
  
  
  # Filter the data
  # criteria 1: remove rows containing only 0's
  temp_cr_1 = (!rowSums(counts(dds) != 0) == 0)
  
  # criteria 2:remove rows if the gene expression is less than 10 (default) or user-specified in "filter_exp_across_n_sample" args
  temp_cr_2 = rowSums(counts(dds) >= gene_filter_threshold) == ncol(counts(dds))
  
  # combine the two criteria
  keep = (temp_cr_1 & temp_cr_2)
  
  # Subset the deseq object to keep only those genes that satisfy the filtering criteria
  dds <- dds[keep, ]
  
  # Estimate the size factor of the dds object
  dds <- estimateSizeFactors(dds)
  
  # Running the differential expression pipeline
  dds <- DESeq(dds)
  
  # Extract and store the differentially expressed genes result in a variable
  res <- results(dds)
  
  # Create an list to store the DESeq object and DEG results file
  temp_list = list(res = res, dds = dds)
  
  return(temp_list)
}

# Function to extract differential gene expression table
extract_DEG_table <- function(res_DEseq){
  
  # Create a dataframe with the genes that passed the filtering criteria
  temp_df = data.frame(ID = rownames(res_DEseq), log2FC = res_DEseq$log2FoldChange, pvalue = res_DEseq$pvalue, padj = res_DEseq$padj, lfcSE = res_DEseq$lfcSE, stat = res_DEseq$stat)
  
  # Assign rownames to gene IDs
  rownames(temp_df) <- temp_df$ID
  
  # Assign group or condition tags to the foldchange direction
  
  # Get group/condition information
  temp_string <- mcols(res_DEseq)[2, 2]
  
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
  
  return(temp_df)
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2]
output_deseq_object <- args[3]
output_deseq_result <- args[4]
output_deg_table <- args[5]

dds_res <- perform_deg_deseq(file_paths, config_file_path)

# get and store deseq2 object
deseq_object <- dds_res$dds
save(deseq_object, file = output_deseq_object)

# get and store deseq2 result
deseq_res <- dds_res$res
save(deseq_res, file = output_deseq_result)

# get and store deg table
deg_table <- extract_DEG_table(deseq_res)
# write the table in tsv format
write.table(deg_table, file = output_deg_table, row.names = FALSE, sep = "\t")
