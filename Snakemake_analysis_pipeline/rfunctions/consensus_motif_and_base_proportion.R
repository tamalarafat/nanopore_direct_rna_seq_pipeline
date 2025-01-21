# Load necessary libraries
library(stringr)
library(yaml)
library(Biostrings)
library(biomaRt)
library(seqinr)
library(dplyr)  # For mutating columns
library(readr)  # for type_convert()
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggseqlogo)
library(reshape2)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Function to load .RData files and return the loaded data
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to annotate the mod sites
motif_data_per_condition <- function(file_paths, 
                                     config_file_path){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
  temp_analyte_id <- config$analyte_type
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Set the first condition in effect_order is the reference (e.g., "control")
  reference_group <- effect_order[1]
  
  temp_diff_mod_cols = colnames(temp_df)[grep("diff_mod_rate_", colnames(temp_df))]
  
  temp_motif_cols = paste0("motif_", temp_analyte_id)
  
  temp_group_list = list()
  
  for (i in 2:length(effect_order)){
    
    ## inside the loop
    temp_treat_group <- effect_order[i]
    
    # Define the new column name for the differential calculation
    diff_col_name <- paste0("diff_mod_rate_", reference_group, "_to_", temp_treat_group)
    
    # Motif sequences in the reference condition
    temp_ref_df <- temp_df[temp_df[[diff_col_name]] > 0, ]
    
    # how many k-mers were detected
    temp_ref_df = temp_ref_df[!is.na(temp_ref_df[[temp_motif_cols]]), ]
    
    # Motif sequences in the treatment condition
    temp_trt_df <- temp_df[temp_df[[diff_col_name]] < 0, ]
    
    # how many k-mers were detected
    temp_trt_df = temp_trt_df[!is.na(temp_trt_df[[temp_motif_cols]]), ]
    
    temp_list <- list(reference_group = temp_ref_df, temp_treat_group = temp_trt_df)
    
    names(temp_list) <- c(reference_group, temp_treat_group)
    
    temp_ind = i - 1
    
    temp_group_list[[temp_ind]] <- temp_list
    
    names(temp_group_list)[temp_ind] <- paste0("conditon_comparison_", reference_group, "_to_", temp_treat_group)
  }
  
  temp_group_list[["col_id"]] <- temp_motif_cols
  
  temp_group_list[["analyte_type"]] <- temp_analyte_id
  
  return(temp_group_list)
}

# Function to generate motif logo
generate_consensus_motif_logo <- function(data_list,
                                          motif_col_id,
                                          store_dir){
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/")
  
  # Get the item names of the list to generate folder with the name
  temp_item_names = names(data_list)
  
  for (i in 1:length(temp_item_names)){
    
    # Let's create a directory to store the GO annotation files
    if (!dir.exists(str_c(temp_dir, temp_item_names[i]))){
      dir.create(str_c(temp_dir, temp_item_names[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_grp_dir <- str_c(temp_dir, temp_item_names[i], "/")
    
    motif_data <- data_list[[i]]
    
    temp_conditions <- names(motif_data)
    
    for (j in 1:length(temp_conditions)) {
      
      # Let's create a directory to store the GO annotation files
      if (!dir.exists(str_c(temp_grp_dir, temp_conditions[j]))){
        dir.create(str_c(temp_grp_dir, temp_conditions[j]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      }
      
      # Assigning the directory path to a local variable
      temp_motif_dir = str_c(temp_grp_dir, temp_conditions[j], "/")
      
      # Let's get the gene ids from the list
      temp_df = motif_data[[j]]
      
      p1 <- ggseqlogo(temp_df[[motif_col_id]], method = "probability") +
        labs(title = paste(nrow(temp_df), "sites;", length(unique(temp_df$ref_seq_id)), 
                           "transcripts;", length(unique(temp_df$ensembl_gene_id)), "genes", sep = " ")) +
        theme_classic() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), # Background of the entire plot
          axis.line = element_blank(),
          axis.ticks.length = unit(0, "cm"), 
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 18),
          legend.position = "top",   # Move legend to the top
          legend.title = element_text(size = 18, color = "black"),
          legend.key.size = unit(2, "line"),
          legend.text = element_text(size = 18, color = "black")
        )
      ggsave(filename = str_c(temp_motif_dir, "consensus_motif_",temp_conditions[j], "_condition.png"), plot = p1, width = 12, height = 6, dpi = 300)
      
      # Save motif information for each condition
      temp_motif_data = data.frame(table(temp_df[[motif_col_id]]))
      temp_motif_data <- temp_motif_data[order(temp_motif_data$Freq, decreasing = TRUE), ]
      write.csv(temp_motif_data, file = str_c(temp_motif_dir, temp_conditions[j], "_condition_motifs.csv"))
    }
  }
}


# Function to generate proportion of bases at each site of themotif
motif_base_proportion_plot <- function(data_list,
                                       motif_col_id,
                                       analyte_type,
                                       store_dir){
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/")
  
  # Get the item names of the list to generate folder with the name
  temp_item_names = names(data_list)
  
  for (i in 1:length(temp_item_names)){
    
    # Let's create a directory to store the GO annotation files
    if (!dir.exists(str_c(temp_dir, temp_item_names[i]))){
      dir.create(str_c(temp_dir, temp_item_names[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_grp_dir <- str_c(temp_dir, temp_item_names[i], "/")
    
    motif_data <- data_list[[i]]
    
    temp_conditions <- names(motif_data)
    
    for (j in 1:length(temp_conditions)) {
      
      # Let's create a directory to store the GO annotation files
      if (!dir.exists(str_c(temp_grp_dir, temp_conditions[j]))){
        dir.create(str_c(temp_grp_dir, temp_conditions[j]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      }
      
      # Assigning the directory path to a local variable
      temp_motif_dir = str_c(temp_grp_dir, temp_conditions[j], "/")
      
      # Let's get the gene ids from the list
      temp_df = motif_data[[j]]
      
      # Your vector of sequences
      temp_motif_seq <-  temp_df[[motif_col_id]]
      
      # Convert the vector into a matrix where each column is a position in the sequences
      seq_matrix <- do.call(rbind, strsplit(temp_motif_seq, ""))
      
      if (analyte_type == "RNA"){
        # Calculate the proportion of each base at each position
        proportions <- apply(seq_matrix, 2, function(x) {
          table(factor(x, levels = c("A", "C", "G", "U"))) / length(x)
        })
      }
      else {
        # Calculate the proportion of each base at each position
        proportions <- apply(seq_matrix, 2, function(x) {
          table(factor(x, levels = c("A", "C", "G", "T"))) / length(x)
        })
      }
      
      
      # Convert the matrix of proportions to a data frame
      proportion_df <- as.data.frame(t(proportions))
      
      proportion_df$position = as.factor(rownames(proportion_df))
      
      # Transform the data into long format for plotting
      df_long = reshape2::melt(proportion_df, variable.name = "position", value.name = "proportion")
      
      # Rename columns to represent positions
      colnames(df_long)[2] <- "bases"
      
      # Define the scale_fill_manual layer based on the analyte_type condition
      fill_scale <- if (analyte_type == "RNA") {
        scale_fill_manual(values = c("A" = grp_col[2], "U" = grp_col[1], "G" = grp_col[3], "C" = grp_col[4]))
      } else {
        scale_fill_manual(values = c("A" = grp_col[2], "T" = grp_col[1], "G" = grp_col[3], "C" = grp_col[4]))
      }
      
      # Build the plot and add the conditional fill_scale
      p1 <- ggplot(data = df_long, aes(x = position, y = proportion)) + 
        geom_bar(stat = "identity", position = "fill", aes(fill = bases)) +
        xlab("Position") + ylab("Proportion of bases") + 
        fill_scale +  # Add the conditional fill scale here
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size = 22, color = "black"),
          axis.ticks.length = unit(0, "cm"), 
          axis.text = element_text(size = 22, colour = "black"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(1, "line"), 
          legend.text = element_text(size = 22)
        ) + 
        guides(color = guide_legend(override.aes = list(size = 4)))
      
      ggsave(filename = str_c(temp_motif_dir, "base_proportion_consensus_motif_",temp_conditions[j], "_condition.png"), plot = p1, width = 12, height = 6, dpi = 300)
    }
  }
}


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
output_dir <- args[3]


motif_data_list <- motif_data_per_condition(file_paths, config_file_path)

motif_data <- motif_data_list[1]

temp_motif_cols = motif_data_list[[2]]

analyte_type = motif_data_list[[3]]

# Get fasta subset and gene transcript table
generate_consensus_motif_logo(motif_data, temp_motif_cols, output_dir)

motif_base_proportion_plot(motif_data, temp_motif_cols, analyte_type, output_dir)

