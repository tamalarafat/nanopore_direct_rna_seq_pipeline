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
cran_packages <- c("ggplot2", "ggrepel", "ggseqlogo", "reshape2")
lapply(cran_packages, install_if_missing)

# Install required Bioconductor packages
bioc_packages <- c("enrichplot", "org.Hs.eg.db", "clusterProfiler")
lapply(bioc_packages, install_bioc_if_missing)

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
  
  # Load all entrez ids of the hg
  hg_entrz = as.character(str_sort(temp_symbol_df$entrezgene_id, numeric = TRUE))
  
  hg_entrz = unique(hg_entrz[!is.na(hg_entrz)])
  
  temp_group_list[["entrezgene_id"]] <- hg_entrz
  
  return(temp_group_list)
}

# Function to perform GSEA analysis
perform_gsea <- function(entrezID_list,
                         hg_gene_set_entrezID,
                         modification_table,
                         store_dir){
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/")
  
  # Get the item names of the list to generate folder with the name
  temp_item_names = names(entrezID_list)
  
  for (i in 1:length(temp_item_names)){
    
    # Let's create a directory to store the GO annotation files
    if (!dir.exists(str_c(temp_dir, temp_item_names[i]))){
      dir.create(str_c(temp_dir, temp_item_names[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_grp_dir <- str_c(temp_dir, temp_item_names[i], "/")
    
    gene_list <- entrezID_list[[temp_item_names[i]]]
    
    temp_conditions <- names(gene_list)
    
    for (j in 1:length(temp_conditions)) {
      
      # Let's create a directory to store the GO annotation files
      if (!dir.exists(str_c(temp_grp_dir, temp_conditions[j]))){
        dir.create(str_c(temp_grp_dir, temp_conditions[j]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      }
      
      # Assigning the directory path to a local variable
      temp_go_dir = str_c(temp_grp_dir, temp_conditions[j], "/")
      
      # Let's get the gene ids from the list
      temp_marker_genes = gene_list[[j]]
      
      # Store transcript, gene id table for the condition
      temp_match_ind = sort(match(temp_marker_genes, modification_table$entrezgene_id))
      temp_match_df = modification_table[temp_match_ind, c("hgnc_symbol", "ref_seq_id", "ensembl_gene_id", "entrezgene_id")]
      
      # All GO terms (significant, non-significant)
      write.csv(temp_match_df, file = str_c(temp_go_dir, temp_conditions[j], "_condition_gene_transcript_ids.csv"))
      
      # performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
      GeneSet_ORA_BP <- enrichGO(gene = temp_marker_genes,
                                 universe = hg_gene_set_entrezID,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "ENTREZID", 
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE,
                                 pool = FALSE)
      
      tryCatch({
        # Lets save the GO annotation results for BP class
        GeneSet_all_BP = GeneSet_ORA_BP@result
        
        # All GO terms (significant, non-significant)
        write.csv(GeneSet_all_BP, file = str_c(temp_go_dir, temp_conditions[j], "_condition_genes_all_biological_processes.csv"))
        
      }, 
      error = function(e){str_c("No biological processe annotation was found for ", temp_conditions[j], " genes.")}
      )
      
      # Figure 1 :: Generate dotplot for all GO terms - BP
      if (dim(GeneSet_ORA_BP)[1] < 1){
        print("No module found")
      } 
      else {
        p1 <- dotplot(GeneSet_ORA_BP, showCategory = 10) + 
          labs(colour = "p.adjust") + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black"),
                axis.title.x = element_text(size = 22, face = "bold", color = "black"), 
                axis.ticks.length = unit(.30, "cm"), 
                axis.text.x = element_text(size = 22, color = "black", face = "bold"),
                axis.text.y = element_text(size = 22, color = "black", face = "bold"),
                legend.key = element_rect(size = 22),
                legend.text = element_text(size = 22),
                legend.title = element_text(size = 22, face = "bold"),
                legend.spacing = unit(2.0, 'cm'),
                legend.key.size = unit(3,"line"),
                legend.position = "none")
        
        ggsave(filename = str_c(temp_go_dir, "GO_bp_",temp_conditions[j], "_condition_genes.png"), plot = p1, width = 14, height = 14, dpi = 300)
      }
      
      
      # Figure 2 :: Generate network of GO terms using all GO terms - BP
      if (dim(GeneSet_ORA_BP)[1] > 0){
        # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. This makes it easier to identify functional modules.
        GeneSet_ORA_BP <- pairwise_termsim(GeneSet_ORA_BP, method = "JC")
        print(dim(GeneSet_ORA_BP@termsim))
        
        if (dim(GeneSet_ORA_BP@termsim)[1] == 1){
          print("No module found")
        } 
        else {
          tryCatch({
            p2 <- emapplot(GeneSet_ORA_BP, color = "qvalue")  + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.title = element_blank(), 
                    text = element_text(size = 24, face = "bold"),
                    axis.line = element_blank(),
                    axis.ticks.length = unit(0, "cm"), 
                    axis.text = element_blank(),
                    legend.key = element_rect(size = 22),
                    legend.spacing = unit(2.0, 'cm'),
                    legend.key.size = unit(3,"line"),
                    legend.position = "none")
            ggsave(filename = str_c(temp_go_dir, "GO_bp_network",temp_conditions[j], "_condition_genes.png"), plot = p2, width = 14, height = 14, dpi = 300)
          }, 
          error = function(e) {
            message <- paste("Error in plotting network for GO terms - BP:", conditionMessage(e))
            cat(message, file = str_c(temp_go_dir, "error_log.txt"), append = TRUE)
            print(message)
          }
          )
        }
      }
      
      else {
        print("No module found")
      }
    }
  }
}

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file_paths <- args[1]  
config_file_path <- args[2] 
gene_transcript_table <- args[3] 
output_dir <- args[4]

annotated_mod_table <- annotate_modification_sites(file_paths, config_file_path, gene_transcript_table)

gene_set <- annotated_mod_table[1]

mod_table <- annotated_mod_table[[2]]

hg_entrez_universe <- annotated_mod_table[[3]]

# Get fasta subset and gene transcript table
perform_gsea(gene_set, hg_entrez_universe, mod_table, output_dir)

