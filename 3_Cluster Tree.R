# Reference: https://mp.weixin.qq.com/s/def1EiaRf-eBrQrKOketkQ

# Clear environment and load required packages
rm(list = ls())
gc() # Clean unused objects in memory

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse//0.Enviroment.R")

dir.create("./3_Clustering_Tree/")

#===============================================================================
# Advanced Hierarchical K-Means Clustering
# Based on mRNA / ncRNA expression matrix + clustering tree results
#===============================================================================
# Load packages
library(tidyverse)
library(factoextra)  # Core package: hkmeans + visualization
library(cluster)
library(ape)
library(data.table)

#===============================================================================
# 1. Load existing data (fully compatible with original code)
#===============================================================================
group_info <- read.csv("../Outdata/1.rawdata/0.sample_info_all.csv", stringsAsFactors = FALSE)
# mRNA expression matrix
exp_matrix <- fread("../Outdata/5.all_data_harmony/combat_all_mRNA.csv")
# # Load ncRNA expression matrix (rows = genes, columns = individual samples)
# eRNA <- fread("../Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# miRNA <- fread("../Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# lncRNA <- fread("../Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# exp_matrix <- rbind(eRNA,miRNA,lncRNA)


rownames(exp_matrix) <- exp_matrix$V1
exp_matrix <- data.frame(exp_matrix[,-1], row.names = rownames(exp_matrix)) # Force retain row names

# Transpose matrix: rows = samples, columns = genes
exp_matrix_t <- t(exp_matrix)
exp_df <- as.data.frame(exp_matrix_t) %>% rownames_to_column("Sample")
exp_batch_df <- merge(exp_df, group_info[, c("Sample", "Batch", "Group")], by = "Sample")

# Calculate mean expression per Batch
batch_group_map <- exp_batch_df %>% select(Batch, Group) %>% distinct()
batch_exp_df <- exp_batch_df %>% 
  group_by(Batch) %>% 
  summarise(across(-c(Sample, Group), mean), .groups = "drop")

batch_exp_df <- merge(batch_exp_df, batch_group_map, by = "Batch")
batch_exp_df <- batch_exp_df %>% mutate(Batch_New = paste0(Batch, "_", Group))

# Construct final matrix: rows = Batch_New, columns = genes
batch_exp_matrix <- batch_exp_df %>% 
  column_to_rownames("Batch_New") %>% 
  select(-Batch, -Group) %>% 
  as.matrix()

# Remove sample GSE233478_Normal: it always clusters with other tumor samples
# Delete row named GSE233478_Normal
batch_exp_matrix <- batch_exp_matrix[ !rownames(batch_exp_matrix) %in% "GSE233478_Normal", ]

#===============================================================================
# 2. Data standardization (mandatory preprocessing for clustering)
#===============================================================================
df <- scale(batch_exp_matrix)  # Standardization: mean = 0, variance = 1

#===============================================================================
# 3. Run Hierarchical K-Means clustering (hkmeans)
#===============================================================================
set.seed(123)  # Fix random seed for reproducible results
k <- 6  # Adjustable cluster number: 2/3/4/5
res.hk <- hkmeans(df, k)  # Core function: Hierarchical K-Means

# Check clustering results
res.hk
table(res.hk$cluster)  # Sample count per cluster
# Cluster 1 2 3 4 5 6 (ncRNA dataset)
# Count  5 3 2 6 2 3 

# Custom color palette:
# 6 colors matching cluster 1-6
custom_palette <- c(
  "#D8C7C2",  # Cluster 1
  "#CECCC7",  # Cluster 2
  "#C8CCD1",  # Cluster 3
  "#E8E7D7",  # Cluster 4
  "#C4D2BE",  # Cluster 5
  "#9199A7"   # Cluster 6
)


#===============================================================================
# 4. Generate publication-style advanced clustering tree with colored boxes
#===============================================================================
pdf("./3_Clustering_Tree/2_hkmeans_clustering_tree2.pdf", width = 10, height = 7)
# pdf("./3_Clustering_Tree/2_hkmeans_clustering_tree_ncRNA2.pdf", width = 10, height = 7)
fviz_dend(
  res.hk, 
  cex = 0.5,                # Label font size
  palette = "jco",          # JCO journal standard color scheme
  # palette = custom_palette,
  rect = TRUE,              # Draw cluster bounding boxes
  rect_border = "jco",      # Box border color
  # rect_border = custom_palette,
  rect_fill = TRUE,         # Fill cluster boxes with color
  show_labels = TRUE,       # Display sample names
  main = paste0("Hierarchical K-Means Clustering (k = ", k, ")")
)
dev.off()

#===============================================================================
# 5. Generate PCA clustering scatter plot for clear group visualization
#===============================================================================
pdf("./3_Clustering_Tree//3_hkmeans_cluster_scatter.pdf", width = 6, height = 5)
# pdf("./3_Clustering_Tree//3_hkmeans_cluster_scatter_ncRNA.pdf", width = 6, height = 5)
fviz_cluster(
  res.hk, 
  palette = "jco",          # Color scheme
  ellipse.type = "convex",  # Convex hull boundary for clusters
  repel = TRUE,             # Avoid overlapping labels
  star.plot = FALSE,        
  show.clust.cent = TRUE,   # Show cluster centroids
  ggtheme = theme_classic(),# Clean plot theme
  main = "Cluster Scatter Plot (Hierarchical K-Means)"
)
dev.off()




# 
# 
# 
# 
# 
# ################################### 2. t-SNE Visualization ##########################################
# # Load required packages
# library(tidyverse)
# library(Rtsne)
# library(data.table)
# 
# # 1. Reuse preloaded raw data to guarantee consistent sample metadata
# # Load sample metadata (Sample, Batch, Group (Phenotype: Normal/AML))
# group_info <- read.csv("./Outdata/1.rawdata/0.sample_info_all.csv", stringsAsFactors = FALSE)
# 
# # # Load mRNA expression matrix (rows = genes, columns = individual samples)
# # exp_matrix <- read.csv("./Outdata/5.all_data_harmony/combat_all_mRNA.csv", row.names = 1, header = T)
# 
# # Load ncRNA expression matrix (rows = genes, columns = individual samples)
# eRNA <- fread("./Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# miRNA <- fread("./Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# lncRNA <- fread("./Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# exp_matrix <- rbind(eRNA,miRNA,lncRNA)
# rownames(exp_matrix) <- exp_matrix$V1
# exp_matrix <- data.frame(exp_matrix[,-1], row.names = rownames(exp_matrix)) # Force retain row names
# 
# # 2. Perform PCA dimensional reduction (t-SNE requires low-dimensional input to avoid curse of dimensionality)
# # Note: prcomp computes PCA by columns (samples). exp_matrix rows=genes, columns=samples, no extra transpose needed
# pca <- prcomp(t(exp_matrix), scale. = TRUE) 
# pca_data <- pca$x[, 1:50]  # Extract top 50 principal components (retain majority variance for t-SNE input)
# 
# # 3. Run t-SNE algorithm
# tsne_result <- Rtsne(pca_data, perplexity = 10, dims = 2, check_duplicates = FALSE)
# 
# # 4. Critical step: Construct t-SNE plotting dataframe matched with sample grouping metadata
# # (1) Extract sample names (1:1 correspondence with PCA/t-SNE output)
# sample_names <- rownames(pca_data)  # Sample names identical to exp_matrix column names
# 
# # (2) Build base t-SNE dataframe
# tsne_df <- data.frame(
#   Sample = sample_names,  # New Sample column for metadata matching
#   TSNE1 = tsne_result$Y[, 1],
#   TSNE2 = tsne_result$Y[, 2],
#   stringsAsFactors = FALSE
# )
# 
# # (3) Merge sample grouping information (Phenotype=Group: Normal/AML, Batch info available)
# tsne_df <- merge(tsne_df, group_info[, c("Sample", "Group", "Batch")], by = "Sample")
# # Rename Group column to Phenotype to match analysis requirement
# tsne_df <- tsne_df %>% rename(Phenotype = Group)
# 
# # (Optional) Concatenate Batch and Phenotype label (e.g. GSE106272_Tumor, consistent with previous workflow)
# tsne_df <- tsne_df %>% 
#   mutate(Batch_New = paste0(Batch, "_", Phenotype))
# 
# # Plotting
# # pdf("./Outplot/1.Sample_Info/3_mRNA_BatchNew_Group_Combine_t-SNE.pdf", width = 12, height = 8, pointsize = 12)
# pdf("./Outplot/1.Sample_Info/4_ncRNA_BatchNew_Group_Combine_t-SNE.pdf", width = 12, height = 8, pointsize = 12)
# ggplot(tsne_df, aes(
#   x = TSNE1, 
#   y = TSNE2, 
#   color = Batch_New,  # Concatenated Batch_New label as color mapping
#   shape = Phenotype   # Different shapes for distinct Phenotype groups
# )) +
#   geom_point(size = 2, alpha = 0.7) +
#   theme_minimal() +
#   # 1. Remove background grid lines
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     # 2. Add plot border: restore panel and axis frame
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # Empty black panel border
#     axis.line = element_line(color = "black", linewidth = 0.8)  # Axis frame (optional for sharper border)
#   ) +
#   labs(
#     # title = "t-SNE Visualization of mRNA Expression: Group & Phenotype",
#     title = "t-SNE Visualization of ncRNA Expression: Group & Phenotype",
#     x = "t-SNE Dimension 1",
#     y = "t-SNE Dimension 2",
#     color = "Batch + Phenotype",
#     shape = "Phenotype (Normal/AML)"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 14, face = "italic"),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 10),
#     legend.position = "right"
#   )
# dev.off()
# 
# # save.image("./Outplot/1.Sample_Info/mRNA-t-SNE.RData")
# # save.image("./Outplot/1.Sample_Info/ncRNA-t-SNE.RData")
