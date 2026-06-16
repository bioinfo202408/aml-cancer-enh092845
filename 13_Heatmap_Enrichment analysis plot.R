# Reference: https://mp.weixin.qq.com/s/dmeB9m5ePK74MXeMOBmpug
# Top Journal Figure Library | Issue 059: Gene Clustering, Enrichment Annotation Heatmap
# Line plot + Heatmap + Enrichment Analysis + Pathway + Single-cell
# Reference: https://mp.weixin.qq.com/s/FvTsDSOUgobz4hjouarN0g

# Clear environment and load required packages
rm(list = ls())
gc() # Clean unused objects in memory

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码//")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")
source("/home/weili/Project/AML/human/AML_combined_analyse//0.Enviroment.R")
dir.create("./13_Heatmap_Enrichment_Analysis_Figures/", showWarnings = F, recursive = T)



########################## Line Plot + Heatmap + Enrichment Analysis + Pathway + Single-cell ######################################
# 1. Load packages (installation guide at the end if not installed)
# devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)  # Human gene annotation database
library(ggsci)
library(ComplexHeatmap)
library(qs)




# 3. Load preprocessed data
load("../scRNA_analyse/20260407_Analysis_Results_GSE1116256/rawdata/Step6.celltype.markers.0.25.Rdata")
head(ct.marker)
seurat.data = qread(file = "../scRNA_analyse/20260407_Analysis_Results_GSE1116256//rawdata/Step3.annotation.qs")



# 4. Extract TOP 20 marker genes per cell cluster
markers <- ct.marker %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

# 5. Generate expression matrix from Seurat object (critical step!)
# Replace scRNA object name with your own Seurat object name!!!
st.data <- prepareDataFromscRNA(
  object = seurat.data,
  diffData = markers,
  showAverage = TRUE
)

# 6. GO enrichment analysis (add functional annotations for each cluster)
enrich <- enrichCluster(
  object = st.data,
  OrgDb = org.Hs.eg.db,
  type = "BP",          # Biological Process; can switch to "MF", "CC", "KEGG"
  organism = "hsa",
  pvalueCutoff = 0.5,
  topn = 5,            # Display top 5 enriched terms per cluster
  seed = 123
)

# 7. Randomly label 40 gene names on plot (adjustable number)
set.seed(123)
markGenes <- sample(unique(markers$gene), 40)


# Line plot only
pdf("./13_Heatmap_ClusterTree/1.pdf")
visCluster(object = st.data,           
           plotType  = "both")
dev.off()


pdf('./13_Heatmap_ClusterTree/SingleCell_Subpopulation_Enrichment_Analysis.pdf', height = 10, width = 16, onefile = F)
visCluster(
  object = st.data,
  plotType = "both",
  column_names_rot = 45,
  showRowNames = FALSE,                
  markGenes = markGenes,
  markGenesSide = "left",             
  annoTermData = enrich,              
  lineSide = "left",                  
  goCol = rep(ggsci::pal_d3()(length(unique(markers$cluster))), each = 5),
  goSize = "pval",                    
  addBar = TRUE,                      
  textbarPos = c(0.8, 0.2)           
)
dev.off()





### Customize subpopulation order
# Check current column order (exclude last two columns)
current_cols <- colnames(st.data$wide.res)[1:(ncol(st.data$wide.res)-2)]
print(current_cols)
# "HSC/MPP"     "Mono/Mac"    "B cell"      "T cell"      "Dendritic"   "Erythro"     "Plasma cell"

# Define target order (names must match original labels exactly)
desired_order <- c("Dendritic","Mono/Mac",  "T cell", "B cell", "Plasma cell", "HSC/MPP", "Erythro")

# Reorder columns of st.data$wide.res
st.data$wide.res <- st.data$wide.res[, c(desired_order, "gene", "cluster")]

# Adjust factor levels of cell_type in st.data$long.res to align line plot grouping order
st.data$long.res$cell_type <- factor(st.data$long.res$cell_type, levels = desired_order)

# Disable column clustering in visCluster to avoid automatic reordering
pdf('./13_Heatmap_ClusterTree/SingleCell_Subpopulation_Enrichment_Analysis_CustomOrder.pdf', height = 10, width = 14, onefile = F)
visCluster(
  object = st.data,
  plotType = "both",
  clusterColumns = FALSE,   # Critical: disable automatic column clustering
  column_names_rot = 45,
  showRowNames = FALSE,                
  markGenes = markGenes,
  markGenesSide = "left",             
  annoTermData = enrich,              
  lineSide = "left",                  
  goCol = rep(ggsci::pal_d3()(length(unique(markers$cluster))), each = 5),
  goSize = "pval",                    
  addBar = TRUE,                      
  textbarPos = c(0.8, 0.2)           
)
dev.off()









######################### FWSE Feature Selection Heatmap (Publication-level ComplexHeatmap) ######################################
# Use original file paths + automatic grouping + auto split up/down-regulated genes + high-resolution PDF
###########################################################################################

# Load packages
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)

# ===================== 1. Import raw expression data (unchanged from original code) 
# mRNA (active running module)
# miRNA_train_expr <- fread("../Outdata/2.train_test_data/1.2.expr_train_mRNA.csv")
# used_gene <- read.table("../Outdata/4.FWSE/4.1_FWSE_mRNA_8_flod.txt")

# miRNA
# miRNA_train_expr <- fread("../Outdata/2.train_test_data/1.2.expr_train_miRNA.csv")
# used_gene <- read.table("../Outdata/4.FWSE/4.1_FWSE_miRNA_8_flod.txt")

# lncRNA
# miRNA_train_expr <- fread("../Outdata/2.train_test_data/1.2.expr_train_lncRNA.csv")
# used_gene <- read.table("../Outdata/4.FWSE/4.1_FWSE_lncRNA_8_flod.txt")

# eRNA
miRNA_train_expr <- fread("../Outdata/2.train_test_data/1.2.expr_train_eRNA.csv")
used_gene <- read.table("../Outdata/4.FWSE/4.1_FWSE_eRNA_8_flod.txt")

{
sample_info_train <- read.csv("../Outdata/2.train_test_data/2.2.train_samples_info.csv", header = TRUE, row.names = 1)

# ===================== 2. Filter FWSE signature genes 
used_gene <- used_gene$V1
names_gene <- miRNA_train_expr$V1
miRNA_train_expr <- miRNA_train_expr[, -1]
miRNA_train_expr <- as.data.frame(miRNA_train_expr)
rownames(miRNA_train_expr) <- names_gene
head(rownames(miRNA_train_expr) )
hmExp <- miRNA_train_expr[used_gene, ]  # Rows: genes, Columns: samples

# Remove gene version suffix
rownames(hmExp) <- sapply(strsplit(rownames(hmExp), "\\."), `[`, 1)

# ===================== 3. Data preprocessing & cleaning
bad_rows <- apply(hmExp, 1, function(row) anyNA(row) | any(is.infinite(row)) | any(is.nan(row)))
hmExp <- hmExp[!bad_rows, , drop = FALSE]
cat("Gene count after cleaning: ", nrow(hmExp), "\n")

# ===================== ✅ Fix: Manual row-wise normalization for ComplexHeatmap 
hmExp <- t(scale(t(hmExp)))  # Gene-wise Z-score normalization (replaces scale="row")

# ===================== 4. Core step: Auto sort samples (Tumor first, Normal second)
sample_info_filter <- sample_info_train[colnames(hmExp), , drop = FALSE]

# Extract Tumor / Normal sample IDs
tumor_samp <- rownames(sample_info_filter)[sample_info_filter$Group == "Tumor"]
normal_samp <- rownames(sample_info_filter)[sample_info_filter$Group == "Normal"]

# Reorder columns: Tumor followed by Normal
hmExp_reordered <- hmExp[, c(tumor_samp, normal_samp)]

# Sample group annotation vector
group_anno <- c(rep("Tumor", length(tumor_samp)), rep("Normal", length(normal_samp)))

# ===================== 5. Plot color palette (matched target figure style) 
# Blue-White-Red gradient
heat_colors <- colorRamp2(c(-2, 0, 2), c("#3A3E96", "#F7F7F7", "#AD3D3E"))

# Group annotation colors
group_colors <- c(Tumor = "#AA899D", Normal = "#50A293")

# ===================== 6. Draw publication-standard ComplexHeatmap 
# Output path (reuse defined folder)
# pdf("./13_Heatmap_Enrichment_Analysis_Figures/mRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
# pdf("./13_Heatmap_Enrichment_Analysis_Figures/miRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
# pdf("./13_Heatmap_Enrichment_Analysis_Figures/lncRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
pdf("./13_Heatmap_Enrichment_Analysis_Figures/eRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)

# Top sample group annotation
column_ha <- HeatmapAnnotation(
  Type = group_anno,
  col = list(Type = group_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

# Generate heatmap (scale="row" removed)
ht <- Heatmap(
  hmExp_reordered,
  name = "z-score",
  cluster_columns = FALSE,        # Disable sample clustering
  cluster_rows = TRUE,            # Enable gene clustering
  show_row_names = FALSE,         # Hide gene labels (too many features)
  show_column_names = FALSE,      # Hide sample labels
  top_annotation = column_ha,     # Attach sample group annotation
  col = heat_colors,              # Expression color gradient
  row_km = 2,                     # Kmeans split into 2 clusters (up / down regulated); set to 1 if all genes share one trend
  column_split = factor(group_anno, levels = c("Tumor", "Normal")),  # Add split lines between groups
  border = TRUE,                  # Draw heatmap outer border
  heatmap_legend_param = list(
    title = "Expression\nz-score",
    title_position = "leftcenter-rot"
  )
)

# Render heatmap object
ht <- draw(ht)

# ===================== 7. Auto extract up-regulated / down-regulated gene clusters
row_clusters <- row_order(ht)

if (length(row_clusters) == 2) {
  cluster1_genes <- rownames(hmExp_reordered)[row_clusters[[1]]]
  cluster2_genes <- rownames(hmExp_reordered)[row_clusters[[2]]]
} else {
  stop("Detected cluster count not equal to 2, please check data")
}

dev.off()
}
# ===================== 8. Export gene lists of two clusters
write.table(cluster1_genes,
            file = "./13_Heatmap_Enrichment_Analysis_Figures/eRNA_up_FWSE_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(cluster2_genes,
            file = "./13_Heatmap_Enrichment_Analysis_Figures/eRNA_down_FWSE_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)




######################### FWSE Feature Heatmap - External Validation Cohort GSE103424 ######################################
# Load packages
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)

#========= 1. Load external validation cohort expression matrix GSE103424 
vali_exp <- fread("../../TCGA_data/all_GSE103424_GTEX_sample_clean_expr.csv")
gene_name <- vali_exp$gene_name
vali_exp <- vali_exp[,c(701:790)]
rownames(vali_exp) <- gene_name

#========= 2. Batch import four types of feature genes and loop plotting (miRNA/mRNA/lncRNA/eRNA)
# Mapping between gene file path and output figure name
gene_list = list(
  miRNA = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_miRNA_8_flod.txt", outname="miRNA"),
  mRNA  = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_mRNA_8_flod.txt", outname="mRNA"),
  lncRNA= list(gfile="../../Outdata/4.FWSE/4.1_FWSE_lncRNA_8_flod_rename.txt", outname="lncRNA"),
  eRNA  = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_eRNA_8_flod.txt", outname="eRNA")
)

# Batch loop for heatmap generation
for(nam in names(gene_list)){
  info = gene_list[[nam]]
  gene_set = read.table(info$gfile)
  hmExp = vali_exp[rownames(vali_exp) %in% gene_set$V1, ]
  
  # 3. Data cleaning: remove rows containing NA / infinite values
  bad_rows <- apply(hmExp, 1, function(row) anyNA(row) | any(is.infinite(row)) | any(is.nan(row)))
  hmExp <- hmExp[!bad_rows, , drop = FALSE]
  cat(info$outname,"Gene count after cleaning: ", nrow(hmExp), "\n")
  
  # Row-wise Z-score normalization for genes
  hmExp <- t(scale(t(hmExp)))
  
  #========= Critical: Auto grouping by sample ID prefix (SRR = Tumor / others = Normal)
  samp_all = colnames(hmExp)
  group_vec = ifelse(grepl("^SRR",samp_all),"Tumor","Normal")
  # Reorder samples: Tumor first, Normal second
  tumor_samp = samp_all[group_vec=="Tumor"]
  normal_samp = samp_all[group_vec=="Normal"]
  hmExp_reordered = hmExp[,c(tumor_samp,normal_samp)]
  group_anno = c(rep("Tumor",length(tumor_samp)),rep("Normal",length(normal_samp)))
  
  # Color scheme
  heat_colors <- colorRamp2(c(-2, 0, 2), c("#3A3E96", "#F7F7F7", "#AD3D3E"))
  group_colors <- c(Tumor = "#AA899D", Normal = "#50A293")
  
  # Sample top annotation
  column_ha <- HeatmapAnnotation(
    Type = group_anno,
    col = list(Type = group_colors),
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  # Open PDF device
  pdf(paste0(".//",info$outname,"_FWSE_ComplexHeatmap_GSE165656.pdf"), width = 5, height = 10)
  
  # Generate heatmap
  ht <- Heatmap(
    hmExp_reordered,
    name = "z-score",
    cluster_columns = FALSE,        
    cluster_rows = TRUE,            
    show_row_names = FALSE,         
    show_column_names = FALSE,      
    top_annotation = column_ha,     
    col = heat_colors,              
    row_km = 2,                     # Kmeans split into up/down regulated clusters
    column_split = factor(group_anno, levels = c("Tumor", "Normal")),
    border = TRUE,                  
    heatmap_legend_param = list(
      title = "Expression\nz-score",
      title_position = "leftcenter-rot"
    )
  )
  ht <- draw(ht)
  dev.off()
  
  # Extract genes from two Kmeans clusters
  row_clusters <- row_order(ht)
  if (length(row_clusters) == 2) {
    cluster1_genes <- rownames(hmExp_reordered)[row_clusters[[1]]]
    cluster2_genes <- rownames(hmExp_reordered)[row_clusters[[2]]]
  } else {
    warning(paste0(info$outname,"Cluster number not equal to 2, skip gene export"))
    next
  }
  
  # Output gene lists
  write.table(cluster1_genes,
              file = paste0("./",info$outname,"_up_FWSE_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  write.table(cluster2_genes,
              file = paste0(".//",info$outname,"_down_FWSE_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
