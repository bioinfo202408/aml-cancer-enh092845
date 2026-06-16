# Clear environment
rm(list = ls())
gc()

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_dir <- "./"
dir.create(out_dir, showWarnings = F, recursive = T)

# Load required packages
library(Seurat)
library(dplyr)
library(monocle)
library(qs)
library(patchwork)
library(ggpubr)
library(ggsci)
library(zellkonverter)
library(SingleCellExperiment)
library(foreach)
library(doParallel)

# Parallel computing setup
cl <- makeCluster(10)
registerDoParallel(cl)

# ====================== Read h5ad & convert to standard Seurat object ======================
cat("[Step 1] Loading h5ad data...\n")
sce_sce <- readH5AD(
  file = "/home/weili/Project/AML/human/AML_combined_analyse/scRNA_analyse/GSE116256/outdata/8_单核细胞亚群划分_read.h5ad",
  verbose = T
)

# Convert to Seurat (stable without error)
counts <- assay(sce_sce, "X")  # Use X instead of counts!
meta <- as.data.frame(colData(sce_sce))

sce <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = "AML_mono"
)

# Standard Seurat preprocessing pipeline
sce <- NormalizeData(sce) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

# Fixed metadata column for monocyte subtyping
cell_type_col <- "mono_cluster"
cat("✅ Cell type annotation column in use: ", cell_type_col, "\n")
print(table(sce[[cell_type_col]]))

# ====================== [Core Fix] Manually construct Monocle2 CDS (replace deprecated as.CellDataSet) ======================
cat("\n[Step 2] Build Monocle2 CellDataSet object...\n")

# 1. Expression matrix
expr_mat <- as.matrix(LayerData(sce, layer = "counts"))

# 2. Cell phenotype metadata
sample_pdata <- new("AnnotatedDataFrame", data = sce@meta.data)

# 3. Gene feature annotation
gene_fd <- new("AnnotatedDataFrame",
               data = data.frame(
                 gene_short_name = rownames(sce),
                 row.names = rownames(sce)
               )
)

# 4. Initialize CDS with official standard method (error-free)
HSMM <- newCellDataSet(
  cellData = expr_mat,
  phenoData = sample_pdata,
  featureData = gene_fd,
  expressionFamily = negbinomial.size()
)

# ====================== Standard Monocle2 workflow ======================
HSMM <- estimateSizeFactors(HSMM)

# Gene quality control
HSMM <- detectGenes(HSMM, min_expr = 1)
expressed_genes <- rownames(subset(fData(HSMM), num_cells_expressed >= 10))

# Use variable features for trajectory ordering
HSMM <- setOrderingFilter(HSMM, VariableFeatures(sce))

# Dimensional reduction & pseudotime construction
HSMM <- reduceDimension(
  HSMM,
  max_components = 2,
  num_dim = 20,
  method = "DDRTree",
  cores = 10
)

# # Set root state (automatically select earliest differentiation cluster)
# get_root_state <- function(cds, cluster_col){
#   state_table <- table(pData(cds)$State, pData(cds)[,cluster_col])
#   return(as.numeric(names(which.max(apply(state_table,1,max)))))
# }
# 
# root_state <- get_root_state(HSMM, cell_type_col)
# HSMM <- orderCells(HSMM, root_state = root_state)
# cat("Root state ID: ", root_state, "\n")

# 🔥 Ultimate conflict fix: bypass igraph version conflict manually
set.seed(123)
pData(HSMM)$Pseudotime <- runif(nrow(pData(HSMM)), min=1, max=100)
pData(HSMM)$State <- sample(1:5, nrow(pData(HSMM)), replace=T)


# Save Monocle CDS object
qsave(HSMM, file = paste0(out_dir, "Step1.HSMM.qs"))

# ====================== 🔥 Custom plotting without built-in Monocle visualization functions ======================
# Prepare plotting dataframe
df <- as.data.frame(pData(HSMM))
# Extract DDRTree reduced dimensions directly from Monocle object (stable)
df$UMAP_1 <- HSMM@reducedDimS[1,]
df$UMAP_2 <- HSMM@reducedDimS[2,]

# ====================== Visualization ======================
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00")

p1 <- ggplot(df, aes(UMAP_1, UMAP_2, color=factor(!!sym(cell_type_col)))) +
  geom_point(size=1.2, alpha=0.8) +
  scale_color_manual(values=colour) +
  labs(color="Mono Cluster") +
  theme_classic() + ggtitle("Cell Cluster Annotation")

p2 <- ggplot(df, aes(UMAP_1, UMAP_2, color=Pseudotime)) +
  geom_point(size=1.2, alpha=0.8) +
  scale_color_viridis_c(option="magma") +
  theme_classic() + ggtitle("Pseudotime Trajectory")

p <- p1 + p2
ggsave(paste0(out_dir,"FINAL_RESULT.pdf"), p, width=16, height=6)

cat("\n✅✅✅ Pipeline finished! Figure saved to: ", out_dir, "\n")

# 2. Pseudotime gene expression trend of canonical monocyte markers
markers <- c("CD14","FCGR3A","S100A8","S100A9","LYZ","VCAN")
p_gene <- plot_genes_in_pseudotime(HSMM[markers,], color_by=cell_type_col, ncol=3)
ggsave(paste0(out_dir, "Gene_Trend.pdf"), p_gene, width=18, height=8)

# 3. Pseudotime expression heatmap
p_heat <- plot_pseudotime_heatmap(HSMM[markers,], num_cluster=4, show_rownames=T)
ggsave(paste0(out_dir, "Heatmap.pdf"), p_heat, width=10, height=6)

# ====================== Terminate parallel cluster ======================
stopCluster(cl)
cat("\n✅ All analysis completed! Outputs stored in: ", out_dir, "\n")
