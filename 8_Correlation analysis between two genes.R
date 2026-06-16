# Reference: https://mp.weixin.qq.com/s/MXKTxL0YYrj5gwvvG1OdEA

# Clear environment and load required packages
rm(list = ls())
gc() # Clean unused objects in memory

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

# Check temporary directory
Sys.getenv("TMPDIR") 

dir.create("./8_Correlation_Analysis_Between_Two_Genes/")

library(ggplot2)
library(ggpubr)
library(ggprism)
library(data.table)
library(psych)  # For corr.test, outputs correlation coefficient and p-value simultaneously
library(dplyr)  # For data wrangling
library(WGCNA)       # corAndPvalue as alternative to corr.test
library(foreach)     # Core package for parallel computing
library(doParallel)  # Enable multi-threading
# Activate 10 computing cores (adjust based on actual server cores; auto uses max available if less than 10)
cl <- makeCluster(10)
registerDoParallel(cl)

### Batch-corrected expression data
# mRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_mRNA.csv")
# miRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# eRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# lncRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# # lncRNA ID conversion
# source("../0.gene_id_to_gene_name_conversion_function.R")
# lncRNA_expr <- convert_gene_id_to_name(
#   data = lncRNA_expr,  # Only modify this line: first column of lncRNA is gene name
#   mapping_file_path = "/home/weili/Project/rawdata/TCGA/geneid_to_genename/lnc_genecode_gene_id_name_mapping1.txt"
# )
# colnames(lncRNA_expr)[1] <- "V1"
# # Merge all RNA types
# allRNA_expr <- rbind(mRNA_expr,miRNA_expr,eRNA_expr,lncRNA_expr)
# allRNA_expr <- as.data.frame(allRNA_expr)
# # Remove duplicated genes
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$V1), ]
# rownames(allRNA_expr) <- allRNA_expr$V1
# head(rownames(allRNA_expr))
# 
# write.csv(allRNA_expr,"../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
#### Load batch-corrected combined RNA expression matrix
allRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
group_info <- read.csv("../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
tumor_sample <- rownames(group_info[group_info$Group == 1,])
tumor_expr <- allRNA_expr[,..tumor_sample]
rownames(tumor_expr) <- allRNA_expr$V1
head(rownames(tumor_expr))


"Enh092845" %in% rownames(tumor_expr)
"CLEC11A" %in% rownames(tumor_expr)

allRNA_expr <- tumor_expr



# #### TCGA-LAML mRNA dataset containing all mRNA genes
# allRNA_expr <- fread("/home/weili/Project/multi_omic/Pathformer/data_TCGA/1.raw_data/TCGA-LAML/TCGA-LAML.mRNA.csv")
# # Convert gene ID to gene symbol
# # Install and load annotation package (install on first use)
# if (!require("biomaRt")) {
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("biomaRt")
# }
# library(biomaRt)
# # Connect to Ensembl database (human hg38)
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # Input ID list (directly use allRNA_expr$V1)
# ensembl_ids_with_version <- allRNA_expr$V1
# # Remove version suffix (numbers after dot)
# ensembl_ids <- gsub("\\..*", "", ensembl_ids_with_version)
# # Batch convert Ensembl ID → gene symbol
# annot <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name"),
#   filters = "ensembl_gene_id",
#   values = ensembl_ids,
#   mart = mart
# )
# # Merge gene symbols back to expression dataframe
# allRNA_expr$gene_name <- annot$external_gene_name[match(ensembl_ids, annot$ensembl_gene_id)]
# # Check output
# head(allRNA_expr[, c("V1", "gene_name")])
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$gene_name), ]
# # Remove rows with NA gene symbols (root cause of errors)
# allRNA_expr <- allRNA_expr[!is.na(allRNA_expr$gene_name), ]
# # Deduplicate again to eliminate duplicate gene symbols
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$gene_name), ]
# # Assign row names (no error guaranteed)
# rownames(allRNA_expr) <- allRNA_expr$gene_name
# # Preview row names
# head(rownames(allRNA_expr))
# 
# "CLEC11A" %in% rownames(allRNA_expr)
# "LRRC4B" %in% rownames(allRNA_expr)
# "SYT3" %in%  rownames(allRNA_expr)
# "C19orf81" %in% rownames(allRNA_expr)
# "SHANK1" %in% rownames(allRNA_expr)
# "GPR32" %in%  rownames(allRNA_expr)
# "SMIM47" %in% rownames(allRNA_expr)
# "ACP4" %in% rownames(allRNA_expr)






# 
# # Combined expression matrix from GSE137851 containing all RNA types
# allRNA_expr <- fread("../../mapping/AML_100_GSE137851/3.combined_all_RNA_common_cols.csv")
# # Assign row names (no error guaranteed)
# rownames(allRNA_expr) <- allRNA_expr$gene_name
# # Preview row names
# head(rownames(allRNA_expr))
# 
# "Enh092845" %in% rownames(allRNA_expr)
# "CLEC11A" %in% rownames(allRNA_expr)
# "LRRC4B" %in% rownames(allRNA_expr)
# "SYT3" %in%  rownames(allRNA_expr)
# "C19orf81" %in% rownames(allRNA_expr)
# "SHANK1" %in% rownames(allRNA_expr)
# "GPR32" %in%  rownames(allRNA_expr)
# "SMIM47" %in% rownames(allRNA_expr)
# "ACP4" %in% rownames(allRNA_expr)




### Start plotting
# ===================== Batch Correlation Scatter Plot Pipeline =====================
library(ggplot2)
library(ggpubr)
library(ggside)
library(tidyverse)
library(ggprism)


# Gene list
# genes <- c("CLEC11A", "LRRC4B", "SYT3", "C19orf81",
#            "SHANK1", "GPR32", "SMIM47", "ACP4")
genes <- c("CLEC11A")
query_gene <- "Enh092845"

# Force all expression values greater than 0
# ===================== Fixed script: Ensure all values in allRNA_expr > 0 =====================
# Only process numeric columns (skip character gene ID columns)
numeric_cols <- sapply(allRNA_expr, is.numeric)
if(sum(numeric_cols) > 0){
  # Replace NA with 0
  allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) {
    x[is.na(x)] <- 0
    x
  }), .SDcols = numeric_cols]
  
  # Calculate global minimum value across numeric columns
  min_val <- min(as.matrix(allRNA_expr[, ..numeric_cols]), na.rm = TRUE)
  
  # Global shift to make all values positive
  if(min_val <= 0) {
    shift <- abs(min_val) + 0.001
    allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) x + shift), .SDcols = numeric_cols]
  }
  
  # Fallback: force all values ≤ 0 to 0.001
  allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) {
    x[x <= 0] <- 0.001
    x
  }), .SDcols = numeric_cols]
}
# ==========================================================================

# Data transformation
expr_df <- as.data.frame(allRNA_expr)
expr_df <- log(expr_df+1)
rownames(expr_df) <- rownames(allRNA_expr)

# Gene existence validation
cat("======= Gene Existence Check =======\n")
cat(query_gene, ": ", query_gene %in% rownames(expr_df), "\n")
for(g in genes){ cat(g, ": ", g %in% rownames(expr_df), "\n") }

genes_exist <- genes[genes %in% rownames(expr_df)]
colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
            "#9467bd","#8c564b","#e377c2","#7f7f7f")

if(!dir.exists("./8_Correlation_Analysis_Between_Two_Genes")) dir.create("./8_Correlation_Analysis_Between_Two_Genes", recursive=T)
result_df <- data.frame()

# Batch loop with zero warning handling
for(i in seq_along(genes_exist)){
  gene <- genes_exist[i]
  col  <- colors[i]
  
  data <- t(expr_df[c(query_gene, gene), ])
  data <- as.data.frame(data)
  data <- data[!rownames(data) %in% c("V1","gene_name"), ]
  
  data[[query_gene]] <- as.numeric(data[[query_gene]])
  data[[gene]] <- as.numeric(data[[gene]])
  data <- na.omit(data)
  
  tryCatch({
    cort <- cor.test(data[[query_gene]], data[[gene]])
    r <- round(cort$estimate, 3)
    p <- cort$p.value
    p_label <- if(is.numeric(p)) sprintf("%.2e", p) else "NA"
  }, error = function(e){
    r <<- NA
    p <<- NA
    p_label <<- "NA"
  })
  
  result_df <- rbind(result_df, data.frame(gene2=gene, pearson_r=r, p_value=p))
  
  pdf(paste0("./8_Correlation_Analysis_Between_Two_Genes/",query_gene,"_",gene,".pdf"), width=6.5, height=5.5)
  
  p <- ggplot(data, aes(x=.data[[query_gene]], y=.data[[gene]])) +
    geom_point(color=col, size=3.2) +
    geom_smooth(method="lm", se=F, color="black", linewidth=0.8) +
    annotate("text", x=min(data[[query_gene]]), y=max(data[[gene]])*0.95,
             label=paste0("r = ", r, "\np = ", p_label),
             hjust=0, size=6) +
    labs(x=query_gene, y=gene) +
    
    # Top marginal density histogram
    geom_xsidehistogram(
      aes(y = after_stat(density)),
      binwidth = 0.05, fill = "#1f77b4"
    ) +
    geom_xsidedensity(
      aes(y = after_stat(density)),
      linewidth = 1.2, color = "#CD6453"
    ) +
    scale_xsidey_continuous(labels = NULL)+
    
    # Right marginal density histogram
    geom_ysidehistogram(
      aes(x = after_stat(density)),
      binwidth = 0.05, fill = "#1f77b4"
    ) +
    geom_ysidedensity(
      aes(x = after_stat(density)),
      linewidth = 1.2, color = "#CD6453"
    ) +
    scale_ysidex_continuous(labels = NULL)+
    
    theme_bw() +
    theme(
      panel.border = element_rect(linewidth=0.5, fill=NA),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14),
      panel.grid = element_blank(),
      ggside.panel.scale = 0.2,
      ggside.axis.line = element_blank(),
      ggside.axis.ticks = element_blank()
    )
  
  print(p)
  dev.off()
  
  message(gene, "  → r = ", r, "  p = ", p_label)
}

write.csv(result_df, "./8_Correlation_Analysis_Between_Two_Genes/Enh092845_correlation_summary2.csv", row.names=F)
message("\n🎉 All analysis finished!")
