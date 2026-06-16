# Reference: https://figureya.online/FigureYa38PCA/FigureYa38PCA.html

# ==============================================================
# Step0 Environment Initialization & Path Configuration
# ==============================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("Initial working directory: ", getwd(), "\n")
# Load global configuration script
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_folder <- "./16_PCA/"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("Current output directory: ", getwd(), "\n")

# =============================================================
#=====================0. Load packages, color palette & ellipse function=====================

options(stringsAsFactors = FALSE)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(sva)   # Required for ComBat batch correction

# Fixed color palette
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C",
           "#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

# Confidence ellipse drawing function
add_ellipase <- function(p, x="PC1", y="PC2", group="Group",
                         ellipase_pro = 0.95, linetype="solid", colour="black", lwd=1){
  obs <- p$data[,c(x,y,group)];colnames(obs)<-c("x","y","group")
  theta <- seq(-pi, pi, length.out = 50)
  circle <- cbind(cos(theta),sin(theta))
  ell <- ddply(obs,"group",function(dd){
    if(nrow(dd)<=2) return(NULL)
    mu <- colMeans(dd[,1:2])
    sig <- var(dd[,1:2])
    ed <- sqrt(qchisq(ellipase_pro,df=2))
    data.frame(sweep(circle%*%chol(sig)*ed,2,mu,"+"))
  })
  names(ell)[2:3] <- c("x","y")
  ell <- ddply(ell,"group",function(dd) dd[chull(dd$x,dd$y),])
  p + geom_polygon(data=ell,aes(x=x,y=y,group=group),fill=NA,colour=colour,lwd=lwd,linetype=linetype)
}

#=====================1. Import lncRNA expression matrix & sample metadata=====================
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_lncRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_miRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_mRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_eRNA.csv")



{
# Extract gene names and construct expression matrix (genes as rows, samples as columns)
gene_names <- before_pca_train$V1
before_pca_train <- before_pca_train[,-1]
rownames(before_pca_train) <- gene_names
expr_raw <- as.data.frame(before_pca_train)
rownames(expr_raw) <- gene_names

# Sample metadata
sample_info <- read.csv("../../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
sample_info$Group <- ifelse(sample_info$Group == 1, "Tumor", "Normal")
sample_info <- sample_info[sample_info$Set %in% "train", ]

# Intersect sample IDs
common_samp <- intersect(colnames(expr_raw), rownames(sample_info))
cat("Total overlapping samples: ",length(common_samp),"\n")

# Raw expression matrix (before batch correction)
expr_before <- expr_raw[, common_samp, drop = FALSE]
expr_before <- apply(expr_before,2,as.numeric)
rownames(expr_before) <- gene_names

# Filter metadata to matched samples
meta_df <- sample_info[common_samp,,drop=F]
batch_vec <- meta_df$Batch    # Batch vector required for ComBat
group_vec <- meta_df$Group

#=====================2. ComBat batch correction [Core Step]=====================
# ComBat input format: genes as rows, samples as columns; batch = batch information
expr_after <- ComBat(dat = expr_before, batch = batch_vec, mod = NULL)

# 2. Filter top variable genes from TPM matrix (remove thousands of zero-expression non-informative genes, key optimization)
gene_sd <- apply(expr_after, 1, sd)
# Select top 1200 most variable genes; adjust to 1000/1500 for fine-tuning
top_var_genes <- names(sort(gene_sd, decreasing = TRUE)[1:1200])
expr_filter <- expr_after[top_var_genes, ]

# 3. PCA parameter adaptation for TPM: lncRNA expression spans wide dynamic range, set scale=F to retain original abundance differences
# PCA setting inside wrapped function: center=T, scale=F
#=====================PCA plotting wrapper function=====================
plot_pca_wrap <- function(expr_in, save_prefix){
  # Critical step: remove genes with zero variance to avoid singular matrix SVD error
  g_sd <- apply(expr_in,1,sd)
  expr_in <- expr_in[g_sd > 1e-6, , drop=F]
  if(nrow(expr_in)<50) stop("Insufficient valid genes for PCA")
  
  pca.res <- prcomp(t(expr_in), center = TRUE, scale. = FALSE)
  pca.pv <- summary(pca.res)$importance[2,]
  pc1_lab <- paste0("PC1 (",round(pca.pv[1]*100,2),"%)")
  pc2_lab <- paste0("PC2 (",round(pca.pv[2]*100,2),"%)")
  
  low_df <- as.data.frame(pca.res$x[,c("PC1","PC2")])
  low_df$Group <- factor(meta_df$Group)
  low_df$Batch <- meta_df$Batch
  
  # PCA colored by disease group
  p1 <- ggplot(low_df)+
    geom_point(aes(PC1,PC2,color=Group),size=2,alpha=0.6,shape=20)+
    scale_color_manual(values=mycol[1:nlevels(low_df$Group)])+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.position=c(0.01,0.99),legend.justification=c(0,1))+
    labs(x=pc1_lab,y=pc2_lab)
  p1_fin <- add_ellipase(p1, ellipase_pro=0.95, colour="dimgray", linetype=1, lwd=1)
  ggsave(paste0(save_prefix,"_PCA_Manual_Group.pdf"),p1_fin,width=5,height=5)
  
  # PCA colored by batch
  n_batch <- length(unique(low_df$Batch))
  col_batch <- colorRampPalette(mycol)(n_batch)
  p2 <- ggplot(low_df)+
    geom_point(aes(PC1,PC2,color=Batch),size=2,alpha=0.6)+
    scale_color_manual(values=col_batch)+
    theme_bw()+theme(panel.grid=element_blank())+
    labs(x=pc1_lab,y=pc2_lab)
  ggsave(paste0(save_prefix,"_PCA_BatchCheck.pdf"),p2,width=8,height=5)
}
}


#=====================4. Generate 4 plots for data before and after batch correction=====================
# Before ComBat correction
plot_pca_wrap(expr_before, save_prefix = "lncrna_Before_trian")
# After ComBat correction
plot_pca_wrap(expr_filter, save_prefix = "lncrna_After_trian")

# Before ComBat correction
plot_pca_wrap(expr_before, save_prefix = "mirna_Before_trian")
# After ComBat correction
plot_pca_wrap(expr_filter, save_prefix = "mirna_After_trian")

# Before ComBat correction
plot_pca_wrap(expr_before, save_prefix = "mrna_Before_trian")
# After ComBat correction
plot_pca_wrap(expr_filter, save_prefix = "mrna_After_trian")

# Before ComBat correction
plot_pca_wrap(expr_before, save_prefix = "erna_Before_trian")
# After ComBat correction
plot_pca_wrap(expr_filter, save_prefix = "erna_After_trian")











# # The following section only processes miRNA data
# 
# 
# 
# # ==============================================================
# # Step0 Environment Initialization & Path Configuration
# # ==============================================================
# rm(list = ls())
# gc()
# options(stringsAsFactors = FALSE)
# 
# setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
# cat("Initial working directory: ", getwd(), "\n")
# source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")
# 
# out_folder <- "./16_PCA/"
# dir.create(out_folder, showWarnings = F, recursive = T)
# setwd(out_folder)
# cat("Current output directory: ", getwd(), "\n")
# 
# #=====================0. Load packages, color palette & ellipse function=====================
# library(ggplot2)
# library(plyr)
# library(dplyr)
# library(data.table)
# library(sva)
# 
# mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C",
#            "#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
# 
# add_ellipase <- function(p, x="PC1", y="PC2", group="Group",
#                          ellipase_pro = 0.95, linetype="solid", colour="black", lwd=1){
#   obs <- p$data[,c(x,y,group)];colnames(obs)<-c("x","y","group")
#   theta <- seq(-pi, pi, length.out = 50)
#   circle <- cbind(cos(theta),sin(theta))
#   ell <- ddply(obs,"group",function(dd){
#     if(nrow(dd)<=2) return(NULL)
#     mu <- colMeans(dd[,1:2])
#     sig <- var(dd[,1:2])
#     ed <- sqrt(qchisq(ellipase_pro,df=2))
#     data.frame(sweep(circle%*%chol(sig)*ed,2,mu,"+"))
#   })
#   names(ell)[2:3] <- c("x","y")
#   ell <- ddply(ell,"group",function(dd) dd[chull(dd$x,dd$y),])
#   p + geom_polygon(data=ell,aes(x=x,y=y,group=group),fill=NA,colour=colour,lwd=lwd,linetype=linetype)
# }
# 
# #=====================PCA plotting wrapper function (built-in zero-variance removal to avoid SVD error)=====================
# plot_pca_wrap <- function(expr_in, save_prefix){
#   g_sd <- apply(expr_in,1,sd)
#   expr_in <- expr_in[g_sd > 1e-6, , drop=F]
#   if(nrow(expr_in)<50) stop("Insufficient valid genes for PCA")
#   
#   pca.res <- prcomp(t(expr_in), center = TRUE, scale. = FALSE)
#   pca.pv <- summary(pca.res)$importance[2,]
#   pc1_lab <- paste0("PC1 (",round(pca.pv[1]*100,2),"%)")
#   pc2_lab <- paste0("PC2 (",round(pca.pv[2]*100,2),"%)")
#   
#   low_df <- as.data.frame(pca.res$x[,c("PC1","PC2")])
#   low_df$Group <- factor(meta_df$Group)
#   low_df$Batch <- meta_df$Batch
#   
#   # PCA colored by disease group
#   p1 <- ggplot(low_df)+
#     geom_point(aes(PC1,PC2,color=Group),size=2,alpha=0.6,shape=20)+
#     scale_color_manual(values=mycol[1:nlevels(low_df$Group)])+
#     theme_bw()+
#     theme(panel.grid = element_blank(),legend.position=c(0.01,0.99),legend.justification=c(0,1))+
#     labs(x=pc1_lab,y=pc2_lab)
#   p1_fin <- add_ellipase(p1, ellipase_pro=0.95, colour="dimgray", linetype=1, lwd=1)
#   ggsave(paste0(save_prefix,"_PCA_Manual_Group.pdf"),p1_fin,width=5,height=5)
#   
#   # PCA colored by batch
#   n_batch <- length(unique(low_df$Batch))
#   col_batch <- colorRampPalette(mycol)(n_batch)
#   p2 <- ggplot(low_df)+
#     geom_point(aes(PC1,PC2,color=Batch),size=2,alpha=0.6)+
#     scale_color_manual(values=col_batch)+
#     theme_bw()+theme(panel.grid=element_blank())+
#     labs(x=pc1_lab,y=pc2_lab)
#   ggsave(paste0(save_prefix,"_PCA_BatchCheck.pdf"),p2,width=8,height=5)
# }
# 
# #=====================1. Import lncRNA expression data=====================
# before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_miRNA.csv")
# 
# gene_names <- before_pca_train$V1
# before_pca_train <- before_pca_train[,-1]
# rownames(before_pca_train) <- gene_names
# expr_raw <- as.data.frame(before_pca_train)
# rownames(expr_raw) <- gene_names
# 
# # Sample metadata
# sample_info <- read.csv("../../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
# sample_info$Group <- ifelse(sample_info$Group == 1, "Tumor", "Normal")
# sample_info <- sample_info[sample_info$Set %in% "train", ]
# 
# common_samp <- intersect(colnames(expr_raw), rownames(sample_info))
# cat("Total overlapping samples: ",length(common_samp),"\n")
# 
# expr_before <- expr_raw[, common_samp, drop = FALSE]
# expr_before <- apply(expr_before,2,as.numeric)
# rownames(expr_before) <- gene_names
# 
# meta_df <- sample_info[common_samp,,drop=F]
# batch_vec <- meta_df$Batch
# 
# #=====================2. ComBat batch correction=====================
# expr_after <- ComBat(dat = expr_before, batch = batch_vec, mod = NULL)
# 
# # Only select variable genes from corrected matrix [core fix for subscript out of bounds]
# gene_sd <- apply(expr_after, 1, sd)
# select_num <- min(1200, length(gene_sd))
# top_var_genes <- names(sort(gene_sd, decreasing = TRUE)[1:select_num])
# expr_filter <- expr_after[top_var_genes, ]
# 
# #=====================3. Generate PCA figures=====================
# plot_pca_wrap(expr_before, save_prefix = "lncrna_Before_trian")
# plot_pca_wrap(expr_filter, save_prefix = "mirna_After_trian")
