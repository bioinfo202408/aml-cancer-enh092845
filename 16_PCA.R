# 参考：https://figureya.online/FigureYa38PCA/FigureYa38PCA.html

# ==============================================================
# Step0 环境初始化 + 路径配置
# ==============================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# 工作目录
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("初始工作路径：", getwd(), "\n")
# 加载全局配置
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出文件夹
out_folder <- "./16_PCA/"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("当前输出目录：", getwd(), "\n")

# =============================================================
#=====================0.加载包、配色、椭圆函数=====================

options(stringsAsFactors = FALSE)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(sva)   # ComBat去批次必需

#固定配色
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C",
           "#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#置信椭圆函数
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

#=====================1.读取lncRNA表达+样本信息=====================
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_lncRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_miRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_mRNA.csv")
before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_eRNA.csv")



{
#提取基因名，构建表达矩阵【基因行，样本列】
gene_names <- before_pca_train$V1
before_pca_train <- before_pca_train[,-1]
rownames(before_pca_train) <- gene_names
expr_raw <- as.data.frame(before_pca_train)
rownames(expr_raw) <- gene_names

#样本注释
sample_info <- read.csv("../../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
sample_info$Group <- ifelse(sample_info$Group == 1, "Tumor", "Normal")
sample_info <- sample_info[sample_info$Set %in% "train", ]

#样本取交集
common_samp <- intersect(colnames(expr_raw), rownames(sample_info))
cat("共有样本：",length(common_samp),"\n")

#原始表达矩阵（未去批次）
expr_before <- expr_raw[, common_samp, drop = FALSE]
expr_before <- apply(expr_before,2,as.numeric)
rownames(expr_before) <- gene_names

#注释同步筛选
meta_df <- sample_info[common_samp,,drop=F]
batch_vec <- meta_df$Batch    #ComBat所需批次向量
group_vec <- meta_df$Group

#=====================2.ComBat去除批次效应【关键】=====================
#ComBat输入：基因行、样本列；batch=批次
expr_after <- ComBat(dat = expr_before, batch = batch_vec, mod = NULL)

#2.TPM筛选高变异lncRNA（剔除大量全0无信息基因，关键优化）
gene_sd <- apply(expr_after, 1, sd)
#优先筛选top1200高变异基因，可改1000/1500微调
top_var_genes <- names(sort(gene_sd, decreasing = TRUE)[1:1200])
expr_filter <- expr_after[top_var_genes, ]

#3.PCA参数适配TPM：lncRNA表达跨度悬殊，scale=F保留原始丰度差异
#修改封装函数内PCA：center=T, scale=F
#=====================PCA绘图函数=====================
plot_pca_wrap <- function(expr_in, save_prefix){
  # 关键：剔除整行方差=0的基因，杜绝奇异矩阵SVD报错
  g_sd <- apply(expr_in,1,sd)
  expr_in <- expr_in[g_sd > 1e-6, , drop=F]
  if(nrow(expr_in)<50) stop("有效基因过少")
  
  pca.res <- prcomp(t(expr_in), center = TRUE, scale. = FALSE)
  pca.pv <- summary(pca.res)$importance[2,]
  pc1_lab <- paste0("PC1 (",round(pca.pv[1]*100,2),"%)")
  pc2_lab <- paste0("PC2 (",round(pca.pv[2]*100,2),"%)")
  
  low_df <- as.data.frame(pca.res$x[,c("PC1","PC2")])
  low_df$Group <- factor(meta_df$Group)
  low_df$Batch <- meta_df$Batch
  
  #分组PCA
  p1 <- ggplot(low_df)+
    geom_point(aes(PC1,PC2,color=Group),size=2,alpha=0.6,shape=20)+
    scale_color_manual(values=mycol[1:nlevels(low_df$Group)])+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.position=c(0.01,0.99),legend.justification=c(0,1))+
    labs(x=pc1_lab,y=pc2_lab)
  p1_fin <- add_ellipase(p1, ellipase_pro=0.95, colour="dimgray", linetype=1, lwd=1)
  ggsave(paste0(save_prefix,"_PCA_Manual_Group.pdf"),p1_fin,width=5,height=5)
  
  #批次PCA
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


#=====================4.分别画【去批次前、去批次后】4张图=====================
#去批次前 Before
plot_pca_wrap(expr_before, save_prefix = "lncrna_Before_trian")
#ComBat去批次后 After
plot_pca_wrap(expr_filter, save_prefix = "lncrna_After_trian")


#去批次前 Before
plot_pca_wrap(expr_before, save_prefix = "mirna_Before_trian")
#ComBat去批次后 After
plot_pca_wrap(expr_filter, save_prefix = "mirna_After_trian")

#去批次前 Before
plot_pca_wrap(expr_before, save_prefix = "mrna_Before_trian")
#ComBat去批次后 After
plot_pca_wrap(expr_filter, save_prefix = "mrna_After_trian")

#去批次前 Before
plot_pca_wrap(expr_before, save_prefix = "erna_Before_trian")
#ComBat去批次后 After
plot_pca_wrap(expr_filter, save_prefix = "erna_After_trian")











# #下面只针对miRNA
# 
# 
# 
# # ==============================================================
# # Step0 环境初始化 + 路径配置
# # ==============================================================
# rm(list = ls())
# gc()
# options(stringsAsFactors = FALSE)
# 
# setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
# cat("初始工作路径：", getwd(), "\n")
# source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")
# 
# out_folder <- "./16_PCA/"
# dir.create(out_folder, showWarnings = F, recursive = T)
# setwd(out_folder)
# cat("当前输出目录：", getwd(), "\n")
# 
# #=====================0.加载包、配色、椭圆函数=====================
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
# #=====================PCA绘图函数（内置去零方差，避免SVD报错）=====================
# plot_pca_wrap <- function(expr_in, save_prefix){
#   g_sd <- apply(expr_in,1,sd)
#   expr_in <- expr_in[g_sd > 1e-6, , drop=F]
#   if(nrow(expr_in)<50) stop("有效基因不足，无法PCA")
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
#   #分组PCA
#   p1 <- ggplot(low_df)+
#     geom_point(aes(PC1,PC2,color=Group),size=2,alpha=0.6,shape=20)+
#     scale_color_manual(values=mycol[1:nlevels(low_df$Group)])+
#     theme_bw()+
#     theme(panel.grid = element_blank(),legend.position=c(0.01,0.99),legend.justification=c(0,1))+
#     labs(x=pc1_lab,y=pc2_lab)
#   p1_fin <- add_ellipase(p1, ellipase_pro=0.95, colour="dimgray", linetype=1, lwd=1)
#   ggsave(paste0(save_prefix,"_PCA_Manual_Group.pdf"),p1_fin,width=5,height=5)
#   
#   #批次PCA
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
# #=====================1.读取lncRNA数据=====================
# before_pca_train <- fread("../../Outdata/2.train_test_data/1.2.expr_train_miRNA.csv")
# 
# gene_names <- before_pca_train$V1
# before_pca_train <- before_pca_train[,-1]
# rownames(before_pca_train) <- gene_names
# expr_raw <- as.data.frame(before_pca_train)
# rownames(expr_raw) <- gene_names
# 
# #样本注释
# sample_info <- read.csv("../../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
# sample_info$Group <- ifelse(sample_info$Group == 1, "Tumor", "Normal")
# sample_info <- sample_info[sample_info$Set %in% "train", ]
# 
# common_samp <- intersect(colnames(expr_raw), rownames(sample_info))
# cat("共有样本：",length(common_samp),"\n")
# 
# expr_before <- expr_raw[, common_samp, drop = FALSE]
# expr_before <- apply(expr_before,2,as.numeric)
# rownames(expr_before) <- gene_names
# 
# meta_df <- sample_info[common_samp,,drop=F]
# batch_vec <- meta_df$Batch
# 
# #=====================2.ComBat校正=====================
# expr_after <- ComBat(dat = expr_before, batch = batch_vec, mod = NULL)
# 
# # 只在校正后现存基因里筛选高变异【修复下标越界核心】
# gene_sd <- apply(expr_after, 1, sd)
# select_num <- min(1200, length(gene_sd))
# top_var_genes <- names(sort(gene_sd, decreasing = TRUE)[1:select_num])
# expr_filter <- expr_after[top_var_genes, ]
# 
# #=====================3.出图=====================
# plot_pca_wrap(expr_before, save_prefix = "lncrna_Before_trian")
# plot_pca_wrap(expr_filter, save_prefix = "mirna_After_trian")
