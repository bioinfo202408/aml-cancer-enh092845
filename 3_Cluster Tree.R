# 参考：https://mp.weixin.qq.com/s/def1EiaRf-eBrQrKOketkQ

# 清空环境并加载必要的包
rm(list = ls())
gc() #清理内存中不再使用的对象

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse//0.Enviroment.R")

dir.create("./3_聚类树/")

#===============================================================================
# 高级分层 K-Means 聚类（Hierarchical K-Means）
# 基于你的 mRNA / ncRNA 表达矩阵 + 聚类树结果
#===============================================================================
# 加载包
library(tidyverse)
library(factoextra)  # 核心：hkmeans + 可视化
library(cluster)
library(ape)
library(data.table)

#===============================================================================
# 1. 读取你已有的数据（完全沿用你的代码）
#===============================================================================
group_info <- read.csv("../Outdata/1.rawdata/0.sample_info_all.csv", stringsAsFactors = FALSE)
# mRNA
exp_matrix <- fread("../Outdata/5.all_data_harmony/combat_all_mRNA.csv")
# # 读取ncRNA表达矩阵（行=基因，列=单个样本）
# eRNA <- fread("../Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# miRNA <- fread("../Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# lncRNA <- fread("../Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# exp_matrix <- rbind(eRNA,miRNA,lncRNA)



rownames(exp_matrix) <- exp_matrix$V1
exp_matrix <- data.frame(exp_matrix[,-1], row.names = rownames(exp_matrix)) #rownames.force=TRUE 强制保留行名

# 转置：行=样本，列=基因
exp_matrix_t <- t(exp_matrix)
exp_df <- as.data.frame(exp_matrix_t) %>% rownames_to_column("Sample")
exp_batch_df <- merge(exp_df, group_info[, c("Sample", "Batch", "Group")], by = "Sample")

# 按 Batch 求平均
batch_group_map <- exp_batch_df %>% select(Batch, Group) %>% distinct()
batch_exp_df <- exp_batch_df %>% 
  group_by(Batch) %>% 
  summarise(across(-c(Sample, Group), mean), .groups = "drop")

batch_exp_df <- merge(batch_exp_df, batch_group_map, by = "Batch")
batch_exp_df <- batch_exp_df %>% mutate(Batch_New = paste0(Batch, "_", Group))

# 构建最终矩阵：行=Batch_New，列=基因
batch_exp_matrix <- batch_exp_df %>% 
  column_to_rownames("Batch_New") %>% 
  select(-Batch, -Group) %>% 
  as.matrix()

# 去除·GSE233478_Normal这个样本：因为他总是和其它tumor聚合在一起
# 删除行名为 GSE233478_Normal 的行
batch_exp_matrix <- batch_exp_matrix[ !rownames(batch_exp_matrix) %in% "GSE233478_Normal", ]

#===============================================================================
# 2. 数据标准化（必须！聚类前提）
#===============================================================================
df <- scale(batch_exp_matrix)  # 标准化：均值=0，方差=1

#===============================================================================
# 3. 执行 分层 K-Means 聚类（hkmeans）
#===============================================================================
set.seed(123)  # 固定随机数，结果可重复
k <- 6  # 你可以自由修改：2/3/4/5 类
res.hk <- hkmeans(df, k)  # 核心函数：分层 K-Means

# 查看结果
res.hk
table(res.hk$cluster)  # 每类样本数量
# 1 2 3 4 5 6   ncRNA
# 5 3 2 6 2 3 

# 自定义颜色：
# 你给的 6个颜色（对应簇1-6）
custom_palette <- c(
  "#D8C7C2",  # 簇1
  "#CECCC7",  # 簇2
  "#C8CCD1",  # 簇3
  "#E8E7D7",  # 簇4
  "#C4D2BE",   # 簇5
  "#9199A7"    #6
)


#===============================================================================
# 4. 绘制高级聚类树（带颜色方框，顶刊风格）
#===============================================================================
pdf("./3_聚类树/2_hkmeans_聚类树2.pdf", width = 10, height = 7)
# pdf("./3_聚类树/2_hkmeans_聚类树_ncRNA2.pdf", width = 10, height = 7)
fviz_dend(
  res.hk, 
  cex = 0.5,                # 字体大小
  palette = "jco",          # JCO 顶刊配色
  # palette = custom_palette,
  rect = TRUE,              # 画聚类方框
  rect_border = "jco",      # 方框颜色
  # rect_border = custom_palette,      # 方框颜色
  rect_fill = TRUE,         # 方框填充
  show_labels = TRUE,       # 显示样本名
  main = paste0("Hierarchical K-Means Clustering (k = ", k, ")")
)
dev.off()

#===============================================================================
# 5. 绘制聚类散点图（PCA降维，清晰展示分组）
#===============================================================================
pdf("./3_聚类树//3_hkmeans_聚类散点图.pdf", width = 6, height = 5)
# pdf("./3_聚类树//3_hkmeans_聚类散点图_ncRNA.pdf", width = 6, height = 5)
fviz_cluster(
  res.hk, 
  palette = "jco",          # 配色
  ellipse.type = "convex",  # 凸包包围
  repel = TRUE,             # 防止标签重叠
  star.plot = FALSE,        
  show.clust.cent = TRUE,   # 显示聚类中心
  ggtheme = theme_classic(),# 干净主题
  main = "Cluster Scatter Plot (Hierarchical K-Means)"
)
dev.off()




# 
# 
# 
# 
# 
# ################################### 2.T-SNE可视化##########################################
# # 加载必需包
# library(tidyverse)
# library(Rtsne)
# library(data.table)
# 
# # 1. 先加载并复用之前的原始数据（确保样本信息一致）
# # 读取样本元数据（含Sample、Batch、Group（Phenotype：正常/AML）对应关系）
# group_info <- read.csv("./Outdata/1.rawdata/0.sample_info_all.csv", stringsAsFactors = FALSE)
# 
# # # 读取mRNA表达矩阵（行=基因，列=单个样本）
# # exp_matrix <- read.csv("./Outdata/5.all_data_harmony/combat_all_mRNA.csv", row.names = 1, header = T)
# 
# # 读取ncRNA表达矩阵（行=基因，列=单个样本）
# eRNA <- fread("./Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# miRNA <- fread("./Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# lncRNA <- fread("./Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# exp_matrix <- rbind(eRNA,miRNA,lncRNA)
# rownames(exp_matrix) <- exp_matrix$V1
# exp_matrix <- data.frame(exp_matrix[,-1], row.names = rownames(exp_matrix)) #rownames.force=TRUE 强制保留行名
# 
# # 2. 对表达矩阵进行PCA降维（t-SNE的输入需低维数据，避免高维灾难）
# # 注意：prcomp默认按列（样本）进行PCA，exp_matrix行=基因、列=样本，无需额外转置
# pca <- prcomp(t(exp_matrix), scale. = TRUE) 
# pca_data <- pca$x[, 1:50]  # 取前50个主成分（保留大部分变异，作为t-SNE输入）
# 
# # 3. 运行t-SNE
# tsne_result <- Rtsne(pca_data, perplexity = 10, dims = 2, check_duplicates = FALSE)
# 
# # 4. 关键步骤：构建t-SNE可视化数据框（匹配样本分组信息）
# # （1）提取样本名（与PCA/t-SNE结果一一对应）
# sample_names <- rownames(pca_data)  # 样本名与exp_matrix列名一致
# 
# # （2）构建基础t-SNE数据框
# tsne_df <- data.frame(
#   Sample = sample_names,  # 新增Sample列，用于匹配分组信息
#   TSNE1 = tsne_result$Y[, 1],
#   TSNE2 = tsne_result$Y[, 2],
#   stringsAsFactors = FALSE
# )
# 
# # （3）匹配样本的分组信息（Phenotype=Group：正常/AML，同时可匹配Batch）
# tsne_df <- merge(tsne_df, group_info[, c("Sample", "Group", "Batch")], by = "Sample")
# # 重命名Group列为Phenotype（与你的需求一致）
# tsne_df <- tsne_df %>% rename(Phenotype = Group)
# 
# # （可选）给Batch拼接Phenotype信息（如GSE106272_Tumor，复用你之前的需求）
# tsne_df <- tsne_df %>% 
#   mutate(Batch_New = paste0(Batch, "_", Phenotype))
# 
# # 画图
# # pdf("./Outplot/1.样本信息/3_mRNA_BatchNew_Group_Combine_t-SNE.pdf", width = 12, height = 8, pointsize = 12)
# pdf("./Outplot/1.样本信息/4_ncRNA_BatchNew_Group_Combine_t-SNE.pdf", width = 12, height = 8, pointsize = 12)
# ggplot(tsne_df, aes(
#   x = TSNE1, 
#   y = TSNE2, 
#   color = Batch_New,  # 拼接后的Batch_New（如GSE106272_Tumor）作为颜色
#   shape = Phenotype   # 不同Group用不同形状
# )) +
#   geom_point(size = 2, alpha = 0.7) +
#   theme_minimal() +
#   # 1. 去除背景网格线（核心设置：关闭主网格和次网格）
#   theme(
#     panel.grid.major = element_blank(),  # 去除主要网格线
#     panel.grid.minor = element_blank(),  # 去除次要网格线
#     # 2. 添加方框（图形边框）：恢复面板边框和坐标轴边框
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 面板方框（无填充，黑色边框）
#     axis.line = element_line(color = "black", linewidth = 0.8)  # 坐标轴边框（可选，增强方框质感）
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
# # save.image("./Outplot/1.样本信息/mRNA-t-SNE.RData")
# # save.image("./Outplot/1.样本信息/ncRNA-t-SNE.RData")
