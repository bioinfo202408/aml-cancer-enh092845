# 参考：https://mp.weixin.qq.com/s/dmeB9m5ePK74MXeMOBmpug  顶刊图社 | 第059期：基因聚类富集注释热图  折线+热图+富集分析+通路+单细胞
# 参考：https://mp.weixin.qq.com/s/FvTsDSOUgobz4hjouarN0g

# 清空环境并加载必要的包
rm(list = ls())
gc() #清理内存中不再使用的对象

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码//")
cat("【初始化】当前工作路径：", getwd(), "\n\n")
source("/home/weili/Project/AML/human/AML_combined_analyse//0.Enviroment.R")
dir.create("./13_热图_富集分析_图/", showWarnings = F, recursive = T)



##########################折线+热图+富集分析+通路+单细胞######################################
# 1. 加载包（如果没装，看文末安装方法）
# devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)  # 人源数据库
library(ggsci)
library(ComplexHeatmap)
library(qs)




# 3. 加载你已经有的数据（你已经运行过了）
load("../scRNA_analyse/20260407分析结果_GSE1116256/rawdata/Step6.celltype.markers.0.25.Rdata")
head(ct.marker)
seurat.data = qread(file = "../scRNA_analyse/20260407分析结果_GSE1116256//rawdata/Step3.annotation.qs")



# 4. 提取每个细胞群 TOP 20 个 marker 基因
markers <- ct.marker %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

# 5. 从你的Seurat对象生成表达矩阵（关键！）
# 把 scRNA 对象替换成你自己的 Seurat 对象名称！！！
st.data <- prepareDataFromscRNA(
  object = seurat.data,
  diffData = markers,
  showAverage = TRUE
)

# 6. GO富集分析（给每个cluster加功能注释）
enrich <- enrichCluster(
  object = st.data,
  OrgDb = org.Hs.eg.db,
  type = "BP",          # 生物过程BP，也可换"MF"、"CC"、"KEGG"
  organism = "hsa",
  pvalueCutoff = 0.5,
  topn = 5,            # 每个cluster展示5个富集条目
  seed = 123
)

# 7. 随机标记40个基因名在图上（可改数量）
set.seed(123)
markGenes <- sample(unique(markers$gene), 40)


# line plot
pdf("./13_热图_聚类树/1.pdf")
visCluster(object = st.data,           
           plotType  = "both")
dev.off()


pdf('./13_热图_聚类树/单细胞不同亚群富集分析.pdf', height = 10, width = 16, onefile = F)
visCluster(
  object = st.data,
  plotType = "both",
  column_names_rot = 45,
  showRowNames = FALSE,                # 注意驼峰
  markGenes = markGenes,
  markGenesSide = "left",             # 改为 markGenesSide
  annoTermData = enrich,              # 改为 annoTermData
  lineSide = "left",                  # 改为 lineSide
  goCol = rep(ggsci::pal_d3()(length(unique(markers$cluster))), each = 5),  # goCol
  goSize = "pval",                    # goSize
  addBar = TRUE,                      # addBar
  textbarPos = c(0.8, 0.2)           # textbarPos
)
dev.off()





###  自定义不同亚群的顺序
# 查看当前列顺序（去掉最后两列）
current_cols <- colnames(st.data$wide.res)[1:(ncol(st.data$wide.res)-2)]
print(current_cols)
# "HSC/MPP"     "Mono/Mac"    "B cell"      "T cell"      "Dendritic"   "Erythro"     "Plasma cell"

# 定义你想要的顺序（注意名称必须与当前完全一致）
desired_order <- c("Dendritic","Mono/Mac",  "T cell", "B cell", "Plasma cell", "HSC/MPP", "Erythro")

# 重新排列 st.data$wide.res 的列
st.data$wide.res <- st.data$wide.res[, c(desired_order, "gene", "cluster")]

# 同时调整 st.data$long.res 中 cell_type 的因子水平，保持折线图分组顺序一致
st.data$long.res$cell_type <- factor(st.data$long.res$cell_type, levels = desired_order)

# 在 visCluster 中关闭列聚类，避免自动重排
pdf('./13_热图_聚类树/单细胞不同亚群富集分析_自定义顺序.pdf', height = 10, width = 14, onefile = F)
visCluster(
  object = st.data,
  plotType = "both",
  clusterColumns = FALSE,   # 关键：禁止对列重新聚类
  column_names_rot = 45,
  showRowNames = FALSE,                # 注意驼峰
  markGenes = markGenes,
  markGenesSide = "left",             # 改为 markGenesSide
  annoTermData = enrich,              # 改为 annoTermData
  lineSide = "left",                  # 改为 lineSide
  goCol = rep(ggsci::pal_d3()(length(unique(markers$cluster))), each = 5),  # goCol
  goSize = "pval",                    # goSize
  addBar = TRUE,                      # addBar
  textbarPos = c(0.8, 0.2)           # textbarPos
)
dev.off()









######################### FWSE特征选择热图（ComplexHeatmap 论文版） ######################################
# 完全使用你的数据路径 + 自动分组 + 自动拆分上下调基因 + 高清PDF
###########################################################################################

# 加载包
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# ===================== 1. 读取数据（和你原来完全一样） 
# mRNA（当前运行）
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

# ===================== 2. 筛选FWSE特征基因 
used_gene <- used_gene$V1
names_gene <- miRNA_train_expr$V1
miRNA_train_expr <- miRNA_train_expr[, -1]
miRNA_train_expr <- as.data.frame(miRNA_train_expr)
rownames(miRNA_train_expr) <- names_gene
head(rownames(miRNA_train_expr) )
hmExp <- miRNA_train_expr[used_gene, ]  # 行：基因，列：样本

# 去除基因版本号
rownames(hmExp) <- sapply(strsplit(rownames(hmExp), "\\."), `[`, 1)

# ===================== 3. 数据清洗
bad_rows <- apply(hmExp, 1, function(row) anyNA(row) | any(is.infinite(row)) | any(is.nan(row)))
hmExp <- hmExp[!bad_rows, , drop = FALSE]
cat("清洗后基因数量：", nrow(hmExp), "\n")

# ===================== ✅ 修复点：ComplexHeatmap 手动行标准化 
hmExp <- t(scale(t(hmExp)))  # 基因行 Z-score 标准化（替代 scale="row"）

# ===================== 4. 核心：自动拆分样本顺序（Tumor 在前，Normal 在后）
sample_info_filter <- sample_info_train[colnames(hmExp), , drop = FALSE]

# 提取 Tumor / Normal 样本
tumor_samp <- rownames(sample_info_filter)[sample_info_filter$Group == "Tumor"]
normal_samp <- rownames(sample_info_filter)[sample_info_filter$Group == "Normal"]

# 重新排序列：Tumor → Normal
hmExp_reordered <- hmExp[, c(tumor_samp, normal_samp)]

# 样本注释
group_anno <- c(rep("Tumor", length(tumor_samp)), rep("Normal", length(normal_samp)))

# ===================== 5. 绘图配色（和你目标图完全一致） 
# 蓝-白-红 配色
heat_colors <- colorRamp2(c(-2, 0, 2), c("#3A3E96", "#F7F7F7", "#AD3D3E"))

# 分组颜色
group_colors <- c(Tumor = "#AA899D", Normal = "#50A293")

# ===================== 6. 绘制 ComplexHeatmap（论文级） 
# 输出路径（直接用你的文件夹）
# pdf("./13_热图_富集分析_图/mRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
# pdf("./13_热图_富集分析_图/miRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
# pdf("./13_热图_富集分析_图/lncRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)
pdf("./13_热图_富集分析_图/eRNA_FWSE_ComplexHeatmap3.pdf", width = 5, height = 10)

# 顶部样本注释
column_ha <- HeatmapAnnotation(
  Type = group_anno,
  col = list(Type = group_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

# 绘制热图（已删除 scale="row"）
ht <- Heatmap(
  hmExp_reordered,
  name = "z-score",
  cluster_columns = FALSE,        # 样本不聚类
  cluster_rows = TRUE,            # 基因聚类
  show_row_names = FALSE,         # 不显示基因名（太多）
  show_column_names = FALSE,      # 不显示样本名
  top_annotation = column_ha,     # 样本分组注释
  col = heat_colors,              # 配色
  row_km = 2,                     # Kmeans 自动分成 2 簇（上调/下调）；如果都是上调或者下调，修改为1
  column_split = factor(group_anno, levels = c("Tumor", "Normal")),  # 分组分割线
  border = TRUE,                  # 热图边框
  heatmap_legend_param = list(
    title = "Expression\nz-score",
    title_position = "leftcenter-rot"
  )
)

# 输出热图
ht <- draw(ht)

# ===================== 7. 自动提取 上调/下调 基因簇
row_clusters <- row_order(ht)

if (length(row_clusters) == 2) {
  cluster1_genes <- rownames(hmExp_reordered)[row_clusters[[1]]]
  cluster2_genes <- rownames(hmExp_reordered)[row_clusters[[2]]]
} else {
  stop("聚类簇数量不是2，请检查")
}

dev.off()
}
# ===================== 8. 导出上下调基因
write.table(cluster1_genes,
            file = "./13_热图_富集分析_图/eRNA_up_FWSE_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(cluster2_genes,
            file = "./13_热图_富集分析_图/eRNA_down_FWSE_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)




######################### FWSE特征选择热图-外部验证集 GSE103424 ######################################
# 加载包
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)

#========= 1.读入外部验证集表达矩阵 GSE103424 
vali_exp <- fread("../../TCGA_data/all_GSE103424_GTEX_sample_clean_expr.csv")
gene_name <- vali_exp$gene_name
vali_exp <- vali_exp[,c(701:790)]
rownames(vali_exp) <- gene_name

#========= 2.批量读取四类特征基因，循环绘图（miRNA/mRNA/lncRNA/eRNA）
# 配置基因文件与输出文件名映射
gene_list = list(
  miRNA = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_miRNA_8_flod.txt", outname="miRNA"),
  mRNA  = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_mRNA_8_flod.txt", outname="mRNA"),
  lncRNA= list(gfile="../../Outdata/4.FWSE/4.1_FWSE_lncRNA_8_flod_rename.txt", outname="lncRNA"),
  eRNA  = list(gfile="../../Outdata/4.FWSE/4.1_FWSE_eRNA_8_flod.txt", outname="eRNA")
)

# 批量循环绘图
for(nam in names(gene_list)){
  info = gene_list[[nam]]
  gene_set = read.table(info$gfile)
  hmExp = vali_exp[rownames(vali_exp) %in% gene_set$V1, ]
  
  # 3.数据清洗：剔除NA/无穷值行
  bad_rows <- apply(hmExp, 1, function(row) anyNA(row) | any(is.infinite(row)) | any(is.nan(row)))
  hmExp <- hmExp[!bad_rows, , drop = FALSE]
  cat(info$outname,"清洗后基因数量：", nrow(hmExp), "\n")
  
  # 行Z-score标准化（基因标准化）
  hmExp <- t(scale(t(hmExp)))
  
  #=========关键：根据样本名自动分组 SRR=Tumor / 其余=Normal
  samp_all = colnames(hmExp)
  group_vec = ifelse(grepl("^SRR",samp_all),"Tumor","Normal")
  # 拆分样本顺序：Tumor在前，Normal在后
  tumor_samp = samp_all[group_vec=="Tumor"]
  normal_samp = samp_all[group_vec=="Normal"]
  hmExp_reordered = hmExp[,c(tumor_samp,normal_samp)]
  group_anno = c(rep("Tumor",length(tumor_samp)),rep("Normal",length(normal_samp)))
  
  # 配色
  heat_colors <- colorRamp2(c(-2, 0, 2), c("#3A3E96", "#F7F7F7", "#AD3D3E"))
  group_colors <- c(Tumor = "#AA899D", Normal = "#50A293")
  
  # 样本顶部注释
  column_ha <- HeatmapAnnotation(
    Type = group_anno,
    col = list(Type = group_colors),
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  # 打开pdf
  pdf(paste0(".//",info$outname,"_FWSE_ComplexHeatmap_GSE165656.pdf"), width = 5, height = 10)
  
  # 绘制热图
  ht <- Heatmap(
    hmExp_reordered,
    name = "z-score",
    cluster_columns = FALSE,        
    cluster_rows = TRUE,            
    show_row_names = FALSE,         
    show_column_names = FALSE,      
    top_annotation = column_ha,     
    col = heat_colors,              
    row_km = 2,                     # Kmeans分成上下调2簇
    column_split = factor(group_anno, levels = c("Tumor", "Normal")),
    border = TRUE,                  
    heatmap_legend_param = list(
      title = "Expression\nz-score",
      title_position = "leftcenter-rot"
    )
  )
  ht <- draw(ht)
  dev.off()
  
  # 提取聚类分簇基因
  row_clusters <- row_order(ht)
  if (length(row_clusters) == 2) {
    cluster1_genes <- rownames(hmExp_reordered)[row_clusters[[1]]]
    cluster2_genes <- rownames(hmExp_reordered)[row_clusters[[2]]]
  } else {
    warning(paste0(info$outname,"聚类不是2簇，跳过导出"))
    next
  }
  
  # 导出上下调基因
  write.table(cluster1_genes,
              file = paste0("./",info$outname,"_up_FWSE_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  write.table(cluster2_genes,
              file = paste0(".//",info$outname,"_down_FWSE_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


