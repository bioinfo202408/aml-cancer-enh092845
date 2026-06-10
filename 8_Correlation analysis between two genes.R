# 参考：https://mp.weixin.qq.com/s/MXKTxL0YYrj5gwvvG1OdEA

# 清空环境并加载必要的包
rm(list = ls())
gc() #清理内存中不再使用的对象

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

#查看临时目录
Sys.getenv("TMPDIR") 

dir.create("./8_两个基因间相关性分析/")

library(ggplot2)
library(ggpubr)
library(ggprism)
library(data.table)
library(psych)  # 用于corr.test，可同时输出相关系数和p值
library(dplyr)  # 用于数据整理
library(WGCNA)       # 用于corAndPvalue替代corr.test
library(foreach)     # 并行计算核心包
library(doParallel)  # 开启多核心
# 开启10个计算核心（根据服务器实际核心数调整，若不足10则自动用最大可用核心）
cl <- makeCluster(10)
registerDoParallel(cl)

### 去批次后数据
# mRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_mRNA.csv")
# miRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_miRNA.csv")
# eRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_eRNA.csv")
# lncRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
# #lncRNA ID 转换
# source("../0.gene_id_to_gene_name转换函数.R")
# lncRNA_expr <- convert_gene_id_to_name(
#   data = lncRNA_expr,  #只需要修改这里即可:lncRNA第一列是基因名
#   mapping_file_path = "/home/weili/Project/rawdata/TCGA/geneid_to_genename/lnc_genecode_gene_id_name_mapping1.txt"
# )
# colnames(lncRNA_expr)[1] <- "V1"
# #合并
# allRNA_expr <- rbind(mRNA_expr,miRNA_expr,eRNA_expr,lncRNA_expr)
# allRNA_expr <- as.data.frame(allRNA_expr)
# #去重复
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$V1), ]
# rownames(allRNA_expr) <- allRNA_expr$V1
# head(rownames(allRNA_expr))
# 
# write.csv(allRNA_expr,"../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
#### 去批次后数据
allRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
group_info <- read.csv("../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv")
tumor_sample <- rownames(group_info[group_info$Group == 1,])
tumor_expr <- allRNA_expr[,..tumor_sample]
rownames(tumor_expr) <- allRNA_expr$V1
head(rownames(tumor_expr))


"Enh092845" %in% rownames(tumor_expr)
"CLEC11A" %in% rownames(tumor_expr)

allRNA_expr <- tumor_expr



# #### TCGA只有AML数，有所有mRNA
# allRNA_expr <- fread("/home/weili/Project/multi_omic/Pathformer/data_TCGA/1.raw_data/TCGA-LAML/TCGA-LAML.mRNA.csv")
# # gene_id转换为gene_name
# # 安装并加载注释包（第一次用需要安装）
# if (!require("biomaRt")) {
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("biomaRt")
# }
# library(biomaRt)
# # 连接 Ensembl 数据库（人 hg38）
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # 你的 ID 列表（直接用你自己的 allRNA_expr$V1）
# ensembl_ids_with_version <- allRNA_expr$V1
# # 去掉版本号（.后面的数字）
# ensembl_ids <- gsub("\\..*", "", ensembl_ids_with_version)
# # 批量转换 ID → gene_name
# annot <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name"),
#   filters = "ensembl_gene_id",
#   values = ensembl_ids,
#   mart = mart
# )
# # 把基因名合并回你的数据框
# allRNA_expr$gene_name <- annot$external_gene_name[match(ensembl_ids, annot$ensembl_gene_id)]
# # 看结果
# head(allRNA_expr[, c("V1", "gene_name")])
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$gene_name), ]
# # 先去掉 gene_name 为 NA 的行（这是报错根源）
# allRNA_expr <- allRNA_expr[!is.na(allRNA_expr$gene_name), ]
# # 再去重（确保没有重复基因名）
# allRNA_expr <- allRNA_expr[!duplicated(allRNA_expr$gene_name), ]
# # 最后设置行名（现在绝对不会报错）
# rownames(allRNA_expr) <- allRNA_expr$gene_name
# # 查看结果
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
# # 新下载的包含所有基因的表达矩阵:GSE137851
# allRNA_expr <- fread("../../mapping/AML_100_GSE137851/3.combined_all_RNA_common_cols.csv")
# # 最后设置行名（现在绝对不会报错）
# rownames(allRNA_expr) <- allRNA_expr$gene_name
# # 查看结果
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




### 开始绘图
# ===================== 批量相关性绘图 =====================
library(ggplot2)
library(ggpubr)
library(ggside)
library(tidyverse)
library(ggprism)


# 基因列表
# genes <- c("CLEC11A", "LRRC4B", "SYT3", "C19orf81",
#            "SHANK1", "GPR32", "SMIM47", "ACP4")
genes <- c("CLEC11A")
query_gene <- "Enh092845"

# 强制所有表达量大于0
# ===================== 修复版：强制 allRNA_expr 全 > 0 =====================
# 只对数值列处理（跳过基因名字符列）
numeric_cols <- sapply(allRNA_expr, is.numeric)
if(sum(numeric_cols) > 0){
  # NA → 0
  allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) {
    x[is.na(x)] <- 0
    x
  }), .SDcols = numeric_cols]
  
  # 计算全局最小值（只算数值列）
  min_val <- min(as.matrix(allRNA_expr[, ..numeric_cols]), na.rm = TRUE)
  
  # 整体平移，让所有值 > 0
  if(min_val <= 0) {
    shift <- abs(min_val) + 0.001
    allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) x + shift), .SDcols = numeric_cols]
  }
  
  # 兜底：所有 <=0 的值强制变成 0.001
  allRNA_expr[, (names(allRNA_expr)[numeric_cols]) := lapply(.SD, function(x) {
    x[x <= 0] <- 0.001
    x
  }), .SDcols = numeric_cols]
}
# ==========================================================================

# 数据转换
expr_df <- as.data.frame(allRNA_expr)
expr_df <- log(expr_df+1)
rownames(expr_df) <- rownames(allRNA_expr)

# 基因存在检查
cat("======= 基因存在检查 =======\n")
cat(query_gene, ": ", query_gene %in% rownames(expr_df), "\n")
for(g in genes){ cat(g, ": ", g %in% rownames(expr_df), "\n") }

genes_exist <- genes[genes %in% rownames(expr_df)]
colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
            "#9467bd","#8c564b","#e377c2","#7f7f7f")

if(!dir.exists("./8_两个基因间相关性分析")) dir.create("./8_两个基因间相关性分析", recursive=T)
result_df <- data.frame()

# 批量循环（零警告版）
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
  
  pdf(paste0("./8_两个基因间相关性分析/",query_gene,"_",gene,".pdf"), width=6.5, height=5.5)
  
  p <- ggplot(data, aes(x=.data[[query_gene]], y=.data[[gene]])) +
    geom_point(color=col, size=3.2) +
    geom_smooth(method="lm", se=F, color="black", linewidth=0.8) +
    annotate("text", x=min(data[[query_gene]]), y=max(data[[gene]])*0.95,
             label=paste0("r = ", r, "\np = ", p_label),
             hjust=0, size=6) +
    labs(x=query_gene, y=gene) +
    
    # 上方分布密度图
    geom_xsidehistogram(
      aes(y = after_stat(density)),
      binwidth = 0.05, fill = "#1f77b4"
    ) +
    geom_xsidedensity(
      aes(y = after_stat(density)),
      linewidth = 1.2, color = "#CD6453"
    ) +
    scale_xsidey_continuous(labels = NULL)+
    
    # 右侧分布密度图
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

write.csv(result_df, "./8_两个基因间相关性分析/Enh092845_相关性汇总2.csv", row.names=F)
message("\n🎉 全部完成！")
