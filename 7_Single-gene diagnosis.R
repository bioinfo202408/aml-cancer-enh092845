# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./7_单基因诊断/"
dir.create(out_dir, showWarnings = F, recursive = T)


# 读取表达量数据
Expr <- fread("../Outdata/5.all_data_harmony/combat_all_RNA_group_info.csv")
Expr <- Expr[,-1]

# 检测是否存在
"IRF1" %in% Expr$V1
"Enh092845" %in% Expr$V1
"CLEC11A" %in% Expr$V1

# 单基因诊断
library(pROC)
library(ggplot2)
library(dplyr)
library(data.table)

# 1. 把 Expr 转为 data.table（你现在的格式）
setDT(Expr)

# 2. 提取三个基因的表达 + 构建分组标签
expr_sub <- Expr[V1 %in% c("IRF1", "Enh092845", "CLEC11A")] %>%
  melt(id.vars = "V1", variable.name = "Sample", value.name = "Expression") %>%
  # 从样本名提取分组：Tumor/Normal
  mutate(Group = ifelse(grepl("_Tumor$", Sample), "Tumor", "Normal")) %>%
  # 转为二分类标签（1=Tumor, 0=Normal）
  mutate(Label = ifelse(Group == "Tumor", 1, 0))

# 3. 定义画ROC曲线的函数（和你示例风格一致）
plot_roc <- function(gene_name, data) {
  # 筛选当前基因数据
  gene_data <- data[V1 == gene_name]
  
  # 计算ROC和AUC
  roc_obj <- roc(gene_data$Label, gene_data$Expression)
  auc_val <- round(auc(roc_obj), 3)
  ci_obj <- ci.auc(roc_obj)
  ci_lower <- round(ci_obj[1], 3)
  ci_upper <- round(ci_obj[3], 3)
  
  # 提取坐标点
  roc_df <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities
  )
  
  # 文件名
  pdf_file <- paste0(out_dir, gene_name, "_ROC.pdf")
  
  # 保存 PDF（关键修复）
  pdf(pdf_file, width = 5, height = 5)
  # 画图（和你示例风格一致）
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#1f77b4", linewidth = 1.2) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = "#1f77b4", alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste0(gene_name, " Diagnostic"),
      x = "1-Specificity (FPR)",
      y = "Sensitivity (TPR)"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "gray90"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    annotate("text", x = 0.7, y = 0.3, 
             label = paste0("AUC: ", auc_val, "\nCI: ", ci_lower, "-", ci_upper),
             size = 5, hjust = 0)
  
  print(p)
  dev.off()
  
  # 返回AUC和CI信息
  return(list(gene = gene_name, auc = auc_val, ci = c(ci_lower, ci_upper)))
}

# 4. 分别画三个基因的ROC曲线
genes <- c("IRF1", "Enh092845", "CLEC11A")
results <- lapply(genes, function(g) plot_roc(g, expr_sub))

# 5. 输出汇总结果
cat("=== 三个基因诊断效能汇总 ===\n")
for (res in results) {
  cat(sprintf("%s: AUC=%.3f, 95%% CI=%.3f-%.3f\n", 
              res$gene, res$auc, res$ci[1], res$ci[2]))
}



# 三个基因表达量小题图

# 加载需要的包
library(ggplot2)
library(dplyr)
library(data.table)
library(ggsignif)

# 1. 转为 data.table
setDT(Expr)

# 2. 提取三个基因的表达 + 分组
expr_sub <- Expr[V1 %in% c("IRF1", "Enh092845", "CLEC11A")] %>%
  melt(id.vars = "V1", variable.name = "Sample", value.name = "Expression") %>%
  # 从样本名提取分组
  mutate(Group = ifelse(grepl("_Tumor$", Sample), "Tumor", "Normal")) %>%
  # 转为因子，保证顺序固定
  mutate(Group = factor(Group, levels = c("Tumor", "Normal")))

# 3. 定义画小提琴图的函数（和示例风格一致）
plot_violin <- function(gene_name, data, out_dir) {
  # 筛选当前基因数据
  gene_data <- data[V1 == gene_name]
  
  # 保存文件名
  pdf_file <- paste0(out_dir, gene_name, "_violin.pdf")
  
  # 画图
  p <- ggplot(gene_data, aes(x = Group, y = Expression, fill = Group)) +
    # 小提琴图
    geom_violin(trim = FALSE, alpha = 0.7) +
    # 箱线图（窄版，放在中间）
    geom_boxplot(width = 0.2, fill = "black", color = "black", outlier.shape = NA) +
    # 散点（可选，你示例里也有散点）
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
    # 显著性标记（Tumor vs Normal）
    geom_signif(comparisons = list(c("Tumor", "Normal")),
                map_signif_level = TRUE,
                y_position = max(gene_data$Expression, na.rm = TRUE) * 1.1,
                annotations = "***",
                tip_length = 0.01) +
    # 颜色（和你示例的蓝/红一致）
    scale_fill_manual(values = c("Tumor" = "#4169E1", "Normal" = "#DC143C")) +
    labs(
      x = "",
      y = paste0(gene_name, " expression")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none"
    )
  
  # 保存为PDF
  pdf(pdf_file, width = 4, height = 5)
  print(p)
  dev.off()
  
  cat("✅ 已保存：", pdf_file, "\n")
  
  # 返回统计结果
  wilcox_test <- wilcox.test(Expression ~ Group, data = gene_data)
  return(list(gene = gene_name, p.value = wilcox_test$p.value))
}

# 4. 批量画图
genes <- c("IRF1", "Enh092845", "CLEC11A")
results <- lapply(genes, function(g) plot_violin(g, expr_sub, out_dir))

# 5. 输出汇总结果
cat("\n=== 三个基因差异表达 Wilcoxon 检验结果 ===\n")
for (res in results) {
  cat(sprintf("%s: p-value = %.2e\n", res$gene, res$p.value))
}

cat("\n🎉 全部完成！图片保存在：", out_dir, "\n")



