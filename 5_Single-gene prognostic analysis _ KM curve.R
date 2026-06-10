# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./5_单基因预后分析_KM曲线/"
dir.create(out_dir, showWarnings = F, recursive = T)

# 加载包
library(survival)
library(survminer)
library(dplyr)
library(data.table)


# 导入数据
GSE1656 <- fread("/home/weili/Project/AML/snakefile/gene_data/homo/all_GSE165656_clean_expr.csv")
# write.csv(GSE1656, "../TCGA_data/GSE165656_expr_all_RNA.csv")
GSE1656_time <- read.csv("../TCGA_data/GSE165656_time_status.csv")
ID_SRR <- read.table("../TCGA_data/GSE165656_ID和SRR对应关系.txt")

# ===================== 1. 读入你的数据 =====================
# 表达矩阵 (行：基因，列：SRR样本)
expr <- GSE1656  
rownames(expr) <- expr$V1  # 第一列是基因名，设为行名
expr <- data.frame(expr[,-1], row.names = rownames(expr)) #rownames.force=TRUE 强制保留行名
head(rownames(expr))


# 生存数据 (ID = 样本数字ID, OS.time, OS)
time <- GSE1656_time  

# 对应关系 (SRR <-> 样本ID)
id_srr <- ID_SRR  
colnames(id_srr) <- c("SRR", "ID")

# ===================== 2. 输入你要分析的基因 =====================
gene <- "CLEC11A"  # 例如：gene <- "FOXO1" Enh092845 IRF1

# ===================== 3. 提取该基因表达量 =====================
if(!gene %in% rownames(expr)){
  stop(paste0(gene, " 不在表达矩阵中！"))
}

gene_exp <- data.frame(
  SRR = colnames(expr),
  exp = as.numeric(expr[gene, ])
)

# ===================== 4. 匹配样本ID + 生存数据 =====================
gene_exp <- gene_exp %>%
  left_join(id_srr, by = "SRR") %>%
  left_join(time[, c("ID", "OS.time", "OS")], by = "ID") %>%
  filter(!is.na(OS.time), !is.na(OS))

# # ===================== 5.1 按表达量中位数分组高低风险 =====================
# gene_exp$risk_group <- ifelse(
#   gene_exp$exp > median(gene_exp$exp, na.rm=T), 
#   "High", "Low"
# )

# ===================== 5.2 按最优截断值分组高低风险 =====================
library(maxstat)
# 自动找最优cutoff
max_stat <- maxstat.test(Surv(OS.time, OS) ~ exp, data = gene_exp, smethod = "LogRank", pmethod = "HL")
cutoff <- max_stat$estimate

# 分组
gene_exp$risk_group <- ifelse(gene_exp$exp > cutoff, "High", "Low")
cat(gene, "最优截断值 =", cutoff, "\n")
# 
# # ===================== 5.3 按均值分组高低风险 =====================
# gene_exp$risk_group <- ifelse(
#   gene_exp$exp > mean(gene_exp$exp, na.rm=T), 
#   "High", "Low"
# )
# 
# # ===================== 5.4 按四分位数分组高低风险 =====================
# # 只保留 最高 25% + 最低 25%，中间 50% 去掉，差异通常最大。
# q1 <- quantile(gene_exp$exp, 0.25, na.rm=T)
# q3 <- quantile(gene_exp$exp, 0.75, na.rm=T)
# 
# gene_exp <- gene_exp[gene_exp$exp <= q1 | gene_exp$exp >= q3, ]
# gene_exp$risk_group <- ifelse(gene_exp$exp >= q3, "High", "Low")


# ===================== 6. 生存分析与绘图 =====================
# 1) 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = gene_exp)

# 2) 发表级主题
pub_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=0.8),
        plot.title = element_text(hjust=0.5, size=14),
        axis.text = element_text(size=11, color="black"))

pdf("./5_单基因预后分析_KM曲线/CLEC11A_KM_curve_cutoff.pdf", width=6, height=5)
# 3) 画KM图
ggsurvplot(
  fit,
  data = gene_exp,
  title = paste0(gene, " - Overall Survival"),
  xlab = "Time (months)",
  ylab = "Survival probability",
  risk.table = TRUE,
  pval = TRUE,
  conf.int = F,
  palette = c("#2E86AB", "#E63946"),
  legend.labs = c(paste0(gene," Low"), paste0(gene," High")),
  ggtheme = pub_theme,
  surv.median.line = "hv"
)
dev.off()
