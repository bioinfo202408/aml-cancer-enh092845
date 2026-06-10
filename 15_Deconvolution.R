# ==============================================================
# 步骤0：环境初始化 & 路径设置
# ==============================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# 主工作目录
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("初始工作路径：", getwd(), "\n")

# 加载全局环境配置
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建并进入反卷积结果目录
out_folder <- "./15_反卷积/"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("当前输出目录：", getwd(), "\n")

# ==============================================================
# 步骤1：安装并加载所有依赖包
# ==============================================================
# if (!require(devtools)) install.packages("devtools")
# # 安装EPIC反卷积包
# devtools::install_github("GfellerLab/EPIC", upgrade = "never")

library(EPIC)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(qs)
library(Seurat)

# ==============================================================
# 步骤2：下载并预处理 TCGA-LAML Bulk RNA-seq 数据
# ==============================================================
## 2.1 下载TCGA-AML数据
query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(query)
data_se <- GDCprepare(query)

## 2.2 FPKM 转 TPM + log2(TPM+1)
fpkm <- assay(data_se, "fpkm_unstrand")

fpkm2tpm <- function(fpkm_mat){
  total_per_sample <- colSums(fpkm_mat)
  tpm_mat <- t(t(fpkm_mat) / total_per_sample) * 1e6
  return(tpm_mat)
}

tpm <- fpkm2tpm(fpkm)
# 过滤低表达基因
keep <- rowSums(tpm > 0) > 0.1 * ncol(tpm)
tpm_filter <- tpm[keep, ]
bulk_expr <- as.matrix(tpm_filter)
bulk_expr_log <- log2(bulk_expr + 1)

## 2.3 匹配AMLFinder风险分组样本
risk_df <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")
bulk_sample_short <- substr(colnames(bulk_expr_log), 1, 12)
shared_sample <- intersect(bulk_sample_short, risk_df$sample)
idx <- which(bulk_sample_short %in% shared_sample)
bulk_expr_filter <- bulk_expr_log[, idx]

## 2.4 Ensembl ID 转换为 Gene Symbol
ens_id <- sub("\\..*", "", rownames(bulk_expr_filter))
symbol <- mapIds(
  org.Hs.eg.db,
  keys = ens_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 去除无注释基因
bulk_symbol <- bulk_expr_filter[!is.na(symbol), ]
rownames(bulk_symbol) <- symbol[!is.na(symbol)]

# ==============================================================
# 步骤3：单细胞数据处理：Marker基因 + 亚群平均表达参考谱
# ==============================================================
## 3.1 读取单细胞Seurat对象
sc_obj <- qread("./GSE116256_单核细胞分高低组亚群+其它亚群.qs")
table(sc_obj$celltype)

## 3.3 筛选高置信Marker
marker_raw <- read.csv("./AML_单细胞不同亚群_marker基因.csv")
marker_sel <- marker_raw %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  dplyr::select(cluster, gene) %>%   # 明确使用dplyr的select
  distinct()
sig_gene_all <- unique(marker_sel$gene)

## 3.4 基因交集（Bulk & Marker）
gene_intersect <- intersect(rownames(bulk_symbol), sig_gene_all)
bulk_sig <- bulk_symbol[gene_intersect, ] %>% as.matrix()

## 3.5 计算单细胞各亚群平均表达矩阵（EPIC参考谱）
sc_expr <- GetAssayData(sc_obj, assay = "RNA", slot = "data")
sc_meta <- data.frame(
  cell = colnames(sc_expr),
  celltype = sc_obj$celltype,
  stringsAsFactors = F
)

# 按细胞类型求均值
sc_expr_t <- as.data.frame(t(sc_expr))
sc_expr_t$celltype <- sc_meta$celltype

sc_mean <- sc_expr_t %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  ungroup()

# 重塑为 EPIC 标准矩阵：行=基因，列=细胞亚群
sc_ref_mat <- sc_mean %>%
  pivot_longer(-celltype, names_to = "gene", values_to = "mean_expr") %>%
  pivot_wider(names_from = celltype, values_from = mean_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

write.csv(as.data.frame(sc_ref_mat), "./AML_单细胞不同亚群_marker基因_的表达矩阵.csv", row.names = T)

# ==============================================================
# 步骤4：基因对齐 + EPIC 自定义反卷积
# ==============================================================
# 三重基因对齐：bulk / 单细胞参考 / Marker
final_gene <- intersect(rownames(bulk_sig), rownames(sc_ref_mat))
bulk_epic_in <- bulk_sig[final_gene, ]
sc_ref_epic_in <- sc_ref_mat[final_gene, ]

# 组装EPIC自定义参考列表
epic_ref <- list(
  refProfiles = sc_ref_epic_in,
  sigGenes = final_gene
)

# 运行反卷积
epic_res <- EPIC(bulk = bulk_epic_in, reference = epic_ref)

# 提取细胞比例结果
cell_frac <- epic_res$cellFractions %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")
write.csv(cell_frac, "./反卷积后结果.csv", row.names = F)

# 查看迭代收敛情况
table(epic_res$fit.gof$convergeCode)

# ==============================================================
# 步骤5：合并风险分组 + 统计检验 + 箱线图可视化
# ==============================================================
# 样本ID匹配
cell_frac$sample_short <- substr(cell_frac$sample_id, 1, 12)
merge_df <- left_join(
  cell_frac,
  risk_df[, c("sample", "risk")],
  by = c("sample_short" = "sample")
)

# 数值型risk转为分类因子
merge_df$risk <- factor(merge_df$risk, levels = c(0,1), labels = c("Low","High"))

# 目标亚群：mono_cell_high
target <- "mono_cell_low"

# Wilcoxon 秩和检验
form <- as.formula(paste0("`", target, "` ~ risk"))
wt_test <- wilcox.test(form, data = merge_df)
cat("\n===== ", target, " 组间检验结果 =====\n")
print(wt_test)

# Wilcoxon rank sum test with continuity correction
# 
# data:  mono_cell_low by risk
# W = 836, p-value = 8.732e-06
# alternative hypothesis: true location shift is not equal to 0

# # 绘制箱线图并保存PDF（服务器无图形界面，只存文件）
# p <- ggplot(merge_df, aes(x = risk, y = .data[[target]], fill = risk)) +
#   geom_boxplot(outlier.shape = NA, width = 0.6) +
#   geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
#   labs(
#     x = "AMLFinder Risk",
#     y = paste0(target, " Cell Fraction"),
#     title = paste0(target, " proportion between risk groups")
#   ) +
#   scale_fill_manual(values = c("Low" = "#6699CC", "High" = "#E25822")) +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# # 文件名规避特殊字符
# pdf_name <- paste0(gsub("/", "_", target), "_fraction_boxplot_mono_cell_low.pdf")
# ggsave(pdf_name, plot = p, width = 6, height = 5, dpi = 300, device = "pdf")
# 


#  画图
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)
library(rlang)

# 路径设置
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/15_反卷积/")
out_dir <- "./"

# 读取数据
cell_frac <- fread("./反卷积后结果_最新.csv", data.table = F)
risk <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv", row.names = NULL)

# 样本ID匹配
cell_frac$sample_short <- substr(cell_frac$sample_id, 1, 12)
merge_df <- dplyr::left_join(
  cell_frac,
  risk[, c("sample", "risk")],
  by = c("sample_short" = "sample")
)

# 分组命名
merge_df$Group <- ifelse(merge_df$risk == 0, "Low Risk", "High Risk")

# 目标细胞
# target_gene <- "mono_cell_low"
target_gene <- "mono_cell_high"


# 关键修复：全部加 dplyr:: 规避冲突，!!sym 兼容列名
plot_df <- merge_df %>% 
  dplyr::select(sample_short, Group, !!sym(target_gene)) %>%
  dplyr::rename(Expression = !!sym(target_gene))

# 过滤极值（保留中间80%）
# plot_df_clean <- plot_df %>%
#   dplyr::group_by(Group) %>%
#   dplyr::mutate(
#     q10 = quantile(Expression, 0.2, na.rm = TRUE),
#     q90 = quantile(Expression, 0.8, na.rm = TRUE)
#   ) %>%
#   dplyr::filter(Expression >= q10 & Expression <= q90) %>%
#   dplyr::ungroup()
plot_df_clean <- plot_df 
# 关键改动：指定x轴分组顺序，让High Risk在右侧
plot_df_clean$Group <- factor(plot_df_clean$Group,
                              levels = c("Low Risk", "High Risk"))

# 绘图参数
my_comparisons <- list(c("Low Risk", "High Risk"))
group_colors <- c("High Risk" = "#AD3D3E", "Low Risk" = "#3A3E96")

# 半小提琴云雨图
p <- ggplot(plot_df_clean, aes(x = Group, y = Expression, color = Group)) +
  geom_half_violin(
    position = position_nudge(x = 0.15),
    side = 'r', trim = F, alpha = 0.8, width = 0.5, linewidth = 1
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black") +
  geom_jitter(size = 1.8, alpha = 0.5, width = 0.15) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0,
    label.y.npc = "top",
    size = 5
  ) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = paste0(target_gene, " Cell Fraction")) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=15, face="bold")
  )

# 保存图片
fname <- paste0(gsub("[ /]", "_", target_gene), "_RiskGroup_NC_云雨图.pdf")
ggsave(fname, plot = p, width = 4.5, height = 6)




# 1. 计算 mono_cell_high 的中位数
median_mono <- median(KM$mono_cell_high, na.rm = TRUE)
cat("mono_cell_high 中位数为：", median_mono, "\n")

# 2. 按中位数分组
KM$mono_group <- ifelse(KM$mono_cell_high >= median_mono, "mono_high", "mono_low")

# 3. 查看分组结果
table(KM$mono_group)
head(KM[, c("sample_id", "mono_cell_high", "mono_group")])




# 加载包（与你原有代码保持一致）
library(survival)
library(ggsurvfit)
library(patchwork)
library(dplyr)
library(data.table)
library(maxstat)
library(scales)
library(survminer)

# ===================== 1. 数据合并：KM(细胞比例) + risk_df(生存信息) =====================
# 按 sample_id 匹配生存、状态、时间
# 根据mono_cell_high 中位数 评分将患者划分为高低评分两组，查看两组的KM曲线差异
KM <- read.csv("./反卷积后结果_最新.csv")
risk_df <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")

KM$sample_short_id <- substr(KM$sample_id, 1, 12)
# 按截断ID匹配
merge_surv <- KM %>%
  dplyr::left_join(risk_df[, c("sample", "OS", "OS.time")], 
                   by = c("sample_short_id" = "sample")) %>%
  dplyr::filter(!is.na(OS.time), !is.na(OS))

# ===================== 2. 基于 mono_cell_high 中位数分组 =====================
# 计算中位数
median_cut <- median(merge_surv$mono_cell_high, na.rm = TRUE)
cat("mono_cell_high 中位数截断值：", median_cut, "\n")

# 分组：高于中位数=High risk，低于=Low risk（和你原图分组逻辑对齐）
merge_surv$risk_group <- ifelse(merge_surv$mono_cell_high > median_cut, "High", "Low")

# 因子顺序 & 标签（完全沿用你原有样式）
merge_surv$risk_group <- factor(
  merge_surv$risk_group,
  levels = c("Low", "High"),
  labels = c("Low risk", "High risk")
)

# ===================== 3. 拟合生存模型、提取P值、HR、置信区间 =====================
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = merge_surv)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = merge_surv)

# 提取HR、95%CI、Log-rank P值
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low  <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# 提取每组样本数、事件数（备用）
tab_surv <- summary(fit)
n_low  <- tab_surv$n[1]
n_high <- tab_surv$n[2]

# ===================== 4. 配色（和你原有配色完全一致） =====================
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# ===================== 5. 绘图（1:1复刻你原有KM图样式） =====================
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # 风险人数表
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # 删失标记
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P值标注
  add_pvalue(
    location = "annotation",
    x = max(merge_surv$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR 及95%CI标注
  annotate("text",
           x = max(merge_surv$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ", hr_val, " (95%CI: ", hr_low, "-", hr_high, ")"),
           hjust = 1, size = 4.2) +
  # 坐标轴标题
  labs(
    title = "mono_cell_high - Overall Survival",
    x = "Time (months)",
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03, 0)) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # 主题完全复用
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3, 0.8, 0.3, 0.3, unit = "cm")
  )

# 上下布局（生存曲线 + 风险表）
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# ===================== 6. 保存图片 =====================
ggsave(
  ".//mono_cell_high_Median_KM.pdf",
  plot = p,
  width = 7, height = 5.5
)

# 查看分组统计
table(merge_surv$risk_group)



