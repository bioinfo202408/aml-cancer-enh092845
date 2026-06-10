# 参考：https://mp.weixin.qq.com/s/pORXLeLI9O_WS06IZTIa-w   半云雨+柱状图
# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./6_小提琴图_半云雨图//"
dir.create(out_dir, showWarnings = F, recursive = T)

# ========================
# 加载包
# ========================
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)  # 🔥 画半小提琴关键包

# ===============================================
###关键基因正常/患病表达    半小提琴云雨图 #######
# ===============================================

# ========================
# 读数据（你的原样）
# ========================
allRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
rownames(allRNA_expr) <- allRNA_expr$V1

group_info <- read.csv("../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv", row.names=1)
group_info$Group <- ifelse(group_info$Group == 1,"Tumor","Normal")

# ========================
# 正确提取表达矩阵（修复版）
# ========================
expr_mat <- as.data.frame(allRNA_expr)
rownames(expr_mat) <- expr_mat$V1
expr_mat <- expr_mat[, -1]  # 去掉V1，只剩表达量

# ========================
# 提取3个基因
# ========================
# genes <- c("Enh092845", "CLEC11A", "IRF1")
# genes <- c("Enh092845")
genes <- c("CLEC11A")
# genes <- c("IRF1")
expr_sub <- expr_mat[genes, ]  # 直接提取
expr_sub <- expr_sub[,-1]


# ========================
# 正确构建绘图矩阵（关键！）
# ========================
plot_df <- as.data.frame(t(expr_sub))
plot_df$sample <- rownames(plot_df)
plot_df$Group <- group_info[plot_df$sample, "Group"]

# 长格式转换
plot_df <- pivot_longer(
  plot_df,
  cols = all_of(genes),
  names_to = "Gene",
  values_to = "Expression"
)

# ========================
# 表达量负值校正（必须）+log2处理
# ========================
plot_df <- plot_df %>%
  group_by(Gene) %>%
  mutate(
    min_exp = min(Expression, na.rm=T),
    Expression = ifelse(min_exp < 0, Expression - min_exp, Expression)
  ) %>%
  ungroup()
plot_df$Expression <- log2(plot_df$Expression+1)

# ======================
# 🔥 关键：剔除每组 top/bottom 10% 极端值
# ======================
plot_df_clean <- plot_df %>%
  group_by(Group, Gene) %>%  # 按分组+基因分别剔除
  mutate(
    q10 = quantile(Expression, 0.3, na.rm = T),
    q90 = quantile(Expression, 0.7, na.rm = T)
  ) %>%
  filter(Expression >= q10 & Expression <= q90) %>%  # 保留中间80%
  ungroup()


# ========================
# 两组比较设置
# ========================
my_comparisons <- list(c("Normal", "Tumor"))

# ========================
# 配色（你参考代码的高级风格）
# ========================
group_colors <- c(
  "Tumor"  = "#AD3D3E",    # 红色
  "Normal" = "#3A3E96"     # 蓝色
)
# '#3A3E96', '#AD3D3E', '#50A293','#E8B75E'


# # ========================
# # ✅ 绘制高质量小提琴图（
# # ========================
# p <- ggplot(plot_df_clean, aes(x = Group, y = Expression, fill = Group, color = Group)) +
#   geom_violin(trim = FALSE, width = 0.9, alpha = 0.5, linewidth = 0.5) +
#   geom_boxplot(
#     width = 0.16,
#     outlier.shape = NA,
#     fill = "white",
#     color = "black",
#     linewidth = 0.4
#   ) +
#   geom_jitter(
#     width = 0.12,
#     size = 1.4,
#     alpha = 0.5,
#     stroke = 0
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons,
#     method = "wilcox.test",
#     label = "p.signif",
#     size = 4
#   ) +
#   scale_fill_manual(values = group_colors) +
#   scale_color_manual(values = group_colors) +
#   labs(
#     title = genes,
#     x = NULL,
#     y = paste0(genes, " Expression log2(TPM+1)")
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#     axis.title.y = element_text(face = "bold", size = 14),
#     axis.text.x = element_text(face = "bold", color = "black", size = 12),
#     axis.text.y = element_text(color = "black", size = 11),
#     legend.position = "none",
#     axis.line = element_line(linewidth = 0.6, color = "black")
#   )
# 
# # ========================
# # 保存图片
# # ========================
# ggsave(
#   filename = paste0(out_dir, genes, "_violin_Tumor_vs_Normal_FWSE数据集_new.pdf"),
#   plot = p,
#   width = 4.8,
#   height = 4.8
# )


# ========================
# ✨ 绘制 NC 风格 半小提琴云雨图
# ========================
p <- ggplot(plot_df_clean, aes(x = Group, y = Expression, color = Group)) +
  # 1. 半小提琴图（右侧）
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  # 2. 箱线图
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1
  ) +
  # 3. 均值点
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  # 4. 抖动散点（雨）
  geom_jitter(
    size = 1.8, alpha = 0.5, width = 0.15
  ) +
  # 5. 显著性检验
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0,
    label.y.npc = "top",
    size = 5
  ) +
  # 6. 颜色
  scale_color_manual(values = group_colors) +
  # 7. 标签
  labs(
    x = NULL,
    y = paste0(genes, " Expression log2(TPM+1)")
  ) +
  # 8. 主题（完全 NC 风格）
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size=14),
    axis.text.x = element_text(angle=0, hjust=0.5, size=14, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    strip.background = element_rect(fill = "transparent"),
    plot.title = element_text(hjust=0.5, size=17, face="bold")
  )

# ========================
# 保存（cairo_pdf 保证符号清晰）
# ========================
ggsave(
  filename = paste0(out_dir, genes, "_NC_云雨图_Tumor_vs_Normal.pdf"),
  plot = p,
  width = 4.5,
  height = 6
)








# ===============================================================================================================
# #### 预后模型不同分组（年龄、性别、细胞遗传学评分，高低风险组的 risk source 半小提琴云雨图 ####
# ===============================================================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix) 
library(ggpubr)
library(gghalves)  # 🔥 必须加载

risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")
clina <- read.csv("../TCGA_data/TCGA_clina_ostime.csv")
clina2 <- read.csv("../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime_new_136_samples.csv")

risk$eln <- clina2$eln[match(risk$sample, clina2$ID)] 
risk$age <- clina$age[match(risk$sample, clina$sample)] 
risk$gender <- clina$gender[match(risk$sample, clina$sample)]

risk <- risk[,c("sample", "OS", "OS.time", "riskscore", "risk", "eln", "age", "gender")]
risk <- risk[risk$eln != "'--", ] 



# ====================== 数据预处理
min_risk <- min(risk$riskscore, na.rm = TRUE)  
risk$riskscore_pos <- risk$riskscore + abs(min_risk)

risk <- risk %>%
  mutate(
    age_group_label = factor(ifelse(age >= 65, ">=65", "<65"), levels = c("<65", ">=65")),
    gender_group_label = factor(gender, levels = c("female", "male"), labels = c("Female", "Male")),
    eln_group_label = factor(eln),
    risk_group_label = factor(ifelse(risk == 0, "Low risk", "High risk"), 
                              levels = c("Low risk", "High risk"))
  )

# ====================== 【关键】剔除每组极端值（上下10%）
filter_extreme <- function(df, group_col, val_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      q10 = quantile(.data[[val_col]], 0.1, na.rm=T),
      q90 = quantile(.data[[val_col]], 0.9, na.rm=T)
    ) %>%
    filter(.data[[val_col]] >= q10 & .data[[val_col]] <= q90) %>%
    ungroup()
}

# ====================== 绘图配色（顶刊红蓝）
colors2 <- c("Low risk"="#3A3E96", "High risk"="#AD3D3E",
             "<65"="#3A3E96", ">=65"="#AD3D3E",
             "Female"="#3A3E96", "Male"="#AD3D3E")


# ======================================================================================
# 🔥 核心：替换成 NC 风格 半小提琴云雨图函数（完全兼容你原来的调用方式）
# ======================================================================================
plot_nc_halfviolin <- function(data, x, y, title, colors) {
  
  # 剔除极端值
  data_clean <- filter_extreme(data, x, y)
  
  ggplot(data_clean, aes(x = .data[[x]], y = .data[[y]], color = .data[[x]])) +
    # 1. 半小提琴（右侧）
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth = 1
    ) +
    # 2. 箱线图
    geom_boxplot(
      outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill = "white"
    ) +
    # 3. 均值点
    stat_summary(
      fun = "mean", geom = "point", shape = 20, size = 3,
      color = "black", fill = "black", alpha = 0.8
    ) +
    # 4. 散点（雨）
    geom_jitter(
      size = 1.8, alpha = 0.5, width = 0.15
    ) +
    # 5. 显著性检验
    stat_compare_means(
      method = "wilcox.test", label = "p.signif",
      tip.length = 0, label.y.npc = "top", size = 5
    ) +
    # 颜色
    scale_color_manual(values = colors) +
    # 标签
    labs(x = "", y = "Risk Score", title = title) +
    # NC 主题
    theme_bw(base_size = 16) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
      strip.background = element_rect(fill = "transparent")
    )
}

# ====================== 绘制 4 张图（调用方式完全不变！）
# 1. 年龄
p1 <- plot_nc_halfviolin(risk, "age_group_label", "riskscore_pos", "Risk by Age Group", colors2)
ggsave("./6_小提琴图/Risk_Age_NC_半云雨图.pdf", p1, width=4.5, height=6)

# 2. 性别
p2 <- plot_nc_halfviolin(risk, "gender_group_label", "riskscore_pos", "Risk by Gender Group", colors2)
ggsave("./6_小提琴图/Risk_Gender_NC_半云雨图.pdf", p2, width=4.5, height=6)

# 3. 高低风险组
p3 <- plot_nc_halfviolin(risk, "risk_group_label", "riskscore_pos", "Risk Group", colors2)
ggsave("./6_小提琴图/Risk_RiskGroup_NC_半云雨图.pdf", p3, width=4.5, height=6)



# 4. ELN 分组（NC 半云雨图 多组版）
# 1. 先把 ELN 分组按你要的顺序重新排序
risk$eln_group_label <- factor(
  risk$eln,
  levels = c("Favorable", "Intermediate", "Adverse")
)

# 2. 配色按顺序定义
colors3 <- c("Favorable"="#50A293", "Intermediate"="#3A3E96", "Adverse"="#AD3D3E")

# 3. 多组显著性比较
eln_comparisons <- list(
  c("Favorable", "Intermediate"),
  c("Favorable", "Adverse"),
  c("Intermediate", "Adverse")
)

# 4. 绘制 NC 风格半云雨图
p4 <- ggplot(risk, aes(x = eln_group_label, y = riskscore_pos, color = eln_group_label)) +
  # 半小提琴
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  # 箱线图
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill="white"
  ) +
  # 均值点
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  # 散点
  geom_jitter(size = 1.8, alpha = 0.5, width = 0.15) +
  # 多组显著性检验
  stat_compare_means(
    comparisons = eln_comparisons,
    method = "wilcox.test", label = "p.signif",
    tip.length = 0, label.y.npc = "top", size = 4.5
  ) +
  scale_color_manual(values = colors3) +
  labs(x = "", y = "Risk Score", title = "Risk by Eln Group") +
  # NC 主题
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size=13),
    axis.text.x = element_text(angle=0, hjust=0.5, size=13, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    plot.title = element_text(hjust=0.5, size=17, face="bold"),
    strip.background = element_rect(fill = "transparent")
  )

# 保存
ggsave("./6_小提琴图/Risk_ELN_NC_半云雨图.pdf", p4, width=5.5, height=6)




# ===============================================================================================================
# #### 拟时序高低分组 → NC 半云雨图版本 ####
# ===============================================================================================================

# 
# ### 拟时序高低分组 → NC 半云雨图版本
# 
load(("/home/weili/Project/AML/human/AML_combined_analyse//scRNA_analyse/20260407分析结果_GSE1116256//rawdata/Step10_拟时序分析_单核细胞高低组亚群_完成.Rdata"))
# 4.4 拟时序箱线图（含显著性）
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gghalves) # 🔥 必须

input.data = data.frame(group = HSMM$group, Pseudotime = HSMM$Pseudotime)

# ======================
# 1. 去除极端值
# ======================
plot_df_clean <- input.data %>%
  group_by(group) %>%
  mutate(
    q10 = quantile(Pseudotime, 0.3, na.rm = T),
    q90 = quantile(Pseudotime, 0.7, na.rm = T)
  ) %>%
  filter(Pseudotime >= q10 & Pseudotime <= q90) %>%
  ungroup()

# ======================
# 2. 两组比较
# ======================
my_comparisons <- list(c("mono_cell_low", "mono_cell_high"))

# ======================
# 3. ✨ 绘制 NC 风格 半小提琴云雨图
# ======================
p <- ggplot(plot_df_clean, aes(x = group, y = Pseudotime, color = group)) +
  # 半小提琴
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  # 箱线图
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill="white"
  ) +
  # 均值点
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  # 散点
  geom_jitter(size = 1.8, alpha = 0.5, width = 0.15) +
  # 显著性
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", label = "p.signif",
    tip.length = 0, label.y.npc = "top", size = 5
  ) +
  # 顶刊配色
  scale_color_manual(values = c("mono_cell_high"="#AD3D3E", "mono_cell_low"="#3A3E96")) +
  # 标签
  labs(x = "", y = "Pseudotime") +
  # NC 主题
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size=14),
    axis.text.x = element_text(angle=0, hjust=0.5, size=14, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    plot.title = element_text(hjust=0.5, size=17, face="bold")
  )

# ======================
# 4. 保存（高清 PDF）
# ======================
ggsave(
  filename = "./6_小提琴图/Pseudotime_NC_半云雨图.pdf",
  plot = p,
  width = 4.5,
  height = 6
)




# ===============================================================================================================
# #### 基质评分、免疫评分、估计分数和肿瘤纯度: NC 半云雨图版本 ####
# ===============================================================================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(gghalves)
library(ggpubr)

# 读数据
core <- read.csv("../Outdata/12.免疫分析/8_estimate_scores_risk_group.csv")

# 数据预处理
scores_df <- core
scores_df$risk <- factor(
  scores_df$risk, 
  levels = c(0, 1),
  labels = c("Low risk", "High risk")
)

# 极端值过滤
filter_extreme <- function(df, group_col, val_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      q10 = quantile(.data[[val_col]], 0.1, na.r=T),
      q90 = quantile(.data[[val_col]], 0.9, na.r=T)
    ) %>%
    filter(.data[[val_col]] >= q10 & .data[[val_col]] <= q90) %>%
    ungroup()
}

# 统一配色（和你全文一致！）
colors2 <- c("Low risk"="#3A3E96", "High risk"="#AD3D3E")
my_comparisons <- list(c("Low risk", "High risk"))

# ===================== 绘图函数：NC 半云雨图（完全统一风格）
plot_estimate_nc <- function(score_col, title, ylab) {
  data_clean <- filter_extreme(scores_df, "risk", score_col)
  
  ggplot(data_clean, aes(x = risk, y = .data[[score_col]], color = risk)) +
    geom_half_violin(
      position = position_nudge(x=0.15, y=0), side='r', 
      trim=F, alpha=0.8, width=0.5, linewidth=1
    ) +
    geom_boxplot(
      outlier.shape=NA, width=0.35, alpha=0.8, linewidth=1, fill="white"
    ) +
    stat_summary(
      fun="mean", geom="point", shape=20, size=2.5,
      color="black", fill="black", alpha=0.8
    ) +
    geom_jitter(size=1, alpha=0.5, width=0.12) +
    stat_compare_means(
      comparisons = my_comparisons,
      method="wilcox.test", label="p.signif",
      tip.length=0, label.y.npc="top", size=3.5
    ) +
    scale_color_manual(values=colors2) +
    labs(title=title, x="", y=ylab) +
    theme_bw(base_size=12) +
    theme(
      plot.title=element_text(hjust=0.5, size=13, face="bold"),
      axis.title.y=element_text(size=11),
      axis.text=element_text(color="black", size=10),
      axis.text.x=element_text(face="bold", size=9),
      legend.position="none",
      panel.grid=element_blank()
    )
}

# 绘制 4 张图
p_stromal <- plot_estimate_nc("StromalScore", "Stromal Score", "Stromal Score")
p_immune  <- plot_estimate_nc("ImmuneScore", "Immune Score", "Immune Score")
p_estima  <- plot_estimate_nc("ESTIMATEScore", "ESTIMATE Score", "ESTIMATE Score")
p_purity  <- plot_estimate_nc("TumorPurity", "Tumor Purity", "Tumor Purity")


# 保存
ggsave(
  "./6_小提琴图//ESTIMATE_stromal_半云雨图.pdf",
  plot = p_stromal,
  width = 4.5,
  height = 6
)

# 保存
ggsave(
  "./6_小提琴图//ESTIMATE_immune_半云雨图.pdf",
  plot = p_immune,
  width = 4.5,
  height = 6
)

# 保存
ggsave(
  "./6_小提琴图//ESTIMATE_estima_半云雨图.pdf",
  plot = p_estima,
  width = 4.5,
  height = 6
)

# 保存
ggsave(
  "./6_小提琴图//ESTIMATE_purity_半云雨图.pdf",
  plot = p_purity,
  width = 4.5,
  height = 6
)





# ===============================================================================================================
# ####CLEC11A在不同细胞亚群表达量 ####
# ===============================================================================================================

load("../scRNA_analyse//20260407分析结果_GSE1116256//rawdata/Step8.4.CLEC11A在不同细胞亚群表达量.Rdata")
# 1. 加载必须包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggpubr)
library(gghalves)      # 画半小提琴必需
library(RColorBrewer)

# 2. 提取每个细胞的 CLEC11A 表达 + 细胞类型
plot_df <- data.frame(
  celltype = seurat.data2$celltype_new,
  Expression = seurat.data2@assays$RNA@data["CLEC11A", ]
)

# 3. 表达量预处理（负值校正 + log2）
plot_df <- plot_df %>%
  mutate(
    min_exp = min(Expression, na.rm = T),
    Expression = ifelse(min_exp < 0, Expression - min_exp, Expression),
    Expression = log2(Expression + 1)
  )

# 4. 剔除极端值（中间80%，更干净）
plot_df_clean <- plot_df %>%
  group_by(celltype) %>%
  mutate(
    q10 = quantile(Expression, 0.4, na.rm=T),
    q90 = quantile(Expression, 0.6, na.rm=T)
  ) %>%
  filter(Expression >= q10 & Expression <= q90) %>%
  ungroup()

# 5. 按表达均值从高到低排序细胞亚群（关键！）
cell_order <- plot_df_clean %>%
  group_by(celltype) %>%
  summarise(mean_exp = mean(Expression, na.rm=T)) %>%
  arrange(desc(mean_exp)) %>%
  pull(celltype)

plot_df_clean$celltype <- factor(plot_df_clean$celltype, levels = cell_order)

# 6. 配色（高级配色，和你原图风格一致）
n <- length(unique(plot_df_clean$celltype))
cell_colors <- brewer.pal(n, "Set2")  # 你原来的配色方案

# 7. 绘制 NC 半小提琴云雨图（100% 对齐你第一张图）
p <- ggplot(plot_df_clean, aes(x = celltype, y = Expression, color = celltype)) +
  # 半小提琴
  geom_half_violin(
    position = position_nudge(x=0.15, y=0),
    side = 'r', trim=F, alpha=0.8, width=0.5, linewidth=1
  ) +
  # 箱线图
  geom_boxplot(
    outlier.shape = NA, width=0.35, alpha=0.8, linewidth=1
  ) +
  # 均值点
  stat_summary(
    fun = "mean", geom = "point", shape=20, size=3,
    color="black", fill="black", alpha=0.8
  ) +
  # 散点（雨）
  geom_jitter(size=1.8, alpha=0.5, width=0.15) +
  # 配色
  scale_color_manual(values = cell_colors) +
  # 标签
  labs(
    x = "Cell Type",
    y = "CLEC11A Expression log2(TPM+1)",
    title = "CLEC11A Expression Across Cell Subsets"
  ) +
  # NC 主题（完全一致）
  theme_bw(base_size=16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(color="black", size=14),
    axis.text.x = element_text(angle=45, hjust=1, size=14, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    strip.background = element_rect(fill="transparent"),
    plot.title = element_text(hjust=0.5, size=17, face="bold")
  )

# 8. 如果你需要：给 Mono_high vs Mono_low 加显著性
# 提取两组数据
two_group <- plot_df_clean %>%
  filter(celltype %in% c("Mono_cell_high", "Mono_cell_low"))

my_comp <- list(c("Mono_cell_low", "Mono_cell_high"))

p <- p +
  stat_compare_means(
    data = two_group,
    comparisons = my_comp,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0,
    size = 5
  )


ggsave(
  "./6_小提琴图_半云雨图//CLEC11A_细胞亚群_半小提琴云雨图.pdf",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)



# ===============================================================================================================
# ####高低风险组的药物敏感性 ####
# ===============================================================================================================

# ========================
# 加载包
# ========================
library(ggplot2)
library(dplyr)

# ========================
# 数据准备
# ========================
df <- data.frame(read.csv("../Outdata/15.药物敏感性/drug_risk_diff_top10.csv"))
# ========================
# 每个药物分配不同颜色（10种鲜明颜色）
# ========================
drug_colors <- c(
  "X5.Fluorouracil" = "#E63946",   # 红色
  "BI.2536"        = "#F77F00",   # 橙色
  "BMS.754807"     = "#FCBF49",   # 黄色
  "Navitoclax"     = "#06D6A0",   # 绿色
  "ABT737"        = "#118AB2",   # 蓝色
  "Dactolisib"     = "#073B4C",   # 深蓝
  "Daporinad"      = "#9B5DE5",   # 紫色
  "PF.4708671"     = "#F15BB5",   # 粉色
  "Trametinib"     = "#00BBF9",   # 天蓝
  "Cisplatin"      = "#8B5A2B"    # 棕色
)

# ========================
# 按diff排序（从小到大）
# ========================
df <- df %>%
  mutate(Drug = factor(Drug, levels = Drug[order(diff)]))

# ========================
# 绘制棒棒糖图
# ========================
p <- ggplot(df, aes(x = diff, y = Drug)) +
  # 1. 绘制线段（从0到diff）
  geom_segment(
    aes(x = 0, xend = diff, y = Drug, yend = Drug, color = Drug),
    linewidth = 1.5, alpha = 0.7
  ) +
  # 2. 绘制圆点（在diff位置）
  geom_point(
    aes(color = Drug),
    size = 5, stroke = 1.5, fill = "white"
  ) +
  # 3. 颜色映射
  scale_color_manual(values = drug_colors) +
  # 4. 添加零线
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  # 5. 添加显著性标记
  geom_text(
    aes(label = sig, 
        hjust = ifelse(diff >= 0, -0.3, 1.3)),
    size = 4.5, fontface = "bold", color = "black"
  ) +
  # 6. 标签
  labs(
    title = "Drug Sensitivity: High vs Low Expression Group",
    x = "Difference (mean_high - mean_low)",
    y = NULL
  ) +
  # 7. 主题设置
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  # 8. 调整x轴范围，给显著性标记留空间
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))

# ========================
# 保存图片
# ========================
ggsave(
  filename = "./6_小提琴图_半云雨图/高低风险组药物预测_drug_lollipop_plot.pdf",
  plot = p,
  width = 10,
  height = 5
)






# ========================
# 加载包
# ========================
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)  # 半小提琴关键包

# ========================
# 读取数据（长格式）
# ========================
drug_sens_target <- read.csv("../Outdata/15.药物敏感性/7种AML常用药物在高低风险组的差异.csv")

# 确保Risk_Group是因子，顺序为 Low 在前
drug_sens_target$Risk_Group <- factor(drug_sens_target$Risk_Group, levels = c("Low", "High"))

# ========================
# 配色（高低风险组）
# ========================
group_colors <- c(
  "High" = "#AD3D3E",   # 红色 - 高风险
  "Low"  = "#3A3E96"    # 蓝色 - 低风险
)

# ========================
# 输出目录
# ========================
out_dir <- "./6_小提琴图_半云雨图/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========================
# 获取所有药物列表
# ========================
drug_list <- unique(drug_sens_target$Drug)
cat("共有", length(drug_list), "种药物：", paste(drug_list, collapse = ", "), "\n")

# ========================
# 循环绘制每种药物，单独保存
# ========================
for(drug in drug_list) {
  
  # 提取当前药物数据
  df_sub <- drug_sens_target %>% filter(Drug == drug)
  
  # 计算统计值（用于标签）
  mean_low <- mean(df_sub$Sensitivity[df_sub$Risk_Group == "Low"], na.rm = TRUE)
  mean_high <- mean(df_sub$Sensitivity[df_sub$Risk_Group == "High"], na.rm = TRUE)
  diff <- mean_high - mean_low
  
  # Wilcoxon检验
  test_result <- wilcox.test(Sensitivity ~ Risk_Group, data = df_sub)
  p_val <- test_result$p.value
  
  # 显著性标记
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))
  
  cat("\n【", drug, "】", "Low=", round(mean_low, 3), 
      " High=", round(mean_high, 3), 
      " Diff=", round(diff, 3), 
      " p=", format(p_val, digits = 3), sig, "\n")
  
  # ========================
  # 绘制半小提琴云雨图
  # ========================
  p <- ggplot(df_sub, aes(x = Risk_Group, y = Sensitivity, color = Risk_Group)) +
    # 1. 半小提琴图（右侧）
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth = 1
    ) +
    # 2. 箱线图
    geom_boxplot(
      outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1
    ) +
    # 3. 均值点
    stat_summary(
      fun = "mean", geom = "point", shape = 20, size = 3,
      color = "black", fill = "black", alpha = 0.8
    ) +
    # 4. 抖动散点（雨）
    geom_jitter(
      size = 2, alpha = 0.5, width = 0.15
    ) +
    # 5. 显著性检验
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.03,
      label.y.npc = "top",
      size = 5
    ) +
    # 6. 颜色
    scale_color_manual(values = group_colors) +
    # 7. 标签
    labs(
      title = drug,
      subtitle = paste0("Low=", round(mean_low, 3), 
                        " | High=", round(mean_high, 3),
                        " | Diff=", round(diff, 3),
                        " | p", sig),
      x = NULL,
      y = "Drug Sensitivity (IC50)"
    ) +
    # 8. 主题（NC风格）
    theme_bw(base_size = 16) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      panel.border = element_rect(color = "black", linewidth = 1)
    )
  
  # ========================
  # 保存图片
  # ========================
  # PDF
  ggsave(
    filename = paste0(out_dir, drug, "_半云雨图_高低风险组.pdf"),
    plot = p,
    width = 4.5,
    height = 6
  )
  
  
}

# Cytarabine 】 Low= 7.434  High= 9.255  Diff= 1.821  p<0.05  1
# 【 Venetoclax 】 Low= 9.439  High= 9.862  Diff= 0.422  p<0.05 1

# 【 Cyclophosphamide 】 Low= 178.59  High= 177.065  Diff= -1.525  p<0.05  2
# 
# 【 Mitoxantrone 】 Low= 2.193  High= 1.902  Diff= -0.292  p= 0.1 ns   2

# 
# 【 Gemcitabine 】 Low= 0.855  High= 0.643  Diff= -0.211  p= 0.289 ns 
# 
# 【 Epirubicin 】 Low= 0.438  High= 0.4  Diff= -0.038  p= 0.279 ns 
# 

# 
# 【 Vincristine 】 Low= 0.351  High= 0.222  Diff= -0.129  p= 0.749 ns 
# 


