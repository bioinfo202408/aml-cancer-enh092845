# 参考：https://mp.weixin.qq.com/s/9bugJ-uy63PrpvIC6Ha0GA  


load("./20260407分析结果/rawdata/Step8.4.CLEC11A在不同细胞亚群表达量.Rdata")




# 清空环境并加载必要的包
rm(list = ls())
gc() #清理内存中不再使用的对象

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")


dir.create("./2_柱状图_蜂群_折线图//")




##########################不同类型RNA过滤前后的柱状图#################################
library(tidyverse)
library(scales)

# 数据准备（与你的原始数据一致）
df <- tibble(
  RNA_type = rep(c("mRNA","miRNA","lncRNA","eRNA"), each = 2),
  stage    = rep(c("Rawdata","filtered"), 4),
  count    = c(20046, 9274, 10473, 707, 192138, 4364, 524986, 25757)
) %>%
  mutate(
    stage = factor(stage, levels = c("Rawdata", "filtered")),
    RNA_type = factor(RNA_type, levels = c("mRNA","miRNA","lncRNA","eRNA"))
  )

# 定义你提供的配色，对应4种RNA类型
rna_colors <- c(
  "mRNA"    = "#839FBF",   # 蓝色
  "miRNA"   = "#93A89B",   # 灰绿
  "lncRNA"  = "#CAE0CA",   # 浅绿
  "eRNA"    = "#FFC89D"    # 浅橙
)

pdf("./2_柱状图_蜂群_折线图//1.RNA过滤前后对比图2.pdf",width = 9, height = 5)

ggplot(df, aes(x = RNA_type, y = count, fill = interaction(RNA_type, stage), alpha = stage)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_y_log10(labels = comma) +
  
  # 核心：用你提供的颜色给每种RNA上色，用透明度区分Rawdata/filtered
  scale_fill_manual(
    values = c(
      "mRNA.Rawdata" = "#839FBF",
      "mRNA.filtered" = "#839FBF",
      "miRNA.Rawdata" = "#93A89B",
      "miRNA.filtered" = "#93A89B",
      "lncRNA.Rawdata" = "#CAE0CA",
      "lncRNA.filtered" = "#CAE0CA",
      "eRNA.Rawdata" = "#FFC89D",
      "eRNA.filtered" = "#FFC89D"
    )
  ) +
  scale_alpha_manual(
    values = c("Rawdata" = 1, "filtered" = 0.6),  # 过滤后的柱子更淡一点
    guide = guide_legend(title = "Data Stage")
  ) +
  
  geom_text(
    aes(label = comma(count)),
    position = position_dodge(width = 0.9),
    vjust = -0.3, size = 3.5
  ) +
  
  labs(
    title = "RNA Type: Rawdata Counts vs Filtered Counts",
    x = "RNA Type", 
    y = "Number of RNA genes (log10)",
    fill = "RNA Type",
    alpha = "Data Stage"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  )

dev.off()



##########################hub ncRNA综合得分柱状图#################################
library(ggplot2)
library(dplyr)

# 1. 读入数据并整理
hub_score <- read.csv("../Outdata/25_AML关键基因查找_new/3_top10_ncRNA_hubsorce.csv")
hub_score <- hub_score[, c("node", "hub_score")]

# 2. 按hub_score从小到大排序，方便水平柱状图展示
hub_score_plot <- hub_score %>% 
  arrange(hub_score) %>%
  mutate(node = factor(node, levels = node))

# 3. 用你给的图配色（8个颜色循环给10个柱子）
color_palette <- c(
  "#F28C66",  # Erythro（灰）
  "#E5719A",
  "#E6C29A",  # Plasma cell（浅棕）
  "#974C01",
  "#F7D354",  # B cell（黄）
  "#89C451",  # T cell（绿）
  "#D673A3",  # HSC/MPP（粉紫）
  "#7B9CC2",  # Dendritic（蓝）
  "#999999",  # Mono_cell_low（橙）
  "#4CB8A1"   # Mono_cell_high（青）
)

# 循环复用颜色，适配10个柱子
color_palette <- rep(color_palette, length.out = nrow(hub_score_plot))

# 4. 绘图（水平柱状图，无网格线）
p <- ggplot(hub_score_plot, aes(x = node, y = hub_score)) +
  geom_col(aes(fill = node), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = color_palette) +
  coord_flip() +  # 水平柱状图
  labs(
    x = "Key ncRNA",
    y = "Hub score",
    title = "Top10 key ncRNA (Hub score)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),  # 去掉网格线
    axis.line = element_line(linewidth = 0.8, color = "black")
  )

# 5. 保存图片
pdf("./2_柱状图_蜂群_折线图/TOP10_ncRNA_barplot.pdf", width = 8, height = 5)
print(p)
dev.off()




##########################motif富集柱状图#################################
# 加载包
library(ggplot2)
library(dplyr)

# ===================== 1. 读取数据（你原来的代码不变） =====================
tf_enrich <- read.table(
  file = "../motif_analyse/3.hub_gene_motif_enrich/hub_gene_motif_enriched_TF_list_MI0022577.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)
colnames(tf_enrich) <- c("TF_ID", "TF_Name", "p_value")
tf_enrich <- tf_enrich[!duplicated(tf_enrich$TF_Name), ]

tf_enrich <- tf_enrich %>%
  mutate(
    p_value_corrected = ifelse(p_value == 0, 1e-300, p_value),
    log10_p = -log10(p_value_corrected)
  ) %>% arrange(desc(log10_p))

tf_enrich_Enh092845 <- tf_enrich

# AML 相关 TF 列表
AML_TF_list <- c("BD4", "LYL1", "RUNX1","RUNX2", "ELF1", "EVT6", "CEBPA", "MYC", "GATA2", "SPI1", "PU.1", "PU1", "HOXA9", 
                 "MEIS1", "RBMX", "L1", "CEBPA", "NPM1", "CBF", "TP53", "MDM2", "RUNX1T1", "GATA1",
                 "TAL1", "SCL", "WT1", "IKZF1", "IRF1", "ELF1", "TWIST1",
                 "MCL1", "BCL2")

tf_enrich_Enh092845_AML <- tf_enrich_Enh092845[tf_enrich_Enh092845$TF_Name %in% AML_TF_list, ]

# ===================== 2. 绘图：你要的水平柱状图 + 8色配色 + 无网格线 =====================
# 你图片里的8个标准配色
color_palette <- c(
  "#4CB8A1", "#E6C29A", "#F7D354", "#999999", 
  "#D673A3", "#7B9CC2", "#F28C66", "#89C451"
)

# ===================== 关键：从大到小排序 =====================
plot_data <- plot_data %>%
  arrange(desc(log10_p)) %>%  # 从大到小
  mutate(TF_Name = factor(TF_Name, levels = TF_Name))  # 固定排序


# 颜色循环适配TF数量
color_use <- rep(color_palette, length.out = nrow(plot_data))


# ===================== 3. 绘图：TF在横坐标 =====================
pdf("./2_柱状图_蜂群_折线图/TF_enrich_bar_xaxis.pdf", width=8, height=5)

ggplot(plot_data, aes(x = TF_Name, y = log10_p)) +
  geom_col(aes(fill = TF_Name), 
           color = "black", linewidth=0.5, width=0.7) +
  scale_fill_manual(values = color_use) +
  
  # 显著性阈值线
  geom_hline(yintercept = 1.301, 
             linetype="dashed", color="black", linewidth=0.6) +
  
  labs(
    x = "Transcription Factor (Motif)",
    y = "-log10(P-value)",
    title = "TF Motif Enrichment"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=14, face="bold"),
    axis.title = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1, size=10, color="black"),
    axis.text.y = element_text(size=10, color="black"),
    legend.position = "none",
    panel.grid = element_blank()
  )

dev.off()
