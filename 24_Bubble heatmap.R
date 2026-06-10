# 参考：https://mp.weixin.qq.com/s?__biz=MzU5NjE2ODU0OQ==&mid=2247483869&idx=1&sn=f35b5747cffd6ee98d7fcf0148dc4752&chksm=fe679bd7c91012c1846dfcb264ede176e6d9d689fb97afad07d1b71e0078dadb64e26acb4a6f&scene=21#wechat_redirect

# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./24_气泡热图//"
dir.create(out_dir, showWarnings = F, recursive = T)


# 加载包
library(ggplot2)
library(dplyr)
library(patchwork)


# 1. 读取数据（你的原始数据路径）
GSDC_top10 <- read.csv("../Outdata/18.Drug_predict/9_key_gene-drug_sensity_predicted_result.csv")

# 保留：至少在 1 个药物中 FDR ≤ 0.05 的基因
GSDC_top10_clean <- GSDC_top10 %>%
  group_by(symbol) %>%
  filter(any(fdr <= 0.05)) %>%  # 关键：只要有一个药物显著，就保留整个基因
  ungroup()

# 查看清理前后行数对比
cat("清理前行数：", nrow(GSDC_top10), "\n")
cat("清理后行数：", nrow(GSDC_top10_clean), "\n")

# 检查相关系数的范围
range(GSDC_top10_clean$cor, na.rm = TRUE)
# [1] -0.4843097  0.3388798

pdf("./24_气泡热图/1_drug_predicted_TOP10_daozhaun.pdf", width = 19, height = 6)

ggplot(GSDC_top10_clean, aes(x = symbol, y = drug)) +  # 交换 x 和 y
  geom_point(
    aes(
      size = -log10(fdr),
      fill = cor,
      color = fdr_label
    ),
    shape = 21,
    stroke = 1.0
  ) +
  scale_fill_gradient2(
    limits = c(-0.5, 0.5),
    low = "#3A3E96",
    mid = "white",
    high = "#AD3D3E",
    midpoint = 0,
    name = "Correlation"
  ) +
  scale_color_manual(
    values = c("FDR <= 0.05" = "black", "FDR > 0.05" = "gray95"),
    name = "FDR",
    guide = guide_legend(
      override.aes = list(fill = "white", size = 4)
    )
  ) +
  scale_size_continuous(
    range = c(4, 16),
    breaks = c(2.5, 5.0, 7.5, 10.0),
    name = "-Log10(FDR)"
  ) +
  labs(
    title = "Correlation between drug sensitivity and gene expression",
    x = "Gene Symbol",  # 交换 x 和 y 的标签
    y = "Drug"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 基因名旋转45度避免重叠
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid = element_line(color = "gray90"),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 10)
  )

dev.off()

write.csv(GSDC_top10, "./24_气泡热图//top10drup——prdiect_result.csv")
