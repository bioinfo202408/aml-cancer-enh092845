# Reference: https://mp.weixin.qq.com/s?__biz=MzU5NjE2ODU0OQ==&mid=2247483869&idx=1&sn=f35b5747cffd6ee98d7fcf0148dc4752&chksm=fe679bd7c91012c1846dfcb264ede176e6d9d689fb97afad07d1b71e0078dadb64e26acb4a6f&scene=21#wechat_redirect

# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plotting_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output directory
out_dir <- "./24_Bubble_Heatmap//"
dir.create(out_dir, showWarnings = F, recursive = T)


# Load packages
library(ggplot2)
library(dplyr)
library(patchwork)


# 1. Import data (raw data path)
GSDC_top10 <- read.csv("../Outdata/18.Drug_predict/9_key_gene-drug_sensity_predicted_result.csv")

# Retain genes with FDR ≤ 0.05 in at least one drug
GSDC_top10_clean <- GSDC_top10 %>%
  group_by(symbol) %>%
  filter(any(fdr <= 0.05)) %>%  # Key: Keep the whole gene if significant in any single drug
  ungroup()

# Compare row counts before and after filtering
cat("Row count before filtering: ", nrow(GSDC_top10), "\n")
cat("Row count after filtering: ", nrow(GSDC_top10_clean), "\n")

# Check correlation coefficient range
range(GSDC_top10_clean$cor, na.rm = TRUE)
# [1] -0.4843097  0.3388798

pdf("./24_Bubble_Heatmap/1_drug_predicted_TOP10_reversed_axis.pdf", width = 19, height = 6)

ggplot(GSDC_top10_clean, aes(x = symbol, y = drug)) +  # Swap x and y axis
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
    x = "Gene Symbol",  # Axis labels after swapping x and y
    y = "Drug"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate gene names by 45° to avoid overlap
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid = element_line(color = "gray90"),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 10)
  )

dev.off()

write.csv(GSDC_top10, "./24_Bubble_Heatmap//top10_drug_predicted_result.csv")
