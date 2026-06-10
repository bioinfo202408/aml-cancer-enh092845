# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

# 创建输出目录
out_dir <- "./12_混淆矩阵/"
dir.create(out_dir, showWarnings = F, recursive = T)


# 清空环境
rm(list = ls())
gc()

# 工作路径（可按你的修改）
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
out_dir <- "./12_混淆矩阵/"
dir.create(out_dir, showWarnings = F, recursive = T)

# 加载包
library(ggplot2)
library(patchwork)

# ========================
# 1. 数据准备
# ========================
# 左图数据
cm1 <- matrix(c(
  35, 11,  # 预测0: 参考0/1
  15, 38   # 预测1: 参考0/1
), nrow = 2, byrow = TRUE)

# 右图数据
cm2 <- matrix(c(
  26, 15,
  9, 37
), nrow = 2, byrow = TRUE)

# 转为长格式数据框
cm_to_df <- function(cm) {
  df <- expand.grid(
    Prediction = factor(c(0,1), levels = c(0,1)),
    Reference  = factor(c(0,1), levels = c(0,1))
  )
  df$count <- as.vector(t(cm))
  df
}

df1 <- cm_to_df(cm1)
df2 <- cm_to_df(cm2)

# 给需要标注比例的格子加上数值（按你图中的标注）
df1$label <- paste0(df1$count)
df1$label[df1$Prediction == 1 & df1$Reference == 1] <- "38\n(0.77)"
df1$label[df1$Prediction == 0 & df1$Reference == 0] <- "35\n(0.70)"

df2$label <- paste0(df2$count)
df2$label[df2$Prediction == 1 & df2$Reference == 1] <- "37\n(0.71)"
df2$label[df2$Prediction == 0 & df2$Reference == 0] <- "26\n(0.74)"

# 配色（你指定的两个颜色）
colors <- c("#b9d7e5", "#8d96c0")
# ========================
# 2. 绘图函数（复刻你图中的样式）
# ========================
plot_cm <- function(data, colors) {
  ggplot(data, aes(x = Reference, y = Prediction, fill = factor(as.numeric(Reference) == as.numeric(Prediction)))) +
    geom_tile(color = "black", linewidth = 1.2) +  # 黑色边框
    geom_text(aes(label = label), color = "white", size = 12, fontface = "bold") +
    scale_fill_manual(values = colors) +
    scale_y_discrete(limits = rev(levels(data$Prediction))) +  # 让1在上方，和原图一致
    labs(x = "Reference", y = "Prediction") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank()
    )
}

# ========================
# 3. 绘制两张图并拼接
# ========================
p1 <- plot_cm(df1, colors)
p2 <- plot_cm(df2, colors)


# ========================
# 4. 保存高清PDF
# ========================
ggsave(
  paste0(out_dir, "外部训练集1.pdf"),
  plot = p1,
  width = 5,
  height = 5
  )
ggsave(
  paste0(out_dir, "外部训练集2.pdf"),
  plot = p2,
  width = 5,
  height = 5
)

cat("✅ 两张混淆矩阵已绘制完成！\n")