# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

# Create output folder
out_dir <- "./12_ConfusionMatrix/"
dir.create(out_dir, showWarnings = F, recursive = T)


# Clear environment again
rm(list = ls())
gc()

# Working directory (modify according to your path)
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
out_dir <- "./12_ConfusionMatrix/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Load packages
library(ggplot2)
library(patchwork)

# ========================
# 1. Data preparation
# ========================
# Left plot confusion matrix data
cm1 <- matrix(c(
  35, 11,  # Predict 0: Reference 0 / Reference 1
  15, 38   # Predict 1: Reference 0 / Reference 1
), nrow = 2, byrow = TRUE)

# Right plot confusion matrix data
cm2 <- matrix(c(
  26, 15,
  9, 37
), nrow = 2, byrow = TRUE)

# Convert matrix to long-format data frame
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

# Add labels with count and proportion for target cells
df1$label <- paste0(df1$count)
df1$label[df1$Prediction == 1 & df1$Reference == 1] <- "38\n(0.77)"
df1$label[df1$Prediction == 0 & df1$Reference == 0] <- "35\n(0.70)"

df2$label <- paste0(df2$count)
df2$label[df2$Prediction == 1 & df2$Reference == 1] <- "37\n(0.71)"
df2$label[df2$Prediction == 0 & df2$Reference == 0] <- "26\n(0.74)"

# Custom fill colors
colors <- c("#b9d7e5", "#8d96c0")
# ========================
# 2. Plot function to replicate figure style
# ========================
plot_cm <- function(data, colors) {
  ggplot(data, aes(x = Reference, y = Prediction, fill = factor(as.numeric(Reference) == as.numeric(Prediction)))) +
    geom_tile(color = "black", linewidth = 1.2) +  # Black tile border
    geom_text(aes(label = label), color = "white", size = 12, fontface = "bold") +
    scale_fill_manual(values = colors) +
    scale_y_discrete(limits = rev(levels(data$Prediction))) +  # Place 1 on top to match original figure
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
# 3. Generate two plots
# ========================
p1 <- plot_cm(df1, colors)
p2 <- plot_cm(df2, colors)


# ========================
# 4. Export high-resolution PDF
# ========================
ggsave(
  paste0(out_dir, "External_Cohort_1.pdf"),
  plot = p1,
  width = 5,
  height = 5
)
ggsave(
  paste0(out_dir, "External_Cohort_2.pdf"),
  plot = p2,
  width = 5,
  height = 5
)

cat("✅ Two confusion matrix figures generated successfully!\n")
