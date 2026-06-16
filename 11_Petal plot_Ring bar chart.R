# Reference: https://mp.weixin.qq.com/s/udbFm4TTaVe9ZXjinizKNA
# Petal plot (radial bar chart) / Radial stacked bar chart
# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

# Create output directory
out_dir <- "./11_PetalPlot_RadialBarChart/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(patchwork)

# ========================
# Input data (3 validation cohorts)
# ========================
metrics <- c("Accuracy", "AUC", "F1", "Sensitivity", "Specificity")

Validation_1 <- c(0.7328, 0.7463, 0.7364, 0.7708, 0.7000)
Validation_2 <- c(0.7241, 0.7022, 0.7551, 0.7115, 0.7429)
Validation_3 <- c(0.6500, 0.7177, 0.6834, 0.8500, 0.4900)

df <- data.frame(
  Metric = metrics,
  Validation_1,
  Validation_2,
  Validation_3
)

# ========================
# Custom color palette (5 colors matched with metrics)
# ========================
custom_palette <- c(
  "Accuracy"      = "#b9d7e5",
  "AUC"           = "#8d96cc",
  "F1"            = "#fdcf9b",
  "Sensitivity"   = "#f89a7f",
  "Specificity"   = "#d4e3ae"
)

# ========================
# Petal plotting function (with legend)
# ========================
plot_single_petal <- function(data, group_name){
  
  petals <- 5
  petal_angle <- 360 / petals
  
  plot_data <- data %>%
    mutate(
      petal = row_number(),
      theta0 = petal * petal_angle
    ) %>%
    reframe(
      theta = theta0 + c(0, -petal_angle/2, 0, petal_angle/2, 0),
      r     = value * c(0, 0.6, 1, 0.6, 0),
      .by = c(Metric, value, petal, theta0)
    )
  
  label_data <- plot_data %>%
    group_by(Metric) %>%
    slice_max(r, n=1) %>%
    ungroup()
  
  ggplot(plot_data, aes(theta, r, group = petal, fill = Metric)) +
    ggforce::stat_bspline(geom = "area", n = 1000) +
    geom_text(data = label_data, 
              aes(label = sprintf("%.3f", value)), 
              size = 3.5, fontface="bold") +
    scale_fill_manual(values = custom_palette) +
    coord_radial() +
    labs(title = group_name, fill="Indicator") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "right",
      legend.text = element_text(size=11),
      legend.title = element_text(size=12, face="bold")
    )
}

# ========================
# Generate 3 separate petal plots
# ========================
p1 <- df %>% select(Metric, value=Validation_1) %>% plot_single_petal("Validation 1")
p2 <- df %>% select(Metric, value=Validation_2) %>% plot_single_petal("Validation 2")
p3 <- df %>% select(Metric, value=Validation_3) %>% plot_single_petal("Validation 3")

# Arrange plots in one row with shared legend
p_all <- p1 + p2 + p3 + plot_layout(ncol=3, guides = "collect") & 
  theme(legend.position = "right")

# ========================
# Export high-resolution figures
# ========================
ggsave(
  paste0(out_dir, "Model_Performance_Smooth_PetalPlot_3Groups_WithLegend.pdf"),
  plot = p_all,
  width = 14,
  height = 5)

ggsave(
  paste0(out_dir, "Model_Performance_Smooth_PetalPlot_1_WithLegend.pdf"),
  plot = p1,
  width = 5,
  height = 5)
ggsave(
  paste0(out_dir, "Model_Performance_Smooth_PetalPlot_2_WithLegend.pdf"),
  plot = p2,
  width = 5,
  height = 5)
ggsave(
  paste0(out_dir, "Model_Performance_Smooth_PetalPlot_3_WithLegend.pdf"),
  plot = p3,
  width = 5,
  height = 5)
