# Reference: https://mp.weixin.qq.com/s/pORXLeLI9O_WS06IZTIa-w   Half-violin + dot plot
# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialize] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output directory
out_dir <- "./6_Violin_HalfViolinPlot//"
dir.create(out_dir, showWarnings = F, recursive = T)

# ========================
# Load packages
# ========================
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)  # Package for half-violin plot

# ===============================================
### Key gene expression in Normal / Tumor samples: Half-violin dot plot #######
# ===============================================

# ========================
# Load raw data
# ========================
allRNA_expr <- fread("../Outdata/5.all_data_harmony/combat_allRNA_20260520.csv")
rownames(allRNA_expr) <- allRNA_expr$V1

group_info <- read.csv("../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv", row.names=1)
group_info$Group <- ifelse(group_info$Group == 1,"Tumor","Normal")

# ========================
# Extract expression matrix
# ========================
expr_mat <- as.data.frame(allRNA_expr)
rownames(expr_mat) <- expr_mat$V1
expr_mat <- expr_mat[, -1]  # Remove gene ID column, retain expression values

# ========================
# Select target genes
# ========================
# genes <- c("Enh092845", "CLEC11A", "IRF1")
# genes <- c("Enh092845")
genes <- c("CLEC11A")
# genes <- c("IRF1")
expr_sub <- expr_mat[genes, ]
expr_sub <- expr_sub[,-1]


# ========================
# Construct plotting data frame
# ========================
plot_df <- as.data.frame(t(expr_sub))
plot_df$sample <- rownames(plot_df)
plot_df$Group <- group_info[plot_df$sample, "Group"]

# Convert to long format
plot_df <- pivot_longer(
  plot_df,
  cols = all_of(genes),
  names_to = "Gene",
  values_to = "Expression"
)

# ========================
# Correct negative expression values + log2 transformation
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
# Filter extreme values: remove top/bottom 10% per group
# ======================
plot_df_clean <- plot_df %>%
  group_by(Group, Gene) %>%
  mutate(
    q10 = quantile(Expression, 0.3, na.rm = T),
    q90 = quantile(Expression, 0.7, na.rm = T)
  ) %>%
  filter(Expression >= q10 & Expression <= q90) %>%
  ungroup()


# ========================
# Group comparison setting
# ========================
my_comparisons <- list(c("Normal", "Tumor"))

# ========================
# Color palette
# ========================
group_colors <- c(
  "Tumor"  = "#AD3D3E",
  "Normal" = "#3A3E96"
)
# '#3A3E96', '#AD3D3E', '#50A293','#E8B75E'


# # ========================
# # Standard violin plot (backup)
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
# # Save figure
# # ========================
# ggsave(
#   filename = paste0(out_dir, genes, "_violin_Tumor_vs_Normal_FWSE_dataset_new.pdf"),
#   plot = p,
#   width = 4.8,
#   height = 4.8
# )


# ========================
# Nature-style half-violin dot plot
# ========================
p <- ggplot(plot_df_clean, aes(x = Group, y = Expression, color = Group)) +
  # 1. Right-sided half violin
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  # 2. Boxplot
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1
  ) +
  # 3. Mean point
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  # 4. Jitter dots
  geom_jitter(
    size = 1.8, alpha = 0.5, width = 0.15
  ) +
  # 5. Statistical significance
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0,
    label.y.npc = "top",
    size = 5
  ) +
  # 6. Color mapping
  scale_color_manual(values = group_colors) +
  # 7. Axis labels
  labs(
    x = NULL,
    y = paste0(genes, " Expression log2(TPM+1)")
  ) +
  # 8. Nature journal theme
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
# Save high-resolution PDF
# ========================
ggsave(
  filename = paste0(out_dir, genes, "_Nature_halfviolin_Tumor_vs_Normal.pdf"),
  plot = p,
  width = 4.5,
  height = 6
)








# ===============================================================================================================
# #### Risk score half-violin plot stratified by age, gender, ELN cytogenetics and risk group ####
# ===============================================================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix) 
library(ggpubr)
library(gghalves)

risk <- read.csv("../Outdata/7.Prognostic_analysis/8_train_lasso_risk.csv")
clina <- read.csv("../TCGA_data/TCGA_clina_ostime.csv")
clina2 <- read.csv("../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime_new_136_samples.csv")

risk$eln <- clina2$eln[match(risk$sample, clina2$ID)] 
risk$age <- clina$age[match(risk$sample, clina$sample)] 
risk$gender <- clina$gender[match(risk$sample, clina$sample)]

risk <- risk[,c("sample", "OS", "OS.time", "riskscore", "risk", "eln", "age", "gender")]
risk <- risk[risk$eln != "'--", ] 



# ====================== Data preprocessing
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

# ====================== Filter extreme values function
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

# ====================== Color palette
colors2 <- c("Low risk"="#3A3E96", "High risk"="#AD3D3E",
             "<65"="#3A3E96", ">=65"="#AD3D3E",
             "Female"="#3A3E96", "Male"="#AD3D3E")


# ======================================================================================
# Core function: Nature-style half-violin plot (compatible with original calls)
# ======================================================================================
plot_nc_halfviolin <- function(data, x, y, title, colors) {
  
  data_clean <- filter_extreme(data, x, y)
  
  ggplot(data_clean, aes(x = .data[[x]], y = .data[[y]], color = .data[[x]])) +
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth = 1
    ) +
    geom_boxplot(
      outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill = "white"
    ) +
    stat_summary(
      fun = "mean", geom = "point", shape = 20, size = 3,
      color = "black", fill = "black", alpha = 0.8
    ) +
    geom_jitter(
      size = 1.8, alpha = 0.5, width = 0.15
    ) +
    stat_compare_means(
      method = "wilcox.test", label = "p.signif",
      tip.length = 0, label.y.npc = "top", size = 5
    ) +
    scale_color_manual(values = colors) +
    labs(x = "", y = "Risk Score", title = title) +
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

# ====================== Generate 4 plots
# 1. Age stratification
p1 <- plot_nc_halfviolin(risk, "age_group_label", "riskscore_pos", "Risk by Age Group", colors2)
ggsave("./6_ViolinPlot/Risk_Age_Nature_halfviolin.pdf", p1, width=4.5, height=6)

# 2. Gender stratification
p2 <- plot_nc_halfviolin(risk, "gender_group_label", "riskscore_pos", "Risk by Gender Group", colors2)
ggsave("./6_ViolinPlot/Risk_Gender_Nature_halfviolin.pdf", p2, width=4.5, height=6)

# 3. Risk group comparison
p3 <- plot_nc_halfviolin(risk, "risk_group_label", "riskscore_pos", "Risk Group", colors2)
ggsave("./6_ViolinPlot/Risk_RiskGroup_Nature_halfviolin.pdf", p3, width=4.5, height=6)



# 4. ELN cytogenetic subgroups (multi-group half-violin)
# Reorder ELN levels
risk$eln_group_label <- factor(
  risk$eln,
  levels = c("Favorable", "Intermediate", "Adverse")
)

# ELN color palette
colors3 <- c("Favorable"="#50A293", "Intermediate"="#3A3E96", "Adverse"="#AD3D3E")

# Multi-group comparison pairs
eln_comparisons <- list(
  c("Favorable", "Intermediate"),
  c("Favorable", "Adverse"),
  c("Intermediate", "Adverse")
)

# Plot
p4 <- ggplot(risk, aes(x = eln_group_label, y = riskscore_pos, color = eln_group_label)) +
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill="white"
  ) +
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  geom_jitter(size = 1.8, alpha = 0.5, width = 0.15) +
  stat_compare_means(
    comparisons = eln_comparisons,
    method = "wilcox.test", label = "p.signif",
    tip.length = 0, label.y.npc = "top", size = 4.5
  ) +
  scale_color_manual(values = colors3) +
  labs(x = "", y = "Risk Score", title = "Risk by ELN Cytogenetic Group") +
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

# Save
ggsave("./6_ViolinPlot/Risk_ELN_Nature_halfviolin.pdf", p4, width=5.5, height=6)




# ===============================================================================================================
# #### Pseudotime comparison between monocyte high/low subgroups: Nature half-violin plot ####
# ===============================================================================================================

load("/home/weili/Project/AML/human/AML_combined_analyse/scRNA_analyse/20260407_analysis_GSE1116256/rawdata/Step10_Pseudotime_analysis_mono_highlow_subgroup_final.Rdata")
# Pseudotime boxplot with significance
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gghalves)

input.data = data.frame(group = HSMM$group, Pseudotime = HSMM$Pseudotime)

# ======================
# Filter extreme values
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
# Two-group comparison
# ======================
my_comparisons <- list(c("mono_cell_low", "mono_cell_high"))

# ======================
# Nature-style half-violin plot
# ======================
p <- ggplot(plot_df_clean, aes(x = group, y = Pseudotime, color = group)) +
  geom_half_violin(
    position = position_nudge(x = 0.15, y = 0),
    side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth=1
  ) +
  geom_boxplot(
    outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1, fill="white"
  ) +
  stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3,
    color = "black", fill = "black", alpha = 0.8
  ) +
  geom_jitter(size = 1.8, alpha = 0.5, width = 0.15) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", label = "p.signif",
    tip.length = 0, label.y.npc = "top", size = 5
  ) +
  scale_color_manual(values = c("mono_cell_high"="#AD3D3E", "mono_cell_low"="#3A3E96")) +
  labs(x = "", y = "Pseudotime") +
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
# Save high-res PDF
# ======================
ggsave(
  filename = "./6_ViolinPlot/Pseudotime_Nature_halfviolin.pdf",
  plot = p,
  width = 4.5,
  height = 6
)




# ===============================================================================================================
# #### Stromal Score, Immune Score, ESTIMATE Score and Tumor Purity: Nature half-violin plot ####
# ===============================================================================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(gghalves)

# Load data
core <- read.csv("../Outdata/12.Immune_analysis/8_estimate_scores_risk_group.csv")

# Preprocess
scores_df <- core
scores_df$risk <- factor(
  scores_df$risk, 
  levels = c(0, 1),
  labels = c("Low risk", "High risk")
)

# Extreme filter function
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

# Unified color palette
colors2 <- c("Low risk"="#3A3E96", "High risk"="#AD3D3E")
my_comparisons <- list(c("Low risk", "High risk"))

# ===================== Plot function for ESTIMATE metrics
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

# Generate four plots
p_stromal <- plot_estimate_nc("StromalScore", "Stromal Score", "Stromal Score")
p_immune  <- plot_estimate_nc("ImmuneScore", "Immune Score", "Immune Score")
p_estima  <- plot_estimate_nc("ESTIMATEScore", "ESTIMATE Score", "ESTIMATE Score")
p_purity  <- plot_estimate_nc("TumorPurity", "Tumor Purity", "Tumor Purity")


# Save
ggsave(
  "./6_ViolinPlot/ESTIMATE_stromal_halfviolin.pdf",
  plot = p_stromal,
  width = 4.5,
  height = 6
)

ggsave(
  "./6_ViolinPlot/ESTIMATE_immune_halfviolin.pdf",
  plot = p_immune,
  width = 4.5,
  height = 6
)

ggsave(
  "./6_ViolinPlot/ESTIMATE_estimate_halfviolin.pdf",
  plot = p_estima,
  width = 4.5,
  height = 6
)

ggsave(
  "./6_ViolinPlot/ESTIMATE_purity_halfviolin.pdf",
  plot = p_purity,
  width = 4.5,
  height = 6
)





# ===============================================================================================================
# #### CLEC11A expression across distinct cell subpopulations ####
# ===============================================================================================================

load("../scRNA_analyse/20260407_analysis_GSE1116256/rawdata/Step8.4.CLEC11A_expression_across_cell_subsets.Rdata")
# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggpubr)
library(gghalves)
library(RColorBrewer)

# Extract CLEC11A expression and cell type metadata
plot_df <- data.frame(
  celltype = seurat.data2$celltype_new,
  Expression = seurat.data2@assays$RNA@data["CLEC11A", ]
)

# Expression normalization: negative correction + log2 transform
plot_df <- plot_df %>%
  mutate(
    min_exp = min(Expression, na.rm = T),
    Expression = ifelse(min_exp < 0, Expression - min_exp, Expression),
    Expression = log2(Expression + 1)
  )

# Filter extreme values (retain central 80%)
plot_df_clean <- plot_df %>%
  group_by(celltype) %>%
  mutate(
    q10 = quantile(Expression, 0.4, na.rm=T),
    q90 = quantile(Expression, 0.6, na.rm=T)
  ) %>%
  filter(Expression >= q10 & Expression <= q90) %>%
  ungroup()

# Order cell types by mean expression (descending)
cell_order <- plot_df_clean %>%
  group_by(celltype) %>%
  summarise(mean_exp = mean(Expression, na.rm=T)) %>%
  arrange(desc(mean_exp)) %>%
  pull(celltype)

plot_df_clean$celltype <- factor(plot_df_clean$celltype, levels = cell_order)

# Color palette for cell populations
n <- length(unique(plot_df_clean$celltype))
cell_colors <- brewer.pal(n, "Set2")

# Nature-style half-violin plot
p <- ggplot(plot_df_clean, aes(x = celltype, y = Expression, color = celltype)) +
  geom_half_violin(
    position = position_nudge(x=0.15, y=0),
    side = 'r', trim=F, alpha=0.8, width=0.5, linewidth=1
  ) +
  geom_boxplot(
    outlier.shape = NA, width=0.35, alpha=0.8, linewidth=1
  ) +
  stat_summary(
    fun = "mean", geom = "point", shape=20, size=3,
    color="black", fill="black", alpha=0.8
  ) +
  geom_jitter(size=1.8, alpha=0.5, width=0.15) +
  scale_color_manual(values = cell_colors) +
  labs(
    x = "Cell Type",
    y = "CLEC11A Expression log2(TPM+1)",
    title = "CLEC11A Expression Across Cell Subsets"
  ) +
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

# Add significance comparison between Mono_high and Mono_low
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
  "./6_ViolinPlot/CLEC11A_cell_subset_halfviolin.pdf",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)



# ===============================================================================================================
# #### Drug sensitivity comparison between high and low risk groups ####
# ===============================================================================================================

# ========================
# Load packages
# ========================
library(ggplot2)
library(dplyr)

# ========================
# Data preparation
# ========================
df <- data.frame(read.csv("../Outdata/15.Drug_sensitivity/drug_risk_diff_top10.csv"))
# Color mapping for each drug
drug_colors <- c(
  "X5.Fluorouracil" = "#E63946",
  "BI.2536"        = "#F77F00",
  "BMS.754807"     = "#FCBF49",
  "Navitoclax"     = "#06D6A0",
  "ABT737"        = "#118AB2",
  "Dactolisib"     = "#073B4C",
  "Daporinad"      = "#9B5DE5",
  "PF.4708671"     = "#F15BB5",
  "Trametinib"     = "#00BBF9",
  "Cisplatin"      = "#8B5A2B"
)

# ========================
# Sort drugs by difference value
# ========================
df <- df %>%
  mutate(Drug = factor(Drug, levels = Drug[order(diff)]))

# ========================
# Lollipop plot for drug sensitivity difference
# ========================
p <- ggplot(df, aes(x = diff, y = Drug)) +
  geom_segment(
    aes(x = 0, xend = diff, y = Drug, yend = Drug, color = Drug),
    linewidth = 1.5, alpha = 0.7
  ) +
  geom_point(
    aes(color = Drug),
    size = 5, stroke = 1.5, fill = "white"
  ) +
  scale_color_manual(values = drug_colors) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(
    aes(label = sig, 
        hjust = ifelse(diff >= 0, -0.3, 1.3)),
    size = 4.5, fontface = "bold", color = "black"
  ) +
  labs(
    title = "Drug Sensitivity: High vs Low Risk Group",
    x = "Mean Difference (High risk - Low risk)",
    y = NULL
  ) +
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
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))

# ========================
# Save figure
# ========================
ggsave(
  filename = "./6_ViolinPlot/HighLowRisk_drug_sensitivity_lollipop.pdf",
  plot = p,
  width = 10,
  height = 5
)






# ========================
# Load packages
# ========================
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)

# ========================
# Load long-format drug sensitivity data
# ========================
drug_sens_target <- read.csv("../Outdata/15.Drug_sensitivity/7_AML_standard_drugs_sensitivity_highlow_risk.csv")

# Set Risk_Group factor order
drug_sens_target$Risk_Group <- factor(drug_sens_target$Risk_Group, levels = c("Low", "High"))

# ========================
# Group color palette
# ========================
group_colors <- c(
  "High" = "#AD3D3E",
  "Low"  = "#3A3E96"
)

# ========================
# Output directory
# ========================
out_dir <- "./6_ViolinPlot/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========================
# Extract all unique drug names
# ========================
drug_list <- unique(drug_sens_target$Drug)
cat("Total", length(drug_list), "drugs:", paste(drug_list, collapse = ", "), "\n")

# ========================
# Loop to generate separate half-violin plot per drug
# ========================
for(drug in drug_list) {
  
  df_sub <- drug_sens_target %>% filter(Drug == drug)
  
  mean_low <- mean(df_sub$Sensitivity[df_sub$Risk_Group == "Low"], na.rm = TRUE)
  mean_high <- mean(df_sub$Sensitivity[df_sub$Risk_Group == "High"], na.rm = TRUE)
  diff <- mean_high - mean_low
  
  test_result <- wilcox.test(Sensitivity ~ Risk_Group, data = df_sub)
  p_val <- test_result$p.value
  
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "ns")))
  
  cat("\n[", drug, "] Low=", round(mean_low, 3), 
      " High=", round(mean_high, 3), 
      " Diff=", round(diff, 3), 
      " p=", format(p_val, digits = 3), sig, "\n")
  
  # ========================
  # Plot half-violin dot plot
  # ========================
  p <- ggplot(df_sub, aes(x = Risk_Group, y = Sensitivity, color = Risk_Group)) +
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      side = 'r', trim = FALSE, alpha = 0.8, width = 0.5, linewidth = 1
    ) +
    geom_boxplot(
      outlier.shape = NA, width = 0.35, alpha = 0.8, linewidth = 1
    ) +
    stat_summary(
      fun = "mean", geom = "point", shape = 20, size = 3,
      color = "black", fill = "black", alpha = 0.8
    ) +
    geom_jitter(
      size = 2, alpha = 0.5, width = 0.15
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.03,
      label.y.npc = "top",
      size = 5
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      title = drug,
      subtitle = paste0("Low=", round(mean_low, 3), 
                        " | High=", round(mean_high, 3),
                        " | Diff=", round(diff, 3),
                        " | p", sig),
      x = NULL,
      y = "Drug Sensitivity (IC50)"
    ) +
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
  # Save PDF
  # ========================
  ggsave(
    filename = paste0(out_dir, drug, "_halfviolin_HighLowRisk.pdf"),
    plot = p,
    width = 4.5,
    height = 6
  )
  
  
}

# Log output reference
# [ Cytarabine ] Low= 7.434  High= 9.255  Diff= 1.821  p<0.05  *
# [ Venetoclax ] Low= 9.439  High= 9.862  Diff= 0.422  p<0.05 *
# [ Cyclophosphamide ] Low= 178.59  High= 177.065  Diff= -1.525  p<0.05  *
# [ Mitoxantrone ] Low= 2.193  High= 1.902  Diff= -0.292  p= 0.1 ns
# [ Gemcitabine ] Low= 0.855  High= 0.643  Diff= -0.211  p= 0.289 ns
# [ Epirubicin ] Low= 0.438  High= 0.4  Diff= -0.038  p= 0.279 ns
# [ Vincristine ] Low= 0.351  High= 0.222  Diff= -0.129  p= 0.749 ns
