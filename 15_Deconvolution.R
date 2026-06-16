# ==============================================================
# Step 0: Environment Initialization & Path Configuration
# ==============================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Main working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("Initial working directory: ", getwd(), "\n")

# Load global environment configuration script
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create and switch to deconvolution output directory
out_folder <- "./15_Deconvolution/"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("Current output directory: ", getwd(), "\n")

# ==============================================================
# Step 1: Install and load all dependent packages
# ==============================================================
# if (!require(devtools)) install.packages("devtools")
# # Install EPIC deconvolution package
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
# Step 2: Download and preprocess TCGA-LAML Bulk RNA-seq data
# ==============================================================
## 2.1 Download TCGA-AML dataset
query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(query)
data_se <- GDCprepare(query)

## 2.2 Convert FPKM to TPM + log2(TPM+1) transformation
fpkm <- assay(data_se, "fpkm_unstrand")

fpkm2tpm <- function(fpkm_mat){
  total_per_sample <- colSums(fpkm_mat)
  tpm_mat <- t(t(fpkm_mat) / total_per_sample) * 1e6
  return(tpm_mat)
}

tpm <- fpkm2tpm(fpkm)
# Filter lowly expressed genes
keep <- rowSums(tpm > 0) > 0.1 * ncol(tpm)
tpm_filter <- tpm[keep, ]
bulk_expr <- as.matrix(tpm_filter)
bulk_expr_log <- log2(bulk_expr + 1)

## 2.3 Match samples with AMLFinder risk stratification information
risk_df <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")
bulk_sample_short <- substr(colnames(bulk_expr_log), 1, 12)
shared_sample <- intersect(bulk_sample_short, risk_df$sample)
idx <- which(bulk_sample_short %in% shared_sample)
bulk_expr_filter <- bulk_expr_log[, idx]

## 2.4 Convert Ensembl ID to Gene Symbol
ens_id <- sub("\\..*", "", rownames(bulk_expr_filter))
symbol <- mapIds(
  org.Hs.eg.db,
  keys = ens_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove genes without gene symbol annotation
bulk_symbol <- bulk_expr_filter[!is.na(symbol), ]
rownames(bulk_symbol) <- symbol[!is.na(symbol)]

# ==============================================================
# Step 3: Single-cell data processing: Marker genes & average expression reference profile for subpopulations
# ==============================================================
## 3.1 Load single-cell Seurat object
sc_obj <- qread("./GSE116256_Monocyte_HighLow_Subgroups_OtherSubpopulations.qs")
table(sc_obj$celltype)

## 3.3 Screen high-confidence marker genes
marker_raw <- read.csv("./AML_SingleCell_Subpopulations_MarkerGenes.csv")
marker_sel <- marker_raw %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  dplyr::select(cluster, gene) %>%
  distinct()
sig_gene_all <- unique(marker_sel$gene)

## 3.4 Intersect genes from Bulk data and marker gene list
gene_intersect <- intersect(rownames(bulk_symbol), sig_gene_all)
bulk_sig <- bulk_symbol[gene_intersect, ] %>% as.matrix()

## 3.5 Calculate average expression matrix for each single-cell subpopulation (EPIC reference signature)
sc_expr <- GetAssayData(sc_obj, assay = "RNA", slot = "data")
sc_meta <- data.frame(
  cell = colnames(sc_expr),
  celltype = sc_obj$celltype,
  stringsAsFactors = F
)

# Compute mean expression per cell type
sc_expr_t <- as.data.frame(t(sc_expr))
sc_expr_t$celltype <- sc_meta$celltype

sc_mean <- sc_expr_t %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  ungroup()

# Reshape to EPIC standard matrix: rows = genes, columns = cell subpopulations
sc_ref_mat <- sc_mean %>%
  pivot_longer(-celltype, names_to = "gene", values_to = "mean_expr") %>%
  pivot_wider(names_from = celltype, values_from = mean_expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

write.csv(as.data.frame(sc_ref_mat), "./AML_SingleCell_Subpopulations_MarkerGene_ExpressionMatrix.csv", row.names = T)

# ==============================================================
# Step 4: Gene alignment + Custom deconvolution with EPIC
# ==============================================================
# Triple gene matching: Bulk data / single-cell reference / marker genes
final_gene <- intersect(rownames(bulk_sig), rownames(sc_ref_mat))
bulk_epic_in <- bulk_sig[final_gene, ]
sc_ref_epic_in <- sc_ref_mat[final_gene, ]

# Assemble custom reference list for EPIC
epic_ref <- list(
  refProfiles = sc_ref_epic_in,
  sigGenes = final_gene
)

# Run deconvolution algorithm
epic_res <- EPIC(bulk = bulk_epic_in, reference = epic_ref)

# Extract cell fraction output
cell_frac <- epic_res$cellFractions %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")
write.csv(cell_frac, "./Deconvolution_Results.csv", row.names = F)

# Check convergence status of iterative fitting
table(epic_res$fit.gof$convergeCode)

# ==============================================================
# Step 5: Merge risk grouping metadata + statistical testing + boxplot visualization
# ==============================================================
# Match sample IDs
cell_frac$sample_short <- substr(cell_frac$sample_id, 1, 12)
merge_df <- left_join(
  cell_frac,
  risk_df[, c("sample", "risk")],
  by = c("sample_short" = "sample")
)

# Convert numeric risk value to categorical factor
merge_df$risk <- factor(merge_df$risk, levels = c(0,1), labels = c("Low","High"))

# Target cell subpopulation: mono_cell_high
target <- "mono_cell_low"

# Wilcoxon rank-sum test
form <- as.formula(paste0("`", target, "` ~ risk"))
wt_test <- wilcox.test(form, data = merge_df)
cat("\n===== ", target, " Inter-group comparison test results =====\n")
print(wt_test)

# Wilcoxon rank sum test with continuity correction
# 
# data:  mono_cell_low by risk
# W = 836, p-value = 8.732e-06
# alternative hypothesis: true location shift is not equal to 0

# # Draw boxplot and export PDF (no GUI on server, save file only)
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
# # Avoid special characters in filename
# pdf_name <- paste0(gsub("/", "_", target), "_fraction_boxplot_mono_cell_low.pdf")
# ggsave(pdf_name, plot = p, width = 6, height = 5, dpi = 300, device = "pdf")
# 

# Visualization workflow
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gghalves)
library(rlang)

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/15_Deconvolution/")
out_dir <- "./"

# Import data
cell_frac <- fread("./Deconvolution_Results_Latest.csv", data.table = F)
risk <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv", row.names = NULL)

# Match sample IDs
cell_frac$sample_short <- substr(cell_frac$sample_id, 1, 12)
merge_df <- dplyr::left_join(
  cell_frac,
  risk[, c("sample", "risk")],
  by = c("sample_short" = "sample")
)

# Define patient risk groups
merge_df$Group <- ifelse(merge_df$risk == 0, "Low Risk", "High Risk")

# Target cell population
# target_gene <- "mono_cell_low"
target_gene <- "mono_cell_high"

# Critical fix: add dplyr:: prefix to resolve namespace conflict, !!sym for flexible column name
plot_df <- merge_df %>% 
  dplyr::select(sample_short, Group, !!sym(target_gene)) %>%
  dplyr::rename(Expression = !!sym(target_gene))

# Filter extreme outliers (retain middle 80% data)
# plot_df_clean <- plot_df %>%
#   dplyr::group_by(Group) %>%
#   dplyr::mutate(
#     q10 = quantile(Expression, 0.2, na.rm = TRUE),
#     q90 = quantile(Expression, 0.8, na.rm = TRUE)
#   ) %>%
#   dplyr::filter(Expression >= q10 & Expression <= q90) %>%
#   dplyr::ungroup()
plot_df_clean <- plot_df 
# Adjust x-axis order to place High Risk group on the right
plot_df_clean$Group <- factor(plot_df_clean$Group,
                              levels = c("Low Risk", "High Risk"))

# Plot comparison settings
my_comparisons <- list(c("Low Risk", "High Risk"))
group_colors <- c("High Risk" = "#AD3D3E", "Low Risk" = "#3A3E96")

# Half-violin rain cloud plot
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

# Export figure
fname <- paste0(gsub("[ /]", "_", target_gene), "_RiskGroup_NC_RainCloudPlot.pdf")
ggsave(fname, plot = p, width = 4.5, height = 6)

# 1. Calculate median value of mono_cell_high
median_mono <- median(KM$mono_cell_high, na.rm = TRUE)
cat("Median value of mono_cell_high: ", median_mono, "\n")

# 2. Stratify patients by median cutoff
KM$mono_group <- ifelse(KM$mono_cell_high >= median_mono, "mono_high", "mono_low")

# 3. Check grouping distribution
table(KM$mono_group)
head(KM[, c("sample_id", "mono_cell_high", "mono_group")])

# Load survival analysis packages (consistent with original script)
library(survival)
library(ggsurvfit)
library(patchwork)
library(dplyr)
library(data.table)
library(maxstat)
library(scales)
library(survminer)

# ===================== 1. Data merging: Deconvolution cell fraction results + survival metadata =====================
# Match survival time, event status and sample ID
# Stratify patients by median value of mono_cell_high, compare Kaplan-Meier survival curves between subgroups
KM <- read.csv("./Deconvolution_Results_Latest.csv")
risk_df <- read.csv("../../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")

KM$sample_short_id <- substr(KM$sample_id, 1, 12)
# Merge data by truncated sample ID
merge_surv <- KM %>%
  dplyr::left_join(risk_df[, c("sample", "OS", "OS.time")], 
                   by = c("sample_short_id" = "sample")) %>%
  dplyr::filter(!is.na(OS.time), !is.na(OS))

# ===================== 2. Patient grouping based on mono_cell_high median cutoff =====================
# Calculate median threshold
median_cut <- median(merge_surv$mono_cell_high, na.rm = TRUE)
cat("Median cutoff value for mono_cell_high: ", median_cut, "\n")

# Group definition: above median = High risk, below median = Low risk (consistent with original logic)
merge_surv$risk_group <- ifelse(merge_surv$mono_cell_high > median_cut, "High", "Low")

# Factor level order & labels (retain original style)
merge_surv$risk_group <- factor(
  merge_surv$risk_group,
  levels = c("Low", "High"),
  labels = c("Low risk", "High risk")
)

# ===================== 3. Fit survival model, extract P value, HR and confidence interval =====================
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = merge_surv)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = merge_surv)

# Extract hazard ratio, 95% CI and Log-rank P value
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low  <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# Extract sample count and event count for each group (for reference)
tab_surv <- summary(fit)
n_low  <- tab_surv$n[1]
n_high <- tab_surv$n[2]

# ===================== 4. Color palette (fully consistent with original scheme) =====================
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# ===================== 5. Generate Kaplan-Meier plot (reproduce original figure style) =====================
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # Risk table below curve
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # Censorship markers
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P value annotation
  add_pvalue(
    location = "annotation",
    x = max(merge_surv$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR and 95% CI text label
  annotate("text",
           x = max(merge_surv$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ", hr_val, " (95%CI: ", hr_low, "-", hr_high, ")"),
           hjust = 1, size = 4.2) +
  # Axis labels
  labs(
    title = "mono_cell_high - Overall Survival",
    x = "Time (months)",
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03, 0)) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # Global theme settings
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

# Vertical layout: survival curve + risk table
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# ===================== 6. Export KM curve PDF =====================
ggsave(
  ".//mono_cell_high_Median_KM.pdf",
  plot = p,
  width = 7, height = 5.5
)

# Print sample count per risk group
table(merge_surv$risk_group)
