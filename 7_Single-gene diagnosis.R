# Clear environment
rm(list = ls())
gc()

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【Initialization】Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output directory
out_dir <- "./7_Single_gene_diagnosis/"
dir.create(out_dir, showWarnings = F, recursive = T)


# Read expression matrix
Expr <- fread("../Outdata/5.all_data_harmony/combat_all_RNA_group_info.csv")
Expr <- Expr[,-1]

# Check gene existence
"IRF1" %in% Expr$V1
"Enh092845" %in% Expr$V1
"CLEC11A" %in% Expr$V1

# Single gene ROC diagnosis analysis
library(pROC)
library(ggplot2)
library(dplyr)
library(data.table)

# 1. Convert Expr to data.table (original format)
setDT(Expr)

# 2. Extract expression of three target genes and construct group label
expr_sub <- Expr[V1 %in% c("IRF1", "Enh092845", "CLEC11A")] %>%
  melt(id.vars = "V1", variable.name = "Sample", value.name = "Expression") %>%
  # Extract group label from sample name: Tumor/Normal
  mutate(Group = ifelse(grepl("_Tumor$", Sample), "Tumor", "Normal")) %>%
  # Convert to binary label (1=Tumor, 0=Normal)
  mutate(Label = ifelse(Group == "Tumor", 1, 0))

# 3. Define ROC plotting function (consistent with example style)
plot_roc <- function(gene_name, data) {
  # Subset data for target gene
  gene_data <- data[V1 == gene_name]
  
  # Calculate ROC curve and AUC value
  roc_obj <- roc(gene_data$Label, gene_data$Expression)
  auc_val <- round(auc(roc_obj), 3)
  ci_obj <- ci.auc(roc_obj)
  ci_lower <- round(ci_obj[1], 3)
  ci_upper <- round(ci_obj[3], 3)
  
  # Extract curve coordinate data
  roc_df <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities
  )
  
  # Output file name
  pdf_file <- paste0(out_dir, gene_name, "_ROC.pdf")
  
  # Save PDF file (key correction)
  pdf(pdf_file, width = 5, height = 5)
  # Plot figure (consistent with example style)
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#1f77b4", linewidth = 1.2) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = "#1f77b4", alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste0(gene_name, " Diagnostic ROC Curve"),
      x = "1-Specificity (FPR)",
      y = "Sensitivity (TPR)"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "gray90"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    annotate("text", x = 0.7, y = 0.3, 
             label = paste0("AUC: ", auc_val, "\n95% CI: ", ci_lower, "-", ci_upper),
             size = 5, hjust = 0)
  
  print(p)
  dev.off()
  
  # Return AUC and confidence interval information
  return(list(gene = gene_name, auc = auc_val, ci = c(ci_lower, ci_upper)))
}

# 4. Batch plot ROC curves for three genes
genes <- c("IRF1", "Enh092845", "CLEC11A")
results <- lapply(genes, function(g) plot_roc(g, expr_sub))

# 5. Output summary statistics
cat("=== Summary of diagnostic performance for three genes ===\n")
for (res in results) {
  cat(sprintf("%s: AUC=%.3f, 95%% CI=%.3f-%.3f\n", 
              res$gene, res$auc, res$ci[1], res$ci[2]))
}



# Violin plots for expression distribution of three genes

# Load required packages
library(ggplot2)
library(dplyr)
library(data.table)
library(ggsignif)

# 1. Convert data to data.table
setDT(Expr)

# 2. Extract expression data and grouping information of three genes
expr_sub <- Expr[V1 %in% c("IRF1", "Enh092845", "CLEC11A")] %>%
  melt(id.vars = "V1", variable.name = "Sample", value.name = "Expression") %>%
  # Extract group label from sample ID
  mutate(Group = ifelse(grepl("_Tumor$", Sample), "Tumor", "Normal")) %>%
  # Convert to factor to fix display order
  mutate(Group = factor(Group, levels = c("Tumor", "Normal")))

# 3. Define violin plot function (consistent with example style)
plot_violin <- function(gene_name, data, out_dir) {
  # Subset data of target gene
  gene_data <- data[V1 == gene_name]
  
  # Output file name
  pdf_file <- paste0(out_dir, gene_name, "_violin.pdf")
  
  # Generate plot
  p <- ggplot(gene_data, aes(x = Group, y = Expression, fill = Group)) +
    # Violin layer
    geom_violin(trim = FALSE, alpha = 0.7) +
    # Narrow boxplot in center
    geom_boxplot(width = 0.2, fill = "black", color = "black", outlier.shape = NA) +
    # Jitter scatter points
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
    # Significance annotation (Tumor vs Normal)
    geom_signif(comparisons = list(c("Tumor", "Normal")),
                map_signif_level = TRUE,
                y_position = max(gene_data$Expression, na.rm = TRUE) * 1.1,
                annotations = "***",
                tip_length = 0.01) +
    # Custom color scheme (consistent with blue/red style)
    scale_fill_manual(values = c("Tumor" = "#4169E1", "Normal" = "#DC143C")) +
    labs(
      x = "",
      y = paste0(gene_name, " expression level")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none"
    )
  
  # Save PDF file
  pdf(pdf_file, width = 4, height = 5)
  print(p)
  dev.off()
  
  cat("✅ Saved figure: ", pdf_file, "\n")
  
  # Return statistical test result
  wilcox_test <- wilcox.test(Expression ~ Group, data = gene_data)
  return(list(gene = gene_name, p.value = wilcox_test$p.value))
}

# 4. Batch generate violin plots
genes <- c("IRF1", "Enh092845", "CLEC11A")
results <- lapply(genes, function(g) plot_violin(g, expr_sub, out_dir))

# 5. Output statistical summary
cat("\n=== Wilcoxon rank-sum test results for differential expression ===\n")
for (res in results) {
  cat(sprintf("%s: p-value = %.2e\n", res$gene, res$p.value))
}

cat("\n🎉 All analysis finished! Figures saved in: ", out_dir, "\n")
