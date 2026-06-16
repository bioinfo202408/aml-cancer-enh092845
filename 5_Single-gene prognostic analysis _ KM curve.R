# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_dir <- "./5_单基因预后分析_KM曲线/"
dir.create(out_dir, showWarnings = F, recursive = T)

# Load packages
library(survival)
library(survminer)
library(dplyr)
library(data.table)


# Import data
GSE1656 <- fread("/home/weili/Project/AML/snakefile/gene_data/homo/all_GSE165656_clean_expr.csv")
# write.csv(GSE1656, "../TCGA_data/GSE165656_expr_all_RNA.csv")
GSE1656_time <- read.csv("../TCGA_data/GSE165656_time_status.csv")
ID_SRR <- read.table("../TCGA_data/GSE165656_ID和SRR对应关系.txt")

# ===================== 1. Load your datasets =====================
# Expression matrix (rows: genes, columns: SRR samples)
expr <- GSE1656  
rownames(expr) <- expr$V1  # Column 1 stores gene names, set as row names
expr <- data.frame(expr[,-1], row.names = rownames(expr)) #rownames.force=TRUE to retain row names
head(rownames(expr))


# Survival data (ID = sample numeric ID, OS.time, OS)
time <- GSE1656_time  

# ID-SRR matching table
id_srr <- ID_SRR  
colnames(id_srr) <- c("SRR", "ID")

# ===================== 2. Input target gene for analysis =====================
gene <- "CLEC11A"  # Example: gene <- "FOXO1" Enh092845 IRF1

# ===================== 3. Extract expression levels of target gene =====================
if(!gene %in% rownames(expr)){
  stop(paste0(gene, " is not found in expression matrix!"))
}

gene_exp <- data.frame(
  SRR = colnames(expr),
  exp = as.numeric(expr[gene, ])
)

# ===================== 4. Match sample IDs and survival information =====================
gene_exp <- gene_exp %>%
  left_join(id_srr, by = "SRR") %>%
  left_join(time[, c("ID", "OS.time", "OS")], by = "ID") %>%
  filter(!is.na(OS.time), !is.na(OS))

# # ===================== 5.1 Stratify by median expression =====================
# gene_exp$risk_group <- ifelse(
#   gene_exp$exp > median(gene_exp$exp, na.rm=T), 
#   "High", "Low"
# )

# ===================== 5.2 Stratify by optimal cutoff value =====================
library(maxstat)
# Automatically calculate optimal cutoff
max_stat <- maxstat.test(Surv(OS.time, OS) ~ exp, data = gene_exp, smethod = "LogRank", pmethod = "HL")
cutoff <- max_stat$estimate

# Group samples
gene_exp$risk_group <- ifelse(gene_exp$exp > cutoff, "High", "Low")
cat(gene, "Optimal cutoff value =", cutoff, "\n")
# 
# # ===================== 5.3 Stratify by mean expression =====================
# gene_exp$risk_group <- ifelse(
#   gene_exp$exp > mean(gene_exp$exp, na.rm=T), 
#   "High", "Low"
# )
# 
# # ===================== 5.4 Stratify by quartiles =====================
# # Keep only top 25% and bottom 25% samples, remove middle 50% for better statistical difference
# q1 <- quantile(gene_exp$exp, 0.25, na.rm=T)
# q3 <- quantile(gene_exp$exp, 0.75, na.rm=T)
# 
# gene_exp <- gene_exp[gene_exp$exp <= q1 | gene_exp$exp >= q3, ]
# gene_exp$risk_group <- ifelse(gene_exp$exp >= q3, "High", "Low")


# ===================== 6. Survival analysis & KM curve plotting =====================
# 1) Fit survival curve
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = gene_exp)

# 2) Publication-level plot theme
pub_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=0.8),
        plot.title = element_text(hjust=0.5, size=14),
        axis.text = element_text(size=11, color="black"))

pdf("./5_单基因预后分析_KM曲线/CLEC11A_KM_curve_cutoff.pdf", width=6, height=5)
# 3) Draw KM curve
ggsurvplot(
  fit,
  data = gene_exp,
  title = paste0(gene, " - Overall Survival"),
  xlab = "Time (months)",
  ylab = "Survival probability",
  risk.table = TRUE,
  pval = TRUE,
  conf.int = F,
  palette = c("#2E86AB", "#E63946"),
  legend.labs = c(paste0(gene," Low"), paste0(gene," High")),
  ggtheme = pub_theme,
  surv.median.line = "hv"
)
dev.off()
