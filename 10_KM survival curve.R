# Reference: https://mp.weixin.qq.com/s/Hv7H1ocl4IgTb2Uf4AqgnQ

# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("【Initialization】Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output directory
out_dir <- "./10_KM_Survival_Curve/"
dir.create(out_dir, showWarnings = F, recursive = T)


# ===================== Top-journal standard KM survival curve Full Version =====================
library(survminer)
library(survival)  # Core survival analysis package, must load together
library(ggsurvfit)
library(patchwork)
library(ggpp)
library(dplyr)
library(scales)

# Load data
# risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv", header = T, row.names = 1)
risk <- read.csv("../Outdata/7.Prognostic analyse/9_test_risk.csv", header = T, row.names = 1)

# Factor setting for grouping
risk$risk_group <- factor(
  risk$risk, 
  levels = c(0, 1), 
  labels = c("Low risk", "High risk")
)

# Fit survival model & Cox regression
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = risk)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = risk)

# Extract HR and 95% confidence interval
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# Color palette
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# Plotting
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # Risk table
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # Censorship mark
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P value annotation
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR value annotation
  annotate("text",
           x = max(risk$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # Axis settings
  labs(
    # title = "Overall survival (training set)",  # Add title here
    title = "Overall survival (test set)",  # Add title here
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # Basic theme
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3,0.8,0.3,0.3,unit = "cm")
  )

# 【Key Fix】Combine into one page
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# Save figure
ggsave(
  # "./10_KM_Survival_Curve/train_KM_curve_full_top_style.pdf",
  "./10_KM_Survival_Curve/test_KM_curve_full_top_style.pdf",
  plot = p,
  width = 7, height = 5.5
)



# ===================== KM curve for age/gender subgroups =====================
### KM curve for age/gender/ELN subgroups
# ===================== Batch plot KM curves for training + test set
dir.create("./Outplot/7.Prognostic_Analysis/KM_plots", recursive = TRUE, showWarnings = FALSE)

library(survival)
library(ggsurvfit)
library(patchwork)
library(ggpp)
library(dplyr)
library(scales)

# ===================== 1. Load data 
train <- read.csv("../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime.csv", header = T, row.names = 1)
test  <- read.csv("../TCGA_data/GSE165656_LAML_879_feature_expr_with_ostime.csv", header = T, row.names = 1)

train_risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")
test_risk  <- read.csv("../Outdata/7.Prognostic analyse/9_test_risk.csv")
colnames(test_risk)[1] <- "sample"

# ===================== 2. Extract clinical information
train_clin <- train[, c("sample","age","gender","stage")]
colnames(train_clin) <- c("sample","Age","Gender","Stage")

test_clin  <- test[, c("sample","Age","Gender","stage")]
colnames(test_clin) <- c("sample","Age","Gender","Stage")

# ===================== 3. Merge risk score + clinical data 
train_risk2 <- train_risk %>% left_join(train_clin, by="sample")
test_risk2  <- test_risk  %>% left_join(test_clin,  by="sample")

# ===================== 4. Generate grouping labels
train_risk2$Age_group <- ifelse(train_risk2$Age < 65, "<65", "≥65")
test_risk2$Age_group  <- ifelse(test_risk2$Age < 65, "<65", "≥65")

train_risk2$Gender_group <- ifelse(train_risk2$Gender == "male", "Male", "Female")
test_risk2$Gender_group  <- ifelse(test_risk2$Gender == "M", "Male", "Female")


write.csv(train_risk2,"./10_KM_Survival_Curve/train_risk_age_gender.csv")
write.csv(test_risk2,"./10_KM_Survival_Curve/test_risk_age_gender.csv")

# train_risk2$Age_group <- ifelse(train_risk2$Age < 65, "<65", ">=65") # Age dichotomization
# test_risk2$Age_group <- ifelse(test_risk2$Age < 65, "<65", ">=65") # Age dichotomization
# train_risk2$Gneder_group <- ifelse(train_risk2$Gender == c("male"), "Male", "Female") # Gender dichotomization

# Filter by age group
# risk <- train_risk2[train_risk2$Age < 65,]
# risk <- train_risk2[train_risk2$Age >= 65,]
# Filter by gender group
# risk <- train_risk2[train_risk2$Gender %in% "female",]
risk <- train_risk2[train_risk2$Gender %in% "male",]
# Filter by ELN group
eln <- read.csv("../TCGA_data/TCGA-ELN.csv")
eln <- eln[,c("ID","ELN")]
risk <- train_risk2
risk$eln <- eln$ELN[match(risk$sample,eln$ID)] 
risk <- risk[ !risk$eln %in% c("", "--", "'--"), ]  # Remove samples without ELN label
risk <- risk[risk$eln %in% "Intermediate",] # Extract intermediate-risk patients
# Filter by TMB group
risk <- read.csv("./10_KM_Survival_Curve/TMB_KM.CSV")


{
# Factor setting for grouping
risk$risk_group <- factor(
  risk$risk, 
  levels = c(0, 1), 
  labels = c("Low risk", "High risk")
)

# Fit survival model & Cox regression
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = risk)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = risk)

# Extract HR and 95% confidence interval
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# Count events and total samples
tab_surv <- summary(fit)
event_low <- tab_surv$n.event[1]
event_high <- tab_surv$n.event[2]
n_low <- tab_surv$n[1]
n_high <- tab_surv$n[2]

# Color palette
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")
}
# Plotting
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # Risk table
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # Censorship mark
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P value annotation
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR value annotation
  annotate("text",
           x = max(risk$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # Axis settings
  labs(
    # title = "Overall survival (training set)",  # Add title here
    # title = "Overall survival (test set)",  # Add title here
    # title = "Age < 65",  # Add title here
    # title = "Age >= 65",  # Add title here
    # title = "Female",  # Add title here
    # title = "Male",  # Add title here
    title = "Intermediate Subgroup",  # Add title here
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # Basic theme
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3,0.8,0.3,0.3,unit = "cm")
  )

# title = "Overall survival (test set)",  # Add title here# 【Key Fix】Combine into one page
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# Save figure
ggsave(
  # "./10_KM_Survival_Curve/train_KM_curve_full_top_style.pdf",
  # "./10_KM_Survival_Curve/test_KM_curve_full_top_style.pdf",
  # "./10_KM_Survival_Curve/age<65_train_KM_curve_full_top_style.pdf",
  # "./10_KM_Survival_Curve/age>=65_train_KM_curve_full_top_style.pdf",
  # "./10_KM_Survival_Curve/female_train_KM_curve_full_top_style.pdf",
  # "./10_KM_Survival_Curve/male_train_KM_curve_full_top_style.pdf",
  "./10_KM_Survival_Curve/Intermediate_train_KM_curve_full_top_style.pdf",
  plot = p,
  width = 7, height = 5.5
)





# ===================== KM curve for four TMB subgroups =====================
# Four subgroups stratified by TMB and risk score
risk <- read.csv("./10_KM_Survival_Curve/TMB_KM.CSV")

fit <- survfit2(Surv(OS.time, OS) ~ group, data = risk)

# ===================== Dedicated 4-color palette for top journals =====================
col_pal <- c(
  "Low risk - Low TMB"    = "#009FC3",    # Blue
  "Low risk - High TMB"   = "#4DAF4A",    # Green
  "High risk - Low TMB"   = "#FF7F00",    # Orange
  "High risk - High TMB"  = "#B30437"     # Red
)

# ===================== Plotting (4 survival curves) =====================
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  
  # Risk table
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4
  ) +
  
  # Censorship mark
  add_censor_mark(shape = 1, size = 1.5, stroke = 1) +
  
  # Log-rank P value (auto-calculated for 4 groups)
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.15,
    hjust = 1, size = 4
  ) +
  
  # Title and axis labels
  labs(
    title = "Overall survival by Risk & TMB",
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  
  scale_x_continuous(breaks = c(0,25,50,75,100), expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  
  # Theme setting
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.5, 0.1),  # Place legend horizontally at bottom center for 4 groups
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.background = element_blank()
  )

# 【Key】Combine main curve and risk table into single page for 4 groups
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# Save output
ggsave(
  "./10_KM_Survival_Curve/TMB_Risk_4group_KM.pdf",
  plot = p,
  width = 7, height = 5.5  # Wider canvas recommended for 4-group plots
)





# ===================== Single-gene KM survival curve =====================

# Load packages (consistent with previous KM scripts)
library(survival)
library(ggsurvfit)
library(patchwork)
library(dplyr)
library(data.table)
library(maxstat)
library(scales)

# Import raw data
GSE1656 <- fread("/home/weili/Project/AML/snakefile/gene_data/homo/all_GSE165656_clean_expr.csv")
GSE1656_time <- read.csv("../TCGA_data/GSE165656_time_status.csv")
ID_SRR <- read.table("../TCGA_data/GSE165656_ID_and_SRR_mapping.txt")

# ===================== 1. Data preprocessing 
expr <- GSE1656  
rownames(expr) <- expr$V1
expr <- data.frame(expr[,-1], row.names = rownames(expr))
time <- GSE1656_time  
id_srr <- ID_SRR  
colnames(id_srr) <- c("SRR", "ID")

# ===================== 2. Target gene for analysis
gene <- "CLEC11A"
# gene <- "IRF1"
# gene <- "Enh092845"  

# ===================== 3. Extract gene expression & match survival data 
gene_exp <- data.frame(
  SRR = colnames(expr),
  exp = as.numeric(expr[gene, ])
)

gene_exp <- gene_exp %>%
  left_join(id_srr, by = "SRR") %>%
  left_join(time[, c("ID", "OS.time", "OS")], by = "ID") %>%
  filter(!is.na(OS.time), !is.na(OS))

# ===================== 4. Optimal cutoff grouping via maxstat test
max_stat <- maxstat.test(Surv(OS.time, OS) ~ exp, data = gene_exp, smethod = "LogRank", pmethod = "HL")
cutoff <- max_stat$estimate
gene_exp$risk_group <- ifelse(gene_exp$exp > cutoff, "Low", "High")

# ===================== 【Critical】Group factor setup identical to above scripts 
gene_exp$risk_group <- factor(
  gene_exp$risk_group, 
  levels = c("Low", "High"), 
  labels = c("Low risk", "High risk")
)

# ===================== 【Fully replicated】Model fitting, HR and event counting 
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = gene_exp)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = gene_exp)

hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

tab_surv <- summary(fit)
event_low <- tab_surv$n.event[1]
event_high <- tab_surv$n.event[2]
n_low <- tab_surv$n[1]
n_high <- tab_surv$n[2]

# Consistent color palette
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# ===================== 【Identical plotting pipeline】
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # Risk table
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # Censorship mark
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P value annotation
  add_pvalue(
    location = "annotation",
    x = max(gene_exp$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR value annotation
  annotate("text",
           x = max(gene_exp$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # Plot title
  labs(
    title = paste0(gene, " - Overall Survival"),
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # Theme fully consistent with prior code
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3,0.8,0.3,0.3,unit = "cm")
  )

# Combine curve and risk table into single PDF page (fix two-page output)
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# Export figure
ggsave(
  "./10_KM_Survival_Curve/SingleGene_CLEC11A_KM.pdf",
  # "./10_KM_Survival_Curve/SingleGene_IRF1_KM.pdf",
  # "./10_KM_Survival_Curve/SingleGene_Enh092845_KM.pdf",
  plot = p,
  width = 7, height = 5.5  
)
