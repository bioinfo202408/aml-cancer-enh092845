# 参考：https://mp.weixin.qq.com/s/Hv7H1ocl4IgTb2Uf4AqgnQ

# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./10_KM生存曲线/"
dir.create(out_dir, showWarnings = F, recursive = T)


# ===================== 顶刊标准KM生存曲线 完整版 =====================
library(survminer)
library(survival)  # 生存分析核心包，必须一起加载
library(ggsurvfit)
library(patchwork)
library(ggpp)
library(dplyr)
library(scales)

# 读取数据
# risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv", header = T, row.names = 1)
risk <- read.csv("../Outdata/7.Prognostic analyse/9_test_risk.csv", header = T, row.names = 1)

# 分组因子设置
risk$risk_group <- factor(
  risk$risk, 
  levels = c(0, 1), 
  labels = c("Low risk", "High risk")
)

# 拟合生存模型与Cox回归
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = risk)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = risk)

# 提取HR与95%置信区间
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# 配色
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# 绘图
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # 风险表
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # 删失标记
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P值标注
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR值标注
  annotate("text",
           x = max(risk$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # 坐标轴设置
  labs(
    # title = "Overall survival (training set)",  # 这里加标题
    title = "Overall survival (test set)",  # 这里加标题
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # 基础主题
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 居中、加粗
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3,0.8,0.3,0.3,unit = "cm")
  )

# 【关键修复】合并成一页
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# 保存图片
ggsave(
  # "./10_KM生存曲线/train_KM_curve_full_top_style.pdf",
  "./10_KM生存曲线/test_KM_curve_full_top_style.pdf",
  plot = p,
  width = 7, height = 5.5
)



# ===================== age/gender亚群 KM曲线 =====================
### age/gender/ELN亚群 KM曲线
# ===================== 批量绘制训练集 + 测试集 KM 曲线
dir.create("./Outplot/7.Prognostic_analyse/KM_plots", recursive = TRUE, showWarnings = FALSE)

library(survival)
library(ggsurvfit)
library(patchwork)
library(ggpp)
library(dplyr)
library(scales)

# ===================== 1. 读取数据 
train <- read.csv("../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime.csv", header = T, row.names = 1)
test  <- read.csv("../TCGA_data/GSE165656_LAML_879_feature_expr_with_ostime.csv", header = T, row.names = 1)

train_risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv")
test_risk  <- read.csv("../Outdata/7.Prognostic analyse/9_test_risk.csv")
colnames(test_risk)[1] <- "sample"

# ===================== 2. 提取临床信息
train_clin <- train[, c("sample","age","gender","stage")]
colnames(train_clin) <- c("sample","Age","Gender","Stage")

test_clin  <- test[, c("sample","Age","Gender","stage")]
colnames(test_clin) <- c("sample","Age","Gender","Stage")

# ===================== 3. 合并风险值 + 临床 
train_risk2 <- train_risk %>% left_join(train_clin, by="sample")
test_risk2  <- test_risk  %>% left_join(test_clin,  by="sample")

# ===================== 4. 生成分组
train_risk2$Age_group <- ifelse(train_risk2$Age < 65, "<65", "≥65")
test_risk2$Age_group  <- ifelse(test_risk2$Age < 65, "<65", "≥65")

train_risk2$Gender_group <- ifelse(train_risk2$Gender == "male", "Male", "Female")
test_risk2$Gender_group  <- ifelse(test_risk2$Gender == "M", "Male", "Female")


write.csv(train_risk2,"./10_KM生存曲线/train_risk_age_gender.csv")
write.csv(test_risk2,"./10_KM生存曲线/test_risk_age_gender.csv")

# train_risk2$Age_group <- ifelse(train_risk2$Age < 65, "<65", ">=65") #年龄二分类
# test_risk2$Age_group <- ifelse(test_risk2$Age < 65, "<65", ">=65") #年龄二分类
# train_risk2$Gneder_group <- ifelse(train_risk2$Gender == c("male"), "Male", "Female") #性别二分类

# 按年龄分组
# risk <- train_risk2[train_risk2$Age < 65,]
# risk <- train_risk2[train_risk2$Age >= 65,]
# 按性别分组
# risk <- train_risk2[train_risk2$Gender %in% "female",]
risk <- train_risk2[train_risk2$Gender %in% "male",]
# 按ELN分组
eln <- read.csv("../TCGA_data/TCGA-ELN.csv")
eln <- eln[,c("ID","ELN")]
risk <- train_risk2
risk$eln <- eln$ELN[match(risk$sample,eln$ID)] 
risk <- risk[ !risk$eln %in% c("", "--", "'--"), ]  #去除没有eln分组的
risk <- risk[risk$eln %in% "Intermediate",] #提取中等风险样本
# 按TMB分组
risk <- read.csv("./10_KM生存曲线/TMB_KM.CSV")


{
# 分组因子设置
risk$risk_group <- factor(
  risk$risk, 
  levels = c(0, 1), 
  labels = c("Low risk", "High risk")
)

# 拟合生存模型与Cox回归
fit <- survfit2(Surv(OS.time, OS) ~ risk_group, data = risk)
cox_fit <- coxph(Surv(OS.time, OS) ~ risk_group, data = risk)

# 提取HR与95%置信区间
hr_val <- round(summary(cox_fit)$conf.int[1], 2)
hr_low <- round(summary(cox_fit)$conf.int[3], 2)
hr_high <- round(summary(cox_fit)$conf.int[4], 2)
logrank_p <- surv_pvalue(fit)$pval

# 统计事件数与总样本
tab_surv <- summary(fit)
event_low <- tab_surv$n.event[1]
event_high <- tab_surv$n.event[2]
n_low <- tab_surv$n[1]
n_high <- tab_surv$n[2]

# 配色
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")
}
# 绘图
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # 风险表
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # 删失标记
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P值标注
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR值标注
  annotate("text",
           x = max(risk$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # 坐标轴设置
  labs(
    # title = "Overall survival (training set)",  # 这里加标题
    # title = "Overall survival (test set)",  # 这里加标题
    # title = "Age < 65",  # 这里加标题
    # title = "Age >= 65",  # 这里加标题
    # title = "Female",  # 这里加标题
    # title = "Male",  # 这里加标题
    title = "Intermediate Subgroup",  # 这里加标题
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # 基础主题
  theme_classic() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 居中、加粗
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.1),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = margin(0.3,0.8,0.3,0.3,unit = "cm")
  )

# title = "Overall survival (test set)",  # 这里加标题# 【关键修复】合并成一页
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# 保存图片
ggsave(
  # "./10_KM生存曲线/train_KM_curve_full_top_style.pdf",
  # "./10_KM生存曲线/test_KM_curve_full_top_style.pdf",
  # "./10_KM生存曲线/age<65_train_KM_curve_full_top_style.pdf",
  # "./10_KM生存曲线/age>=65_train_KM_curve_full_top_style.pdf",
  # "./10_KM生存曲线/female_train_KM_curve_full_top_style.pdf",
  # "./10_KM生存曲线/male_train_KM_curve_full_top_style.pdf",
  "./10_KM生存曲线/Intermediate_train_KM_curve_full_top_style.pdf",
  plot = p,
  width = 7, height = 5.5
)





# ===================== 按TMB分组  四个分组 KM曲
# 按TMB分组  四个分组
risk <- read.csv("./10_KM生存曲线/TMB_KM.CSV")

fit <- survfit2(Surv(OS.time, OS) ~ group, data = risk)

# ===================== 4 组专用配色（顶刊 4 色）
col_pal <- c(
  "Low risk - Low TMB"    = "#009FC3",    # 蓝色
  "Low risk - High TMB"   = "#4DAF4A",    # 绿色
  "High risk - Low TMB"   = "#FF7F00",    # 橙色
  "High risk - High TMB"  = "#B30437"     # 红色
)

# ===================== 绘图（4 条曲线）
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  
  # 风险表
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4
  ) +
  
  # 删失点
  add_censor_mark(shape = 1, size = 1.5, stroke = 1) +
  
  # Log-rank P 值（4 组也能自动算）
  add_pvalue(
    location = "annotation",
    x = max(risk$OS.time)*0.95, y = 0.15,
    hjust = 1, size = 4
  ) +
  
  # 标题与坐标轴
  labs(
    title = "Overall survival by Risk & TMB",
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  
  scale_x_continuous(breaks = c(0,25,50,75,100), expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  
  # 主题
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.5, 0.1),  # 4组图例放底部中间更美观
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.background = element_blank()
  )

# 【关键】4 组合并成一页
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# 保存
ggsave(
  "./10_KM生存曲线/TMB_Risk_4group_KM.pdf",
  plot = p,
  width = 7, height = 5.5  # 4组图适当加宽更好看
)





# ===================== 单基因 KM曲线 =====================

# 加载包（和你上面的KM代码完全一致）
library(survival)
library(ggsurvfit)
library(patchwork)
library(dplyr)
library(data.table)
library(maxstat)
library(scales)

# 导入数据
GSE1656 <- fread("/home/weili/Project/AML/snakefile/gene_data/homo/all_GSE165656_clean_expr.csv")
GSE1656_time <- read.csv("../TCGA_data/GSE165656_time_status.csv")
ID_SRR <- read.table("../TCGA_data/GSE165656_ID和SRR对应关系.txt")

# ===================== 1. 数据整理 
expr <- GSE1656  
rownames(expr) <- expr$V1
expr <- data.frame(expr[,-1], row.names = rownames(expr))
time <- GSE1656_time  
id_srr <- ID_SRR  
colnames(id_srr) <- c("SRR", "ID")

# ===================== 2. 分析基因
gene <- "CLEC11A"
# gene <- "IRF1"
# gene <- "Enh092845"  

# ===================== 3. 提取表达 + 匹配生存 
gene_exp <- data.frame(
  SRR = colnames(expr),
  exp = as.numeric(expr[gene, ])
)

gene_exp <- gene_exp %>%
  left_join(id_srr, by = "SRR") %>%
  left_join(time[, c("ID", "OS.time", "OS")], by = "ID") %>%
  filter(!is.na(OS.time), !is.na(OS))

# ===================== 4. 最优截断值分组
max_stat <- maxstat.test(Surv(OS.time, OS) ~ exp, data = gene_exp, smethod = "LogRank", pmethod = "HL")
cutoff <- max_stat$estimate
gene_exp$risk_group <- ifelse(gene_exp$exp > cutoff, "Low", "High")

# ===================== 【关键】完全和你上面代码一样的分组设置 
gene_exp$risk_group <- factor(
  gene_exp$risk_group, 
  levels = c("Low", "High"), 
  labels = c("Low risk", "High risk")
)

# ===================== 【完全复刻】拟合模型、HR、事件数 
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

# 配色（和你完全一样）
col_pal <- c("Low risk"="#009FC3", "High risk"="#B30437")

# ===================== 【完全一模一样】绘图代码
p <- fit %>%
  ggsurvfit(linewidth = 0.8) +
  # 风险表
  add_risktable(
    risktable_height = 0.28,
    risktable_stats = c("{n.risk}"),
    size = 4.5
  ) +
  # 删失标记
  add_censor_mark(shape = 1, size = 2, stroke = 1) +
  # Log-rank P值标注
  add_pvalue(
    location = "annotation",
    x = max(gene_exp$OS.time)*0.95, y = 0.22,
    hjust = 1, size = 4.2,
    caption = "Log-rank p = {p.value}"
  ) +
  # HR值标注
  annotate("text",
           x = max(gene_exp$OS.time)*0.95, y = 0.32,
           label = paste0("HR = ",hr_val," (95%CI: ",hr_low,"-",hr_high,")"),
           hjust = 1, size = 4.2) +
  # 标题
  labs(
    title = paste0(gene, " - Overall Survival"),
    x = "Time (months)", 
    y = "Survival probability (%)"
  ) +
  scale_x_continuous(expand = c(0.03,0)) +
  scale_y_continuous(limits = c(0,1.05), labels = percent_format()) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  # 主题（和你100%一样）
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

# 合并一页（关键：解决两页PDF）
p <- p + patchwork::plot_layout(nrow = 2, heights = c(3, 1))

# 保存
ggsave(
  "./10_KM生存曲线/单基因——CLEC11A_KM.pdf",
  # "./10_KM生存曲线/单基因——IRF1_KM.pdf",
  # "./10_KM生存曲线/单基因——Enh092845_KM.pdf",
  plot = p,
  width = 7, height = 5.5  
)
