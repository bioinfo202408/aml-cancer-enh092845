# 参考：https://figureya.online/FigureYa18oncoplot_update/FigureYa18oncoplot_update.html
# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./14_突变瀑布图·_oncoplot/"
dir.create(out_dir, showWarnings = F, recursive = T)



########################################################################
# FigureYa oncoplot for TCGA-LAML (Acute Myeloid Leukemia)
# 复刻 Cell 级瀑布图 + 自定义临床信息 + 精美配色
########################################################################

# 加载包
library(TCGAbiolinks)
library(maftools)
library(RColorBrewer)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)


#=======================================================================
#  Step 1: 下载 TCGA-LAML 临床数据
#=======================================================================
clinical <- GDCquery(
  project = "TCGA-LAML", 
  data.category = "Clinical", 
  data.type = "Clinical Supplement", 
  data.format = "BCR XML"
)
GDCdownload(clinical)
cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient")
colnames(cliquery)[1] <- "Tumor_Sample_Barcode"

# 提取 LAML 关键临床特征（FAB分型、性别、染色体风险）
# 这些是 LAML 文章最常用的注释
cliquery$FAB <- cliquery$leukemia_french_american_british_morphology_code
cliquery$Gender <- cliquery$gender
cliquery$Chr_Risk <- cliquery$cytogenetic_abnormalities
cliquery$Age <- cliquery$age
cliquery$Race <- cliquery$race_list

#=======================================================================
#  Step 2: 下载 TCGA-LAML 突变数据 (MAF)
#=======================================================================
query_maf <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_maf)
maf_data <- GDCprepare(query_maf)


#=========================================
# 读取你已经下载好的 MAF 数据
#=========================================
maf <- read.maf(maf = maf_data, clinicalData = cliquery, isTCGA = TRUE)

#=========================================
# 清洗 FAB（完全匹配你的数据）
#=========================================
cliquery$FAB <- cliquery$leukemia_french_american_british_morphology_code
cliquery$FAB <- gsub(" Undifferentiated", "", cliquery$FAB)  # M0 Undifferentiated → M0
cliquery$FAB[cliquery$FAB == "Not Classified"] <- "Unknown"

#=========================================
# 清洗 Gender
#=========================================
cliquery$Gender <- toupper(cliquery$gender)


# 读取之前risk分组信息
risk <- read.csv("../Outdata/7.Prognostic analyse/8_train_lasso_risk.csv", stringsAsFactors = FALSE)
cliquery <- cliquery[cliquery$Tumor_Sample_Barcode %in% risk$sample,]
# 读取eln分组信息
eln <- read.csv("../TCGA_data/TCGA-ELN.csv")
# 先把 eln 变成「ID -> ELN」的字典（快速匹配）
eln_dict <- setNames(eln$ELN, eln$ID)
# 按样本 ID 匹配，自动填充到 cliquery
cliquery$Eln <- eln_dict[cliquery$Tumor_Sample_Barcode]
# 给年龄分组
cliquery$Age <- ifelse(cliquery$age_at_initial_pathologic_diagnosis >= 65, ">=65", "<65")


# 分高低风险组
high_sample <- risk$sample[risk$risk == 1]
low_sample <- risk$sample[risk$risk == 0]
cliquery_high <- cliquery[cliquery$Tumor_Sample_Barcode %in% high_sample,]
cliquery_low <- cliquery[cliquery$Tumor_Sample_Barcode %in% low_sample,]
# 只保留我们需要的列（避免多余空列报错）
cliquery_high <- cliquery_high[, c("Tumor_Sample_Barcode","Gender","FAB","Eln","Age")]
cliquery_low <- cliquery_low[, c("Tumor_Sample_Barcode","Gender","FAB","Eln","Age")]

# 去除没有eln的行-样本
cliquery_high <- cliquery_high[!is.na(cliquery_high$Eln) & cliquery_high$Eln != "'--", ]
table(cliquery_high$Eln, useNA = "always")
cliquery_low <- cliquery_low[!is.na(cliquery_low$Eln) & cliquery_low$Eln != "'--", ]
table(cliquery_low$Eln, useNA = "always")
# Adverse    Favorable Intermediate         <NA> 
#   10           22           36            0 

#=========================================
# 重新构建 maf
#=========================================
maf_high <- read.maf(maf = maf_data, clinicalData = cliquery_high, isTCGA = TRUE)
maf_low <- read.maf(maf = maf_data, clinicalData = cliquery_low, isTCGA = TRUE)

#=========================================
# 配色（完全匹配你的数据）
#=========================================
# 突变类型颜色
mut_colors <- RColorBrewer::brewer.pal(10, "Paired")
names(mut_colors) <- c(
  'Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
  'Frame_Shift_Ins','In_Frame_Ins','Splice_Site',
  'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit'
)

# 临床注释颜色
anno_colors <- list(
  Gender = c("MALE" = "#4682B4", "FEMALE" = "#FF69B4"),
  FAB = c("M0"="#F94144","M1"="#F3722C","M2"="#F9C74F","M3"="#90BE6D",
          "M4"="#43AA8B","M5"="#277DA1","M6"="#6A4C93","M7"="#F9844A","Unknown"="#999999"),
  Eln = c("Favorable"="#38B000","Intermediate"="#FFBE0B","Adverse"="#D00000"),
  Age = c("<65"="#3A3E96", ">=65"="#AD3D3E")
)


#=========================================
# 按 CytRisk 分组排序（文章标准图）
#=========================================
pdf("./14_突变瀑布图_oncoplot//High_risk_LAML_GroupBy_CytRisk.pdf", width=16, height=11)
oncoplot(
  maf = maf_high,
  top = 30,
  clinicalFeatures = c("Gender","FAB","Eln","Age"),
  sortByAnnotation = TRUE,
  colors = mut_colors,
  annotationColor = anno_colors,
  borderCol = "white",
  fontSize = 1.1
)
dev.off()

pdf("./14_突变瀑布图_oncoplot//Low_risk_LAML_GroupBy_CytRisk.pdf", width=16, height=11)
oncoplot(
  maf = maf_low,
  top = 30,
  clinicalFeatures = c("Gender","FAB","Eln","Age"),
  sortByAnnotation = TRUE,
  colors = mut_colors,
  annotationColor = anno_colors,
  borderCol = "white",
  fontSize = 1.1
)
dev.off()






save.image("./14_突变瀑布图_oncoplot/AML.Rdata")


