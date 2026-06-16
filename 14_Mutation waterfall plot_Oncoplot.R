# Reference: https://figureya.online/FigureYa18oncoplot_update/FigureYa18oncoplot_update.html
# Clear environment
rm(list = ls())
gc()

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output directory
out_dir <- "./14_Mutation_Oncoplot/"
dir.create(out_dir, showWarnings = F, recursive = T)



########################################################################
# FigureYa oncoplot for TCGA-LAML (Acute Myeloid Leukemia)
# Reproduce Cell-level oncoplot + customized clinical annotations + refined color scheme
########################################################################

# Load packages
library(TCGAbiolinks)
library(maftools)
library(RColorBrewer)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)


#=======================================================================
#  Step 1: Download TCGA-LAML clinical data
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

# Extract core clinical features for LAML (FAB subtype, gender, cytogenetic risk)
# These are standard annotations used in LAML publications
cliquery$FAB <- cliquery$leukemia_french_american_british_morphology_code
cliquery$Gender <- cliquery$gender
cliquery$Chr_Risk <- cliquery$cytogenetic_abnormalities
cliquery$Age <- cliquery$age
cliquery$Race <- cliquery$race_list

#=======================================================================
#  Step 2: Download TCGA-LAML mutation data (MAF)
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
# Load pre-downloaded MAF data
#=========================================
maf <- read.maf(maf = maf_data, clinicalData = cliquery, isTCGA = TRUE)

#=========================================
# Clean FAB subtype labels (consistent with dataset)
#=========================================
cliquery$FAB <- cliquery$leukemia_french_american_british_morphology_code
cliquery$FAB <- gsub(" Undifferentiated", "", cliquery$FAB)  # Convert M0 Undifferentiated to M0
cliquery$FAB[cliquery$FAB == "Not Classified"] <- "Unknown"

#=========================================
# Clean Gender labels
#=========================================
cliquery$Gender <- toupper(cliquery$gender)


# Load risk stratification results from previous analysis
risk <- read.csv("../Outdata/7.Prognostic_analysis/8_train_lasso_risk.csv", stringsAsFactors = FALSE)
cliquery <- cliquery[cliquery$Tumor_Sample_Barcode %in% risk$sample,]
# Load ELN cytogenetic classification data
eln <- read.csv("../TCGA_data/TCGA-ELN.csv")
# Create ID-to-ELN mapping dictionary for fast matching
eln_dict <- setNames(eln$ELN, eln$ID)
# Match ELN classification to clinical metadata by sample ID
cliquery$Eln <- eln_dict[cliquery$Tumor_Sample_Barcode]
# Age grouping
cliquery$Age <- ifelse(cliquery$age_at_initial_pathologic_diagnosis >= 65, ">=65", "<65")


# Split samples into high-risk and low-risk groups
high_sample <- risk$sample[risk$risk == 1]
low_sample <- risk$sample[risk$risk == 0]
cliquery_high <- cliquery[cliquery$Tumor_Sample_Barcode %in% high_sample,]
cliquery_low <- cliquery[cliquery$Tumor_Sample_Barcode %in% low_sample,]
# Retain only required columns to avoid errors from redundant empty columns
cliquery_high <- cliquery_high[, c("Tumor_Sample_Barcode","Gender","FAB","Eln","Age")]
cliquery_low <- cliquery_low[, c("Tumor_Sample_Barcode","Gender","FAB","Eln","Age")]

# Remove samples with missing ELN classification
cliquery_high <- cliquery_high[!is.na(cliquery_high$Eln) & cliquery_high$Eln != "'--", ]
table(cliquery_high$Eln, useNA = "always")
cliquery_low <- cliquery_low[!is.na(cliquery_low$Eln) & cliquery_low$Eln != "'--", ]
table(cliquery_low$Eln, useNA = "always")
# Adverse    Favorable Intermediate         <NA> 
#   10           22           36            0 

#=========================================
# Reconstruct MAF objects stratified by risk group
#=========================================
maf_high <- read.maf(maf = maf_data, clinicalData = cliquery_high, isTCGA = TRUE)
maf_low <- read.maf(maf = maf_data, clinicalData = cliquery_low, isTCGA = TRUE)

#=========================================
# Color palette (matched to dataset categories)
#=========================================
# Mutation type colors
mut_colors <- RColorBrewer::brewer.pal(10, "Paired")
names(mut_colors) <- c(
  'Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
  'Frame_Shift_Ins','In_Frame_Ins','Splice_Site',
  'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit'
)

# Clinical annotation colors
anno_colors <- list(
  Gender = c("MALE" = "#4682B4", "FEMALE" = "#FF69B4"),
  FAB = c("M0"="#F94144","M1"="#F3722C","M2"="#F9C74F","M3"="#90BE6D",
          "M4"="#43AA8B","M5"="#277DA1","M6"="#6A4C93","M7"="#F9844A","Unknown"="#999999"),
  Eln = c("Favorable"="#38B000","Intermediate"="#FFBE0B","Adverse"="#D00000"),
  Age = c("<65"="#3A3E96", ">=65"="#AD3D3E")
)


#=========================================
# Sort samples by Cytogenetic Risk (standard plotting strategy for publications)
#=========================================
pdf("./14_Mutation_Oncoplot/High_risk_LAML_GroupBy_CytRisk.pdf", width=16, height=11)
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

pdf("./14_Mutation_Oncoplot/Low_risk_LAML_GroupBy_CytRisk.pdf", width=16, height=11)
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



save.image("./14_Mutation_Oncoplot/AML.Rdata")
