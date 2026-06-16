# I. Silhouette coefficients (Brief supplementary version)
# 1. Formula: \(s_i=\frac{b_i-a_i}{\max(a_i,b_i)}\)
# \(a_i\): Average distance from sample i to other samples within the same group (Normal / AML);
# \(b_i\): Average distance from sample i to all samples in the other group;
# Overall score = mean value of all \(s_i\), range \([-1,1]\); higher value indicates better group separation.

# II. JS-score
# Code repository: /home/weili/Project/AML/human/AML_combined_analyse/JScore_code




# ==============================================================
# Step0 Environment initialization & path configuration
# ==============================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plotting_Code/")
cat("Initial working directory: ", getwd(), "\n")
# Load global configuration
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_folder <- "./17_JS-score_Silhouette_coefficients//"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("Current output directory: ", getwd(), "\n")



# ========== JS-score Half-violin plot - Full dataset ==========
library(ggplot2)
library(dplyr)
library(gghalves)
library(ggsignif)

rna_files <- list(
  miRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_miRNA.txt",
  mRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_mRNA.txt",
  eRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_eRNA.txt",
  lncRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_lncRNA.txt"
)

rna_data <- lapply(names(rna_files), function(type){
  df <- read.table(rna_files[[type]])
  df$RNA_ID <- rownames(df)
  df %>%
    dplyr::select(RNA_ID, JSscore) %>%
    mutate(JSscore = as.numeric(JSscore), RNA_Type = type)
}) %>% bind_rows()

rna_data <- rna_data[-1,]

# Data cleaning
# clean_data <- rna_data %>%
#   filter(!is.na(JSscore)) %>%
#   filter(JSscore > 0 & JSscore < 1) %>%
#   group_by(RNA_Type) %>%
#   mutate(
#     q10 = quantile(JSscore,0.1,na.rm=T),
#     q90 = quantile(JSscore,0.9,na.rm=T)
#   ) %>%
#   filter(JSscore >= q10 & JSscore <= q90) %>%
#   ungroup()

clean_data <- rna_data


rna_color <- c("mRNA"="#AD3D3E","eRNA"="#96CEB4","lncRNA"="#4ECDC4","miRNA"="#45B7D1")
# Fix X-axis order: mRNA, eRNA, lncRNA, miRNA
clean_data$RNA_Type <- factor(clean_data$RNA_Type,
                              levels = c("mRNA","eRNA","lncRNA","miRNA"))
# Comparison sets: three ncRNA types vs mRNA respectively
comp_list <- list(c("eRNA","mRNA"),c("lncRNA","mRNA"),c("miRNA","mRNA"))


# Generate plot
pdf("./RNA_JSscore_HalfViolin_Box_Jitter_ncRNA_vs_mRNA.pdf",width=7,height=6)
p <- ggplot(clean_data,aes(x=RNA_Type,y=JSscore,color=RNA_Type,fill=RNA_Type))+
  geom_half_violin(side="r",trim=F,alpha=0.7,width=0.5,position=position_nudge(x=0.12))+
  geom_boxplot(width=0.22,outlier.shape=NA,alpha=0.7,position=position_nudge(x=-0.1))+
  stat_summary(fun="mean",geom="point",shape=20,size=3,color="black")+
  geom_jitter(width=0.1,size=0.8,alpha=0.35)+
  geom_signif(
    comparisons = comp_list,
    test = "wilcox.test",
    map_signif_level=T,
    tip_length=0,
    textsize=4
  )+
  scale_fill_manual(values=rna_color)+
  scale_color_manual(values=rna_color)+
  labs(x=NULL,y="Jensen-Shannon score (JS-score)")+
  theme_bw(base_size =15)+
  theme(
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x=element_text(face="bold",size=13),
    axis.text.y=element_text(size=12),
    axis.title.y=element_text(face="bold",size=14)
  )
print(p)
dev.off()



# ==========【Key: Batch Wilcoxon test, print P-values to console
cat("==================== Wilcox rank-sum test result ====================\n")
for(cc in comp_list){
  g1 <- cc[1]
  g2 <- cc[2]
  val1 <- clean_data$JSscore[clean_data$RNA_Type == g1]
  val2 <- clean_data$JSscore[clean_data$RNA_Type == g2]
  wt <- wilcox.test(val1, val2)
  cat(sprintf("%s vs %s | P-value = %.2e\n",g1,g2,wt$p.value))
}
cat("=====================================================================\n")


# eRNA vs mRNA | P-value = 0.00e+00
# lncRNA vs mRNA | P-value = 0.00e+00
# miRNA vs mRNA | P-value = 2.82e-271






# ========== JS-score Half-violin plot - Random subsample of 1000 genes per group ==========
# Randomly sample 1000 genes for plotting
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plotting_Code/")
cat("Initial working directory: ", getwd(), "\n")
# Load global configuration
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_folder <- "./17_JS-score_Silhouette_coefficients//"
dir.create(out_folder, showWarnings = F, recursive = T)
setwd(out_folder)
cat("Current output directory: ", getwd(), "\n")

library(ggplot2)
library(dplyr)
library(gghalves)
library(ggsignif)

# File paths
rna_files <- list(
  miRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_miRNA.txt",
  mRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_mRNA.txt",
  eRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_eRNA.txt",
  lncRNA = "../../Outdata/1.rawdata/3.JScore_filter/4.2.JSscore_H_lncRNA.txt"
)

# Read data
rna_data <- lapply(names(rna_files), function(type){
  df <- read.table(rna_files[[type]])
  df$RNA_ID <- rownames(df)
  df %>%
    dplyr::select(RNA_ID, JSscore) %>%
    mutate(JSscore = as.numeric(JSscore), RNA_Type = type)
}) %>% bind_rows()
rna_data <- rna_data[-1,]

# Data cleaning
clean_data <- rna_data %>%
  filter(!is.na(JSscore)) %>%
  filter(JSscore > 0 & JSscore < 1) %>%
  group_by(RNA_Type) %>%
  mutate(
    q10 = quantile(JSscore,0.1,na.rm=T),
    q90 = quantile(JSscore,0.9,na.rm=T)
  ) %>%
  filter(JSscore >= q10 & JSscore <= q90) %>%
  ungroup()
# clean_data <- rna_data

# ==========1. Fix X-axis order: mRNA as the first category
clean_data$RNA_Type <- factor(clean_data$RNA_Type,
                              levels = c("mRNA","eRNA","lncRNA","miRNA"))

# ==========2. Randomly sample 1000 genes per RNA type for jitter points
set.seed(123)
sub_data <- clean_data %>%
  group_by(RNA_Type) %>%
  slice_sample(n = 1000) %>%
  ungroup()

# ==========Color palette + comparison sets: three ncRNA types vs mRNA
rna_color <- c("mRNA"="#FF6B6B","eRNA"="#96CEB4","lncRNA"="#4ECDC4","miRNA"="#45B7D1")
comp_list <- list(c("eRNA","mRNA"),c("lncRNA","mRNA"),c("miRNA","mRNA"))

# Generate plot
pdf("./RNA_JSscore_HalfViolin_SubSample1000.pdf",width=7,height=6)
p <- ggplot(clean_data,aes(x=RNA_Type,y=JSscore,color=RNA_Type,fill=RNA_Type))+
  # Right half violin (full dataset)
  geom_half_violin(side="r",trim=F,alpha=0.7,width=0.5,position=position_nudge(x=0.12))+
  # Boxplot (full dataset)
  geom_boxplot(width=0.22,outlier.shape=NA,alpha=0.7,position=position_nudge(x=-0.1))+
  # Mean value point (full dataset)
  stat_summary(fun="mean",geom="point",shape=20,size=3,color="black")+
  # Jitter points only using subsampled data sub_data
  geom_jitter(data=sub_data,width=0.1,size=0.8,alpha=0.35)+
  # Significance annotation
  geom_signif(
    comparisons = comp_list,
    test = "wilcox.test",
    map_signif_level=T,
    tip_length=0,
    textsize=4
  )+
  scale_fill_manual(values=rna_color)+
  scale_color_manual(values=rna_color)+
  labs(x=NULL,y="Jensen-Shannon score (JS-score)")+
  theme_bw(base_size =15)+
  theme(
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x=element_text(face="bold",size=13),
    axis.text.y=element_text(size=12),
    axis.title.y=element_text(face="bold",size=14)
  )
print(p)
dev.off()


# ========== Calculate Silhouette coefficients for four RNA types in batch ==========
library(data.table)
library(cluster)   # Must be loaded in advance to call silhouette()
library(dplyr)
library(ggplot2)

# Wrapper function for expression matrix reading
read_rna_exp <- function(path){
  dt = fread(path)
  gn = dt$V1
  exp_df = as.data.frame(dt[,-1])
  rownames(exp_df) = gn
  return(exp_df)
}

# Read each RNA expression matrix
mRNA   <- read_rna_exp("../../Outdata/5.all_data_harmony/combat_all_mRNA.csv")
lncRNA <- read_rna_exp("../../Outdata/5.all_data_harmony/combat_all_lncRNA.csv")
miRNA  <- read_rna_exp("../../Outdata/5.all_data_harmony/combat_all_miRNA.csv")
eRNA   <- read_rna_exp("../../Outdata/5.all_data_harmony/combat_all_eRNA.csv")

# Integrate into list
exp_list <- list(
  mRNA=mRNA,
  lncRNA=lncRNA,
  miRNA=miRNA,
  eRNA=eRNA
)

# Sample grouping information
group_info  <- read.csv("../../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv",row.names=1)
group_info$Group <- ifelse(group_info$Group==1,"Tumor","Normal")
sample_group <- factor(group_info$Group)

group_df_all <- data.frame(
  sample = rownames(group_info),
  grp    = as.character(sample_group)
)

# Loop to compute silhouette coefficient for each RNA type
sil_result <- data.frame()
for(rna_name in names(exp_list)){
  mat <- exp_list[[rna_name]]
  mat_t <- t(mat)
  mat_t <- na.omit(mat_t)
  keep_sample <- rownames(mat_t)
  
  sub_group <- group_df_all[group_df_all$sample %in% keep_sample,"grp"]
  grp_int <- as.integer(factor(sub_group))
  
  dist_mat <- dist(mat_t,method="euclidean")
  sil_obj <- silhouette(grp_int,dist_mat)
  mean_sil <- mean(sil_obj[,3])
  
  sil_result <- rbind(sil_result,
                      data.frame(RNA_Type=rna_name,Silhouette_Coefficient=mean_sil))
}

print(sil_result)
# RNA_Type Silhouette_Coefficient
# 1    miRNA             0.03600181
# 2   lncRNA             0.02319722
# 3     mRNA             0.01396005
# 4     eRNA             0.03156884
write.csv(sil_result,".//Silhouette_allRNA_result.csv",row.names=F)

# Bar plot for silhouette score
sil_result$RNA_Type <- factor(sil_result$RNA_Type,levels=c("mRNA","eRNA","lncRNA","miRNA"))
rna_color <- c("mRNA"="#AD3D3E","eRNA"="#96CEB4","lncRNA"="#4ECDC4","miRNA"="#45B7D1")
comp_list <- list(c("eRNA","mRNA"),c("lncRNA","mRNA"),c("miRNA","mRNA"))


# Sort ascending: from smallest to largest silhouette value
sil_sort <- sil_result[order(sil_result$Silhouette_Coefficient), ]
# Fix factor order using sorted rows
sil_sort$RNA_Type <- factor(sil_sort$RNA_Type, levels = sil_sort$RNA_Type)


pdf(".//Silhouette_Barplot.pdf",width=5.5,height=5)
p <- ggplot(sil_sort,aes(x=RNA_Type,y=Silhouette_Coefficient,fill=RNA_Type))+
  geom_col(width=0.65)+
  geom_text(aes(label=sprintf("%.3f",Silhouette_Coefficient)),vjust=-0.4,size=4.2)+
  ggsignif::geom_signif(comparisons=comp_list,map_signif_level=T,tip_length=0,textsize=4)+
  scale_fill_manual(values=rna_color)+
  labs(x=NULL,y="Mean Silhouette coefficient")+
  theme_bw(base_size=14)+
  theme(legend.position="none",panel.grid=element_blank(),
        axis.text.x=element_text(face="bold",size=12),
        axis.title.y=element_text(face="bold",size=13))
print(p)
dev.off()
