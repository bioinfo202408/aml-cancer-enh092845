# Reference: https://mp.weixin.qq.com/s/9bugJ-uy63PrpvIC6Ha0GA  

load("./20260407_analysis_result/rawdata/Step8.4.CLEC11A_expression_in_different_cell_subpopulations.Rdata")

# Clear environment and load required packages
rm(list = ls())
gc() # Clean unused objects in memory

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

dir.create("./2_Bar_Swarm_Line_Plots/")

########################## Bar plot for RNA counts before and after filtering #################################
library(tidyverse)
library(scales)

# Data preparation (consistent with original raw data)
df <- tibble(
  RNA_type = rep(c("mRNA","miRNA","lncRNA","eRNA"), each = 2),
  stage    = rep(c("Rawdata","filtered"), 4),
  count    = c(20046, 9274, 10473, 707, 192138, 4364, 524986, 25757)
) %>%
  mutate(
    stage = factor(stage, levels = c("Rawdata", "filtered")),
    RNA_type = factor(RNA_type, levels = c("mRNA","miRNA","lncRNA","eRNA"))
  )

# Custom color palette corresponding to 4 RNA types
rna_colors <- c(
  "mRNA"    = "#839FBF",   # Blue
  "miRNA"   = "#93A89B",   # Gray-green
  "lncRNA"  = "#CAE0CA",   # Light green
  "eRNA"    = "#FFC89D"    # Light orange
)

pdf("./2_Bar_Swarm_Line_Plots/1.RNA_filter_comparison_plot2.pdf",width = 9, height = 5)

ggplot(df, aes(x = RNA_type, y = count, fill = interaction(RNA_type, stage), alpha = stage)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_y_log10(labels = comma) +
  
  # Core: assign fixed colors to each RNA type, distinguish Raw/filtered by transparency
  scale_fill_manual(
    values = c(
      "mRNA.Rawdata" = "#839FBF",
      "mRNA.filtered" = "#839FBF",
      "miRNA.Rawdata" = "#93A89B",
      "miRNA.filtered" = "#93A89B",
      "lncRNA.Rawdata" = "#CAE0CA",
      "lncRNA.filtered" = "#CAE0CA",
      "eRNA.Rawdata" = "#FFC89D",
      "eRNA.filtered" = "#FFC89D"
    )
  ) +
  scale_alpha_manual(
    values = c("Rawdata" = 1, "filtered" = 0.6),  # Filtered bars with lighter transparency
    guide = guide_legend(title = "Data Stage")
  ) +
  
  geom_text(
    aes(label = comma(count)),
    position = position_dodge(width = 0.9),
    vjust = -0.3, size = 3.5
  ) +
  
  labs(
    title = "RNA Type: Raw Counts vs Filtered Counts",
    x = "RNA Type", 
    y = "Number of RNA Genes (log10 scale)",
    fill = "RNA Type",
    alpha = "Data Stage"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  )

dev.off()

########################## Hub ncRNA comprehensive score bar plot #################################
library(ggplot2)
library(dplyr)

# 1. Import and organize data
hub_score <- read.csv("../Outdata/25_AML_key_gene_identification_new/3_top10_ncRNA_hubscore.csv")
hub_score <- hub_score[, c("node", "hub_score")]

# 2. Sort by hub_score ascending for horizontal bar plot display
hub_score_plot <- hub_score %>% 
  arrange(hub_score) %>%
  mutate(node = factor(node, levels = node))

# 3. Custom color palette (10 bars reuse 8 predefined colors cyclically)
color_palette <- c(
  "#F28C66",
  "#E5719A",
  "#E6C29A",
  "#974C01",
  "#F7D354",
  "#89C451",
  "#D673A3",
  "#7B9CC2",
  "#999999",
  "#4CB8A1"
)

# Recycle colors to match total bar count
color_palette <- rep(color_palette, length.out = nrow(hub_score_plot))

# 4. Plot horizontal bar chart without grid lines
p <- ggplot(hub_score_plot, aes(x = node, y = hub_score)) +
  geom_col(aes(fill = node), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = color_palette) +
  coord_flip() +  # Convert to horizontal bar plot
  labs(
    x = "Key ncRNA",
    y = "Hub Score",
    title = "Top 10 Key ncRNA Ranked by Hub Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(linewidth = 0.8, color = "black")
  )

# 5. Save figure
pdf("./2_Bar_Swarm_Line_Plots/TOP10_ncRNA_barplot.pdf", width = 8, height = 5)
print(p)
dev.off()

########################## Motif enrichment bar plot #################################
# Load packages
library(ggplot2)
library(dplyr)

# ===================== 1. Import raw data (original code unchanged) =====================
tf_enrich <- read.table(
  file = "../motif_analysis/3.hub_gene_motif_enrich/hub_gene_motif_enriched_TF_list_MI0022577.txt",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)
colnames(tf_enrich) <- c("TF_ID", "TF_Name", "p_value")
tf_enrich <- tf_enrich[!duplicated(tf_enrich$TF_Name), ]

tf_enrich <- tf_enrich %>%
  mutate(
    p_value_corrected = ifelse(p_value == 0, 1e-300, p_value),
    log10_p = -log10(p_value_corrected)
  ) %>% arrange(desc(log10_p))

tf_enrich_Enh092845 <- tf_enrich

# AML-associated transcription factor list
AML_TF_list <- c("BD4", "LYL1", "RUNX1","RUNX2", "ELF1", "EVT6", "CEBPA", "MYC", "GATA2", "SPI1", "PU.1", "PU1", "HOXA9", 
                 "MEIS1", "RBMX", "L1", "CEBPA", "NPM1", "CBF", "TP53", "MDM2", "RUNX1T1", "GATA1",
                 "TAL1", "SCL", "WT1", "IKZF1", "IRF1", "ELF1", "TWIST1",
                 "MCL1", "BCL2")

tf_enrich_Enh092845_AML <- tf_enrich_Enh092845[tf_enrich_Enh092845$TF_Name %in% AML_TF_list, ]

# ===================== 2. Plot: horizontal bar chart + 8-color palette + no grid lines =====================
# Standard 8-color palette
color_palette <- c(
  "#4CB8A1", "#E6C29A", "#F7D354", "#999999", 
  "#D673A3", "#7B9CC2", "#F28C66", "#89C451"
)

# ===================== Key step: sort values descending =====================
plot_data <- tf_enrich_Enh092845_AML %>%
  arrange(desc(log10_p)) %>%  # Sort from highest to lowest
  mutate(TF_Name = factor(TF_Name, levels = TF_Name))  # Fix display order

# Recycle colors to match TF count
color_use <- rep(color_palette, length.out = nrow(plot_data))

# ===================== 3. Plot with TF names on X-axis =====================
pdf("./2_Bar_Swarm_Line_Plots/TF_motif_enrichment_barplot.pdf", width=8, height=5)

ggplot(plot_data, aes(x = TF_Name, y = log10_p)) +
  geom_col(aes(fill = TF_Name), 
           color = "black", linewidth=0.5, width=0.7) +
  scale_fill_manual(values = color_use) +
  
  # Significance threshold dashed line
  geom_hline(yintercept = 1.301, 
             linetype="dashed", color="black", linewidth=0.6) +
  
  labs(
    x = "Transcription Factor (Motif)",
    y = "-log10(P-value)",
    title = "Transcription Factor Motif Enrichment Analysis"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=14, face="bold"),
    axis.title = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1, size=10, color="black"),
    axis.text.y = element_text(size=10, color="black"),
    legend.position = "none",
    panel.grid = element_blank()
  )

dev.off()
