# Drawing reference: https://mp.weixin.qq.com/s/zI7VN5Q_yYDkzDE_KGhTXg
# Requires completed GO/KEGG enrichment analysis results


# Clear environment and load required packages
rm(list = ls())
gc() # Clean unused objects in memory

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【Initialization】Current working directory: ", getwd(), "\n\n")

# Check temporary directory
Sys.getenv("TMPDIR") 

dir.create("./9_GO_KEGG_enrichment_analysis/")



# ==============================================
# 1. Load packages
# ==============================================
library(ggplot2)
library(dplyr)
library(patchwork)


# ==============================================
# ######KEGG Enrichment Analysis####
# ==============================================
# ==============================================
# 2. Construct your KEGG data frame (directly copy and run)
# ==============================================
kegg_data <- read.csv("./9_GO_KEGG_enrichment_analysis/KEGG_Enh092845.csv")
kegg_data <- read.csv("./9_GO_KEGG_enrichment_analysis/mRNA_enriched_KEGG_terms_related_to_prognostic_ncRNA.csv")


# # Retain pathways of interest
# keypathway <- c("Oxidative phosphorylation","Lysosome biogenesis","Phagosome","Adherens junction",
#                 "Insulin resistance","Thermogenesis","Salmonella infection",
#                 "Human cytomegalovirus infection","Vibrio cholerae infection","Parkinson disease")
# kegg_data <- kegg_data[kegg_data$Description %in% keypathway,]


# Select Top10 (sorted by Count descending)
kegg_data <- kegg_data %>%
  arrange(desc(Count)) %>%
  head(10)


# Calculate -log10(p.adjust)
kegg_data <- kegg_data %>%
  mutate(logP = -log10(p.adjust))

# Convert to factor to keep vertical order consistent with input data
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# 3. Color palette (consistent with GO plot)
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# ==============================================
# 4. Left: Lollipop plot (Count)
# ==============================================
p1 <- ggplot(kegg_data, aes(x = Count, y = subcategory)) +
  geom_segment(aes(y = subcategory, yend = subcategory, x = 0, xend = Count),
               linewidth = 0.8, color = "grey80") +
  geom_point(aes(color = subcategory), size = 6) +
  scale_color_manual(values = Color) +
  scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
  labs(title = "Cell marker numbers", x = NULL, y = NULL) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 5. Left: Bubble plot
# ==============================================
kegg_data$size <- as.numeric(sapply(strsplit(kegg_data$GeneRatio, "/"), "[", 1))
p3 <- ggplot(kegg_data, aes(x = 0, y = subcategory)) +  
  geom_point(aes(color = subcategory, size = size)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    panel.border = element_blank(),    
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_text(hjust = 1,color = "black", face = "plain"),
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "plain")
  )


# Alternative code if above plot has overlapping subcategory labels
# # ==============================================
# # 4. Left: Lollipop plot (each term independent, no merge)
# # ==============================================
# # Create unique ID for each row to avoid automatic grouping
# kegg_data$id <- 1:nrow(kegg_data)
# 
# p1 <- ggplot(kegg_data, aes(x = Count, y = id)) +
#   geom_segment(aes(y = id, yend = id, x = 0, xend = Count),
#                linewidth = 0.8, color = "grey80") +
#   geom_point(aes(color = subcategory), size = 6) +
#   
#   # Force Y-axis labels to display subcategory
#   scale_y_continuous(breaks = 1:nrow(kegg_data),
#                      labels = kegg_data$subcategory) +
#   
#   scale_color_manual(values = Color) +
#   scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
#   labs(title = "Cell marker numbers", x = NULL, y = NULL) +
#   theme_bw(base_size = 18) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
#     axis.text.x = element_text(color = "black", face = "plain"),
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(size = 0.5, fill = NA)
#   )
# # ==============================================
# # 5. Left: Bubble plot (each row independent, Y-axis labeled with subcategory)
# # ==============================================
# kegg_data$size <- as.numeric(sapply(strsplit(kegg_data$GeneRatio, "/"), "[", 1))
# 
# # Use unique ID as Y-axis to force individual rows
# kegg_data$id <- 1:nrow(kegg_data)
# 
# p3 <- ggplot(kegg_data, aes(x = 0, y = id)) +  
#   geom_point(aes(color = subcategory, size = size)) +  
#   scale_color_manual(values = Color, guide = "none") +  
#   scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
#   scale_y_continuous(breaks = 1:nrow(kegg_data), labels = kegg_data$subcategory) +
#   scale_x_discrete(expand = c(0,0)) +  
#   labs(title = "", x = NULL, y = NULL) +  
#   theme_bw(base_size = 18) +  
#   theme(    
#     plot.title = element_text(hjust = 0, size = 18, face = "plain"),
#     panel.border = element_blank(),    
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     axis.ticks = element_blank(),    
#     panel.grid = element_blank(),    
#     legend.position = "right",
#     legend.text = element_text(face = "plain"),
#     legend.title = element_text(face = "plain")
#   )




# ==============================================
# 6. Right: Bar plot (-log10(p.adjust))
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "KEGG enrichment terms", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 7. Combine plots and save PDF
# ==============================================
p <- p1 + p2 + plot_layout(widths = c(1, 1.8))
# Lollipop plot + right bar plot
ggsave("./9_GO_KEGG_enrichment_analysis/KEGG_lollipop_barplot.pdf",
       plot = p,
       width = 10,
       height = 6,
       device = "pdf")

# Bubble plot + right bar plot
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")
ggsave("./9_GO_KEGG_enrichment_analysis/KEGG_bubble_barplot.pdf",
       plot = p,
       width = 10,
       height = 6,
       device = "pdf")









# ==============================================
# ######GO Enrichment Analysis######
# ==============================================
kegg_data <- read.csv("./9_GO_KEGG_enrichment_analysis/mRNA_enriched_GO_terms_related_to_prognostic_ncRNA.csv")

# Calculate GeneRatio: convert "33/739" to numeric value e.g. 0.04465
kegg_data$GeneRatio_value <- sapply(kegg_data$GeneRatio, function(x) eval(parse(text = x)))


# Select Top10 (sorted by GeneRatio_value descending)
kegg_data <- kegg_data %>%
  arrange(desc(GeneRatio_value)) %>%
  head(10)


# Calculate -log10(p.adjust)
kegg_data <- kegg_data %>%
  mutate(logP = -log10(p.adjust))

# Convert to factor to keep vertical order consistent with input data
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# 3. Color palette (consistent with GO plot)
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# # ==============================================
# # 4. Left: Lollipop plot (Count)
# # ==============================================
# p1 <- ggplot(kegg_data, aes(x = Count, y = subcategory)) +
#   geom_segment(aes(y = subcategory, yend = subcategory, x = 0, xend = Count),
#                linewidth = 0.8, color = "grey80") +
#   geom_point(aes(color = subcategory), size = 6) +
#   scale_color_manual(values = Color) +
#   scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
#   labs(title = "Cell marker numbers", x = NULL, y = NULL) +
#   theme_bw(base_size = 18) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
#     axis.text.x = element_text(color = "black", face = "plain"),
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(size = 0.5, fill = NA)
#   )

# ==============================================
# 5. Left: Bubble plot
# ==============================================
p3 <- ggplot(kegg_data, aes(x = 0, y = Description)) +  
  geom_point(aes(color = Description, size = Count)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    panel.border = element_blank(),    
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "plain")
  )

# ==============================================
# 6. Right: Bar plot (-log10(p.adjust))
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "GO enrichment terms", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 7. Combine plots and save PDF
# ==============================================
# p <- p1 + p2 + plot_layout(widths = c(1, 1.8))
# ggsave("./9_GO_KEGG_enrichment_analysis/KEGG_lollipop_barplot.pdf",
#        plot = p,
#        width = 10,
#        height = 6,
#        device = "pdf")

# Bubble plot + right bar plot
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")
ggsave("./9_GO_KEGG_enrichment_analysis/GO_bubble_barplot.pdf",
       plot = p,
       width = 8,
       height = 6,
       device = "pdf")




# ==============================================
# ######IRF1 KO Perturbed Genes - GO Enrichment Analysis######
# ==============================================
library(openxlsx)
library(tidyverse)
library(patchwork)

# Load data
kegg_data <- read.xlsx("./9_GO_KEGG_enrichment_analysis/IRF1_KO_virtual_ko_go_result_ENRICH.xlsx")
kegg_data <- kegg_data[-1,]

# Extract numerator from "numerator/denominator" format as Count column
kegg_data$Count <- as.integer(sapply(strsplit(kegg_data$InTerm_InList, "/"), function(x) x[1]))

# Calculate GeneRatio safely, skip malformed entries like "13/-"
kegg_data$GeneRatio_value <- sapply(kegg_data$InTerm_InList, function(x) {
  if(grepl("^\\d+/\\d+$", x)) {
    eval(parse(text = x))
  } else {
    NA
  }
})

# Remove rows with invalid data for clean plotting
kegg_data <- kegg_data[!is.na(kegg_data$GeneRatio_value), ]


# Sort by LogP ascending, take TOP 10
# ======================================================
kegg_data <- kegg_data %>%
  # arrange(desc(LogP)) %>%   # Sort by significance descending
  arrange(LogP) %>%   # Sort ascending (remove desc for ascending order)
  head(10)

# Compute logP (consistent with previous code)
kegg_data <- kegg_data %>%
  mutate(logP = -LogP)

# Factorize Description to keep vertical order Top1 -> Top10
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# Color palette
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# ==============================================
# Left: Bubble plot
# ==============================================
p3 <- ggplot(kegg_data, aes(x = 0, y = Term)) +
  geom_point(aes(color = Term, size = Count)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,30, 10), limits = c(0,30)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    panel.border = element_blank(),    
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "plain")
  )

# ==============================================
# Right: Bar plot (-log10(p.adjust))
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "GO enrichment terms", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )

# ==============================================
# Combine plots & export figure
# ==============================================
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")

ggsave("./9_GO_KEGG_enrichment_analysis/IRF1_KO_GO_TOP10_LogP_bubble_barplot.pdf",
       plot = p,
       width = 12,
       height = 7,
       device = "pdf",
       bg="white")

cat("✅ GO TOP10 plot sorted by LogP generated successfully!\n")
