# Reference: https://figureya.online/FigureYa292HCCsubtype/FigureYa292HCCsubtype.html
# Clear environment
rm(list = ls())
gc()

# Set working directory
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# Create output folder
out_dir <- "./18_NMF/"
dir.create(out_dir, showWarnings = F, recursive = T)

# ===================== NMF based on 1048 signature genes =====================
library(NMF)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

expr <- fread("../Outdata/5.all_data_harmony/5.1_all_expr_train.csv")
expr <- expr[expr$Group %in% 1,]
sample <- expr$Sample
expr <- expr[,!(1:3)]
expr <- as.data.frame(t(expr))
colnames(expr) <- sample
head(expr)

## Key point for NMF: input matrix must be non-negative
## Log-normalized expression contains negative values, shift all values to positive globally
min_val <- min(as.matrix(expr))
if(min_val < 0){
  expr <- expr - min_val + 1e-4 # Global offset to eliminate negative values, tiny offset to avoid zero
}
expr_mat <- as.matrix(expr)
cat("Expression matrix: Genes =",nrow(expr_mat)," Samples =",ncol(expr_mat),"\n")

#=========1. Select optimal rank from 2~6, plot rank selection figure (silhouette coefficient + cophenetic correlation dual line plot)
result <- nmf(
  expr_mat,        
  rank = 2:6,      
  method = "lee",  
  nrun = 15,      
  seed = 123       
)

pdf("./18_NMF/1_Optimal_Cluster_Selection_TCGA.pdf", width = 8, height = 8)
plot(result)
dev.off()

#=========2. Formal NMF decomposition with rank=3
result1 <- nmf(
  expr_mat,        
  rank = 3,        
  method = "lee",  
  nrun = 15,      
  seed = 123       
)

# Sample subtype information
group <- data.frame(
  cluster = predict(result1),
  sample_id = names(predict(result1)),
  row.names = names(predict(result1))
)
table(group$cluster)

#=========3. Consensus matrix heatmap for samples (consensusmap)
pdf("./18_NMF/2_NMF_Consensus_Matrix_Heatmap-K3.pdf", width = 5, height = 5)
consensusmap(
  result1, 
  labRow = NA,          
  labCol = NA,          
  annCol = group,       
  annColors = list(     
    cluster = c("1" = "#456990", "2" = "#EF767A", "3" = "#48C0AA")
  )
)
dev.off()

#=========6. Extract signature genes for each subtype
# 1. Extract W matrix: gene × 3 subtype weight matrix
W <- basis(result1)
# 2. Assign each gene to subtype with maximum weight (1/2/3)
gene_cluster_id <- apply(W, 1, which.max)

# 3. Split gene symbols by subtype
gene_all <- rownames(expr_mat)
gene_list_all <- split(gene_all, factor(gene_cluster_id,levels = c(1,2,3)))
names(gene_list_all) = c("Cluster1","Cluster2","Cluster3")

# Verify total gene count equals input 1048 genes
sum(sapply(gene_list_all,length)) 

cluster_gene1 <- gene_list_all$Cluster1
cluster_gene2 <- gene_list_all$Cluster2
cluster_gene3 <- gene_list_all$Cluster3

write.table(cluster_gene1,"./18_NMF/cluster_gene1.txt", sep = ",", quote = FALSE, row.names = FALSE)
write.table(cluster_gene2,"./18_NMF/cluster_gene2.txt", sep = ",", quote = FALSE, row.names = FALSE)
write.table(cluster_gene3,"./18_NMF/cluster_gene3.txt", sep = ",", quote = FALSE, row.names = FALSE)

# ==================== GO/KEGG enrichment analysis based on subtype signature genes =====================
# Load enrichment & visualization packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(patchwork)

# Create folder for enrichment output
out_go_dir <- "./18_NMF/GO_Enrich_Cluster/"
dir.create(out_go_dir,recursive = T,showWarnings = F)

# Read signature genes of each cluster
cluster1_gene <- read.table("./18_NMF/cluster_gene1.txt")
gene <- cluster1_gene$V1[-1]
cluster2_gene <- read.table("./18_NMF/cluster_gene2.txt")
gene <- cluster2_gene$V1[-1]
cluster3_gene <- read.table("./18_NMF/cluster_gene3.txt")
gene <- cluster3_gene$V1[-1]

# Convert gene SYMBOL to ENTREZ ID
gene_map <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
entrez_vect <- gene_map$ENTREZID %>% unique()

# GO-ALL enrichment
go_res <- enrichGO(gene = entrez_vect,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   readable = T)

# Convert result to dataframe
go_df <- as.data.frame(go_res)

# KEGG enrichment analysis (run KEGG if GO returns no significant terms)
kegg_res <- enrichKEGG(gene = entrez_vect,
                       organism = "hsa", # Human hsa database
                       pvalueCutoff = 1)
kegg_df <- as.data.frame(kegg_res)
kegg_sig <- subset(kegg_df,pvalue <0.05)

# Save raw enrichment table for current cluster
# write.csv(go_df, "./18_NMF/GO_Enrich_Cluster/cluster1_GO_all_AllResult.csv")
# write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster1_KEGG_AllResult.csv")

# write.csv(go_df, "./18_NMF/GO_Enrich_Cluster/cluster2_GO_all_AllResult.csv")
# write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster2_KEGG_AllResult.csv")

# Cluster3 has no significant GO results, only KEGG terms available
# write.csv(go_df, "./18_NMF/GO_Enrich_Cluster/cluster3_GO_all_AllResult.csv")
write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster3_KEGG_AllResult.csv")

go_df <- kegg_sig

# Extract top 10 enriched pathways sorted by gene count
go_top10 <- go_df %>% arrange(desc(Count)) %>% head(10)
go_top10$logP <- -log10(go_top10$p.adjust)
go_top10$Description <- factor(go_top10$Description,levels = rev(go_top10$Description))
range(go_top10$Count) # Adjust bubble size scale according to value range

# Color scheme for NMF subtypes
nmf_col <- c("Cluster1"="#456990","Cluster2"="#EF767A","Cluster3"="#48C0AA")
col_set <- "#48C0AA"

# Left panel: Bubble plot
p_bubble <- ggplot(go_top10,aes(x=0,y=Description))+
  geom_point(aes(size=Count),color=col_set)+
  scale_size_continuous(name = "Marker Gene Num.", breaks = seq(0,5, 1), limits = c(0,5)) + 
  scale_x_discrete(expand=c(0,0))+
  labs(x=NULL,y=NULL)+
  theme_bw(base_size = 13)+
  theme(panel.border = element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())

# Right panel: Bar plot of -log10(p.adjust)
# =========Fix: add stat="identity"=========
p_bar <- ggplot(go_top10,aes(y=Description))+
  geom_bar(aes(x=logP),fill=col_set,width=0.5,alpha=0.7,stat="identity")+
  geom_text(aes(x=0.05,label=Description),hjust=0,size=3.8)+
  labs(x="-log10(adjusted P-value)",title=paste0("Cluster3 KEGG Top10 Enrichment"))+
  theme_bw(base_size = 13)+
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank())

# Combine plots: left bubble + right bar
p_all <- p_bubble + p_bar + plot_layout(widths = c(0.15,1.2))

# Save PDF figure
ggsave("./18_NMF/GO_Enrich_Cluster/Cluster3_KEGG_Top10_Enrich.pdf",
       plot=p_all,width=11,height=5.5)
