# 参考：:https://figureya.online/FigureYa292HCCsubtype/FigureYa292HCCsubtype.html
# 清空环境
rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")

# 创建输出目录
out_dir <- "./18_NMF/"
dir.create(out_dir, showWarnings = F, recursive = T)

# =====================1048特征基因进行NMF =====================
library(NMF)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)



expr <- fread("../Outdata/5.all_data_harmony/5.1_all_expr_train.csv")
expr <- expr[expr$Group %in% 1,]
smaple <- expr$Sample
expr <- expr[,!(1:3)]
expr <- as.data.frame(t(expr))
colnames(expr) <- smaple
head(expr)


## NMF关键：NMF要求输入矩阵非负，你数据存在负值(校正后log值)，全局最小值平移转正
min_val <- min(as.matrix(expr))
if(min_val < 0){
  expr <- expr - min_val + 1e-4 # 整体偏移消除负数，极小值防0
}
expr_mat <- as.matrix(expr)
cat("表达矩阵：基因",nrow(expr_mat)," 样本",ncol(expr_mat),"\n")



#=========1. 筛选最优rank 2~6，绘制K值筛选图（轮廓系数+cophenetic双折线）
result <- nmf(
  expr_mat,        
  rank = 2:6,      
  method = "lee",  
  nrun = 15,      
  seed = 123       
)

pdf("./18_NMF//1_最佳分簇选择_TCGA.pdf", width = 8, height = 8)
plot(result)
dev.off()

#=========2. rank=3正式NMF分解
result1 <- nmf(
  expr_mat,        
  rank = 3,        
  method = "lee",  
  nrun = 15,      
  seed = 123       
)

# 样本分组信息
group <- data.frame(
  cluster = predict(result1),
  sample_id = names(predict(result1)),
  row.names = names(predict(result1))
)
table(group$cluster)

#=========3. 样本一致性矩阵热图 consensusmap
pdf("./18_NMF//2_NMF一致性矩阵热图-K4.pdf", width = 5, height = 5)
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


#=========6. 提取特征基因
# 1.提取W矩阵：gene×3亚型权重
W <- basis(result1)
# 2.每个基因找权重最大值对应的亚型编号(1/2/3)
gene_cluster_id <- apply(W, 1, which.max)

# 3.按亚型拆分基因名
gene_all <- rownames(expr_mat)
gene_list_all <- split(gene_all, factor(gene_cluster_id,levels = c(1,2,3)))
names(gene_list_all) = c("Cluster1","Cluster2","Cluster3")

# 统计总基因数=总输入基因1048
sum(sapply(gene_list_all,length)) 

clsuter_gene1 <- gene_list_all$Cluster1
clsuter_gene2 <- gene_list_all$Cluster2
clsuter_gene3 <- gene_list_all$Cluster3

write.table(clsuter_gene1,"./18_NMF/clsuter_gene1.txt", sep = ",", quote = FALSE, row.names = FALSE)
write.table(clsuter_gene2,"./18_NMF/clsuter_gene2.txt", sep = ",", quote = FALSE, row.names = FALSE)
write.table(clsuter_gene3,"./18_NMF/clsuter_gene3.txt", sep = ",", quote = FALSE, row.names = FALSE)






# ====================根据每个簇的差异基因进行GO/KEGG富集 =====================
#加载富集&绘图包
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(patchwork)

#创建富集输出文件夹
out_go_dir <- "./18_NMF/GO_Enrich_Cluster/"
dir.create(out_go_dir,recursive = T,showWarnings = F)

# 读取每个簇的特征基因
clsuter1_gene <- read.table("./18_NMF/clsuter_gene1.txt")
gene <- clsuter1_gene$V1[-1]
clsuter2_gene <- read.table("./18_NMF/clsuter_gene2.txt")
gene <- clsuter2_gene$V1[-1]
clsuter3_gene <- read.table("./18_NMF/clsuter_gene3.txt")
gene <- clsuter3_gene$V1[-1]



#基因symbol转ENTREZID
gene_map <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
entrez_vect <- gene_map$ENTREZID %>% unique()

#GO-ALL富集
go_res <- enrichGO(gene = entrez_vect,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   readable = T)

#转成数据框
go_df <- as.data.frame(go_res)

# KEGG富集分析  如果上面的GO没结果，再跑KEGG
kegg_res <- enrichKEGG(gene = entrez_vect,
                       organism = "hsa", #人源hsa
                       pvalueCutoff = 1)
kegg_df <- as.data.frame(kegg_res)
kegg_sig <- subset(kegg_df,pvalue <0.05)

#保存该簇原始富集表格
# write.csv(, "./18_NMF/GO_Enrich_Cluster/cluster1_GO_all_AllResult.csv")
# write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster1_KEGG_AllResult.csv")

# write.csv(go_df, "./18_NMF/GO_Enrich_Cluster/cluster2_GO_all_AllResult.csv")
# write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster3_KEGG_AllResult.csv")

# 3的GO结果没有，只有KEGG
# write.csv(go_df, "./18_NMF/GO_Enrich_Cluster/cluster3_GO_all_AllResult.csv")
write.csv(kegg_sig, "./18_NMF/GO_Enrich_Cluster/cluster3_KEGG_AllResult.csv")



go_df <- kegg_sig
  
  #取TOP10（按Count降序）
  go_top10 <- go_df %>% arrange(desc(Count)) %>% head(10)
  go_top10$logP <- -log10(go_top10$p.adjust)
  go_top10$Description <- factor(go_top10$Description,levels = rev(go_top10$Description))
  range(go_top10$Count)  #根据范围，修改气泡的大小 scale_size_continuous
  
  #绘图配色
  nmf_col <- c("Cluster1"="#456990","Cluster2"="#EF767A","Cluster3","Cluster3"="#48C0AA")
  col_set <- "#456990"
  
  #左：气泡图
  p_bubble <- ggplot(go_top10,aes(x=0,y=Description))+
    geom_point(aes(size=Count),color=col_set)+
    scale_size_continuous(name = "Marker Num.", breaks = seq(0,5, 1), limits = c(0,5)) + 
    scale_x_discrete(expand=c(0,0))+
    labs(x=NULL,y=NULL)+
    theme_bw(base_size = 13)+
    theme(panel.border = element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())
  
  #右：-log10(p.adjust)条形图
  # =========修复：添加 stat="identity"=========
  p_bar <- ggplot(go_top10,aes(y=Description))+
    geom_bar(aes(x=logP),fill=col_set,width=0.5,alpha=0.7,stat="identity")+
    geom_text(aes(x=0.05,label=Description),hjust=0,size=3.8)+
    labs(x="-log10(p.adjust)",title=paste0(cl," GO-BP Top10"))+
    theme_bw(base_size = 13)+
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank())
  
  #拼接图：左气泡+右条形
  p_all <- p_bubble + p_bar + plot_layout(widths = c(0.15,1.2))
  
  #保存PDF
  ggsave(("./18_NMF/GO_Enrich_Cluster/Cluster3_GO_Top10_Enrich.pdf"),
         plot=p_all,width=11,height=5.5)
  
  
  