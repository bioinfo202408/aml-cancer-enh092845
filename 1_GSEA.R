# 清空环境并加载必要的包
rm(list = ls())
gc()

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

### GSEA富集分析
# 加载包
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(cowplot)

dir.create("./画图代码/1_GSEA/")

# ======================1.读取logFC排序文件  ======================
file_name <- "../Outdata/26_虚拟敲除/GMP-LIKE-IRF1基因虚拟敲除显著相关基因.csv"

# 读取差异结果
ko_res <- read.csv(file_name, stringsAsFactors = FALSE)
ko_res <- ko_res[-1,] # 去除自己那行


# 基因名转换
gene_map <- bitr(
  geneID = ko_res$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ko_res <- inner_join(ko_res, gene_map, by = c("gene" = "SYMBOL"))

# 构建 GSEA 基因列表
geneList <- ko_res$log2FC
names(geneList) <- ko_res$ENTREZID  # 这里修复！！！
geneList <- sort(geneList, decreasing = TRUE)

# 去掉无效ID
geneList <- geneList[!is.na(names(geneList))]


# ======================2.开始GSEA富集分析  ======================

# MSigDB C2 (包含 KEGG/Reactome/BioCarta):::::::::::优选跑这个
msig_c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene)  # 必须保留这两列
# 关键！把基因ID转成 字符串格式，防止匹配失败
names(geneList) <- as.character(names(geneList))
gse_c2 <- GSEA(
  geneList  = geneList,
  TERM2GENE = msig_c2,
  scoreType = "pos",
  pvalueCutoff = 1
)
as.data.frame(gse_c2)[, c("ID", "Description", "NES", "pvalue", "p.adjust")]
# REACTOME_INNATE_IMMUNE_SYSTEM                                       REACTOME_INNATE_IMMUNE_SYSTEM                    REACTOME_INNATE_IMMUNE_SYSTEM 2.1834969
# YAGI_AML_FAB_MARKERS                                                         YAGI_AML_FAB_MARKERS                             YAGI_AML_FAB_MARKERS 1.8689319
# SMID_BREAST_CANCER_NORMAL_LIKE_UP                               SMID_BREAST_CANCER_NORMAL_LIKE_UP                SMID_BREAST_CANCER_NORMAL_LIKE_UP 1.8626032
# CHEN_METABOLIC_SYNDROM_NETWORK                                     CHEN_METABOLIC_SYNDROM_NETWORK                   CHEN_METABOLIC_SYNDROM_NETWORK 1.6494176
# KRIGE_RESPONSE_TO_TOSEDOSTAT_6HR_DN                           KRIGE_RESPONSE_TO_TOSEDOSTAT_6HR_DN              KRIGE_RESPONSE_TO_TOSEDOSTAT_6HR_DN 1.4430633
# ROSS_ACUTE_MYELOID_LEUKEMIA_CBF                                   ROSS_ACUTE_MYELOID_LEUKEMIA_CBF                  ROSS_ACUTE_MYELOID_LEUKEMIA_CBF 1.4381418
# WONG_ADULT_TISSUE_STEM_MODULE                                       WONG_ADULT_TISSUE_STEM_MODULE                    WONG_ADULT_TISSUE_STEM_MODULE 1.2835072
# KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN                         KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN             KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN 1.2680660
# ROSS_AML_WITH_AML1_ETO_FUSION                                       ROSS_AML_WITH_AML1_ETO_FUSION                    ROSS_AML_WITH_AML1_ETO_FUSION 1.2428077
# LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN 1.0218979
# JAATINEN_HEMATOPOIETIC_STEM_CELL_UP                           JAATINEN_HEMATOPOIETIC_STEM_CELL_UP              JAATINEN_HEMATOPOIETIC_STEM_CELL_UP 1.0206344
# VERHAAK_AML_WITH_NPM1_MUTATED_DN                                 VERHAAK_AML_WITH_NPM1_MUTATED_DN                 VERHAAK_AML_WITH_NPM1_MUTATED_DN 0.9585043
# YAGI_AML_WITH_T_8_21_TRANSLOCATION                             YAGI_AML_WITH_T_8_21_TRANSLOCATION               YAGI_AML_WITH_T_8_21_TRANSLOCATION 0.8897697
# MULLIGHAN_MLL_SIGNATURE_2_DN                                         MULLIGHAN_MLL_SIGNATURE_2_DN                     MULLIGHAN_MLL_SIGNATURE_2_DN 0.6958040
# TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_DN                   TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_DN          TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_DN 0.6792279
# pvalue  p.adjust
# REACTOME_INNATE_IMMUNE_SYSTEM                    0.01487103 0.1998002
# YAGI_AML_FAB_MARKERS                             0.03996004 0.1998002
# SMID_BREAST_CANCER_NORMAL_LIKE_UP                0.03996004 0.1998002
# CHEN_METABOLIC_SYNDROM_NETWORK                   0.09090909 0.3409091
# KRIGE_RESPONSE_TO_TOSEDOSTAT_6HR_DN              0.18681319 0.4770230
# ROSS_ACUTE_MYELOID_LEUKEMIA_CBF                  0.19080919 0.4770230
# WONG_ADULT_TISSUE_STEM_MODULE                    0.27972028 0.4995005
# KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN             0.29070929 0.4995005
# ROSS_AML_WITH_AML1_ETO_FUSION                    0.29970030 0.4995005
# LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN 0.47252747 0.6498047
# JAATINEN_HEMATOPOIETIC_STEM_CELL_UP              0.47652348 0.6498047
# VERHAAK_AML_WITH_NPM1_MUTATED_DN                 0.52347652 0.6543457
# YAGI_AML_WITH_T_8_21_TRANSLOCATION               0.56743257 0.6547299
# MULLIGHAN_MLL_SIGNATURE_2_DN                     0.72227772 0.7352647
# TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_DN          0.73526474 0.7352647



# GO
gse_go <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed = 123
)

# 查看显著通路
as.data.frame(gse_go)[, c("ID", "Description", "NES", "pvalue", "p.adjust")]
# ID                                                               Description      NES      pvalue   p.adjust
# GO:0098542 GO:0098542                                        defense response to other organism 2.264336 0.000451171 0.03724163
# GO:0045087 GO:0045087                                                    innate immune response 2.179977 0.001306956 0.03724163
# GO:0009607 GO:0009607                                               response to biotic stimulus 2.154016 0.001533479 0.03724163
# GO:0043207 GO:0043207                                      response to external biotic stimulus 2.154016 0.001533479 0.03724163
# GO:0044419 GO:0044419 biological process involved in interspecies interaction between organisms 2.154016 0.001533479 0.03724163
# GO:0051707 GO:0051707                                                response to other organism 2.154016 0.001533479 0.03724163
# GO:0140546 GO:0140546                                              defense response to symbiont 2.129749 0.001092333 0.03724163

# Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

gse_hallmark <- GSEA(
  geneList = geneList,
  TERM2GENE = msig_hallmark,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  scoreType = "pos",
  seed = 123
)
# 查看显著通路
as.data.frame(gse_hallmark)[, c("ID", "Description", "NES", "pvalue")]

# ReactomePA富集
library(ReactomePA)
gse_reactome <- gsePathway(
  geneList     = geneList,
  organism     = "human",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  scoreType    = "pos"
)
as.data.frame(gse_reactome)[, c("ID", "Description", "NES", "pvalue", "p.adjust")]

# MSigDB C6 癌基因特征
msig_c6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
gse_c6 <- GSEA(
  geneList  = geneList,
  TERM2GENE = msig_c6,
  scoreType = "pos",
  pvalueCutoff = 1
)
as.data.frame(gse_c6)[, c("ID", "Description", "NES", "pvalue", "p.adjust")]

# MSigDB C7 免疫特征
msig_c7 <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)
gse_c7 <- GSEA(
  geneList  = geneList,
  TERM2GENE = msig_c7,
  scoreType = "pos",
  pvalueCutoff = 1
)
as.data.frame(gse_c7)[, c("ID", "Description", "NES", "pvalue", "p.adjust")]




# ======================3.画图  ======================
# 你的GSEA结果对象（已存在）
# 这里选择哪种富集分析结果
em <- gse_c2

# ====================== 3.1. 单基因集GSEA图（第1个通路） ======================
title1 <- em$Description[1]
filename1 <- paste0("./1_GSEA/1.", clean_filename(title1), ".pdf")
pdf(filename1, width = 10, height = 6, pointsize = 8)
gseaplot2(em, 
          geneSetID = 1, 
          title = title1,
          color = "pink",
          base_size = 10,
          pvalue_table = TRUE,          # 【关键】显示NES、pvalue、p.adjust
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3)
dev.off()


# ====================== 3.2. 多基因集GSEA叠加图（1-8号通路） ======================
title3 <- "GSEA_Enrichment"
filename3 <- paste0("./1_GSEA/2.", clean_filename(title3), ".pdf")
pdf(filename3, width = 10, height = 7, pointsize = 8)
gseaplot2(em, 
          geneSetID = 1:5, 
          title = title3,
          color = c("pink", "lightblue", "#EFC000"),
          pvalue_table = TRUE,
          base_size = 10,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3)
dev.off()



# ====================== 3.3. 气泡图（Dotplot） ======================
title4 <- "Pathway_Enrichment"
filename4 <- paste0("./1_GSEA//3.", clean_filename(title4), ".pdf")
pdf(filename4, width = 6, height = 8, pointsize = 8)
dotplot(em, 
        showCategory = 15,
        color = "p.adjust",
        font.size = 8,
        title = title4) +
  scale_color_gradientn(colours = c("lightblue", "white", "pink"),
                        name = "Adjusted p-value") +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 7))
dev.off()

# ====================== 3.4. 脊线图（Ridgeplot） ======================
title5 <- "GSEA Enrichment Analysis"
filename5 <- paste0("./1_GSEA//4.2", clean_filename(title5), ".pdf")
pdf(filename5, width = 10, height = 9, pointsize = 8)
ridgeplot(em, 
          showCategory = 15,
          # fill = "NES",  # "pvalue", "p.adjust","NES"都能上色
          # fill = "pvalue",
          fill = "p.adjust",
          core_enrichment = FALSE) +
  scale_fill_gradientn(colours = c("lightblue", "#EFC000", "pink"),
                       name = "p.adjust") +
  labs(title = title5, x = "Rank in Ordered Gene List") + # ✅ 横坐标改成你要的
  xlim(0, 10) +  # 直接截断横坐标，只显示 0~10 的范围
  theme(axis.text.y = element_text(size = 7),
        legend.key.height = unit(1, "cm"))
dev.off()

save.image("./1_GSEA/20260403_gsea.Rdata")

